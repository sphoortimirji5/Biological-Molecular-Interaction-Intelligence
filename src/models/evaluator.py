"""
ModelEvaluator — computes classification and ranking metrics for DTI models.

Produces a structured evaluation report (JSON-serializable dict) with:
- AUROC, AUPRC (threshold-agnostic)
- Precision@k, Recall@k (top-k ranking quality)
- Accuracy, F1 (threshold-based at p=0.5)
"""
import structlog
import numpy as np
from typing import Dict, Any, List
from datetime import datetime, timezone

from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    precision_score,
    recall_score,
    f1_score,
    accuracy_score,
)

logger = structlog.get_logger(__name__)


class ModelEvaluator:
    """Evaluate binary classification and ranking performance."""

    def __init__(self, k_values: List[int] | None = None):
        self.k_values = k_values or [10, 50, 100]

    def evaluate(
        self,
        y_true: np.ndarray,
        y_prob: np.ndarray,
        split_name: str = "test",
    ) -> Dict[str, Any]:
        """
        Compute all metrics for a single evaluation split.

        Args:
            y_true: Ground-truth binary labels (0/1).
            y_prob: Predicted probabilities (0.0-1.0).
            split_name: Label for this evaluation (e.g. 'test', 'val').

        Returns:
            Dict with per-metric scores and metadata.
        """
        y_true = np.asarray(y_true)
        y_prob = np.asarray(y_prob)
        y_pred = (y_prob >= 0.5).astype(int)

        report = {
            "split": split_name,
            "n_samples": len(y_true),
            "n_positive": int(y_true.sum()),
            "n_negative": int((1 - y_true).sum()),
            "threshold_free": {
                "auroc": self._safe_auroc(y_true, y_prob),
                "auprc": self._safe_auprc(y_true, y_prob),
            },
            "threshold_based": {
                "accuracy": float(accuracy_score(y_true, y_pred)),
                "precision": float(precision_score(y_true, y_pred, zero_division=0)),
                "recall": float(recall_score(y_true, y_pred, zero_division=0)),
                "f1": float(f1_score(y_true, y_pred, zero_division=0)),
            },
            "ranking": self._ranking_metrics(y_true, y_prob),
            "evaluated_at": datetime.now(timezone.utc).isoformat(),
        }

        logger.info(
            "model_evaluated",
            split=split_name,
            auroc=report["threshold_free"]["auroc"],
            auprc=report["threshold_free"]["auprc"],
            f1=report["threshold_based"]["f1"],
        )
        return report

    def compare(self, reports: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Compare multiple evaluation reports side-by-side.

        Returns a summary dict keyed by split name with key metrics.
        """
        summary = {}
        for r in reports:
            summary[r["split"]] = {
                "auroc": r["threshold_free"]["auroc"],
                "auprc": r["threshold_free"]["auprc"],
                "f1": r["threshold_based"]["f1"],
                "accuracy": r["threshold_based"]["accuracy"],
                "n_samples": r["n_samples"],
            }
        return {
            "comparison": summary,
            "compared_at": datetime.now(timezone.utc).isoformat(),
        }

    def _ranking_metrics(self, y_true: np.ndarray, y_prob: np.ndarray) -> Dict[str, Any]:
        """Compute precision@k and recall@k for configured k values."""
        # Sort by predicted probability (descending)
        ranked_indices = np.argsort(-y_prob)
        ranked_labels = y_true[ranked_indices]
        total_positives = int(y_true.sum())

        metrics = {}
        for k in self.k_values:
            if k > len(y_true):
                continue
            top_k = ranked_labels[:k]
            tp_at_k = int(top_k.sum())
            metrics[f"precision@{k}"] = round(tp_at_k / k, 4)
            metrics[f"recall@{k}"] = round(tp_at_k / max(total_positives, 1), 4)

        return metrics

    @staticmethod
    def _safe_auroc(y_true: np.ndarray, y_prob: np.ndarray) -> float | None:
        """AUROC — returns None if only one class present."""
        if len(np.unique(y_true)) < 2:
            return None
        return float(roc_auc_score(y_true, y_prob))

    @staticmethod
    def _safe_auprc(y_true: np.ndarray, y_prob: np.ndarray) -> float | None:
        """AUPRC — returns None if only one class present."""
        if len(np.unique(y_true)) < 2:
            return None
        return float(average_precision_score(y_true, y_prob))
