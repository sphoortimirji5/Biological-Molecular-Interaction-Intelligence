"""
TrainingManager — end-to-end training pipeline orchestrator.

Wires together:
  TrainingDataBuilder → XGBoostWrapper → ModelEvaluator → persistence
"""
import json
import io
import structlog
from typing import Dict, Any, Optional
from pathlib import Path
from datetime import datetime, timezone

from src.features.store import FeatureStore
from src.models.training_data import TrainingDataBuilder
from src.models.xgboost_model import XGBoostWrapper
from src.models.evaluator import ModelEvaluator

logger = structlog.get_logger(__name__)


class TrainingManager:
    """Orchestrate the full training lifecycle."""

    def __init__(
        self,
        feature_store: FeatureStore,
        model_params: Dict[str, Any] | None = None,
        model_version: str = "v1",
        compound_feature_version: str = "morgan_v1",
        protein_feature_version: str = "esm2_v1",
    ):
        self.feature_store = feature_store
        self.model_params = model_params
        self.model_version = model_version
        self.compound_feature_version = compound_feature_version
        self.protein_feature_version = protein_feature_version

    async def run(self) -> Dict[str, Any]:
        """
        Execute the full pipeline: build data → train → evaluate → persist.

        Returns:
            Dict with model artifact path and evaluation report.
        """
        logger.info("training_pipeline_started", version=self.model_version)

        # 1. Assemble training data
        builder = TrainingDataBuilder(feature_store=self.feature_store)
        splits = await builder.build(
            compound_feature_version=self.compound_feature_version,
            protein_feature_version=self.protein_feature_version,
        )

        if splits["train"].empty:
            raise ValueError("No training data assembled — check interactions and features.")

        # 2. Train XGBoost
        model = XGBoostWrapper(
            params=self.model_params,
            feature_version=self.model_version,
        )
        model.train(splits["train"])

        # 3. Evaluate on val and test
        evaluator = ModelEvaluator()
        reports = []
        for split_name in ["val", "test"]:
            if splits[split_name].empty:
                continue
            preds = model.predict(splits[split_name])
            report = evaluator.evaluate(
                y_true=preds["label"].values,
                y_prob=preds["probability"].values,
                split_name=split_name,
            )
            reports.append(report)

        comparison = evaluator.compare(reports) if reports else {}

        # 4. Persist model + evaluation report to S3
        model_s3_path = self._persist_model(model)
        report_s3_path = self._persist_report(comparison, reports)

        result = {
            "model_path": model_s3_path,
            "report_path": report_s3_path,
            "model_version": self.model_version,
            "training_metadata": splits["metadata"],
            "evaluation": comparison,
            "completed_at": datetime.now(timezone.utc).isoformat(),
        }

        logger.info(
            "training_pipeline_completed",
            model_path=model_s3_path,
            auroc=comparison.get("comparison", {}).get("test", {}).get("auroc"),
        )
        return result

    def _persist_model(self, model: XGBoostWrapper) -> str:
        """Save model to a temp dir, then upload to S3."""
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            local_path = Path(tmpdir) / "model"
            model.save(local_path)

            s3_prefix = f"models/xgboost/{self.model_version}"

            # Upload model pickle
            with open(local_path / "model.pkl", "rb") as f:
                self.feature_store.storage.put_object(
                    bucket=self.feature_store.bucket,
                    key=f"{s3_prefix}/model.pkl",
                    body=io.BytesIO(f.read()),
                    content_type="application/octet-stream",
                )

            # Upload signature
            with open(local_path / "signature.json", "rb") as f:
                self.feature_store.storage.put_object(
                    bucket=self.feature_store.bucket,
                    key=f"{s3_prefix}/signature.json",
                    body=io.BytesIO(f.read()),
                    content_type="application/json",
                )

        s3_uri = f"s3://{self.feature_store.bucket}/{s3_prefix}/model.pkl"
        logger.info("model_persisted", path=s3_uri)
        return s3_uri

    def _persist_report(
        self,
        comparison: Dict[str, Any],
        reports: list,
    ) -> str:
        """Save evaluation report to S3."""
        report_data = {
            "comparison": comparison,
            "detailed_reports": reports,
        }
        s3_key = f"models/xgboost/{self.model_version}/evaluation.json"
        body = json.dumps(report_data, indent=2, default=str)
        self.feature_store.storage.put_object(
            bucket=self.feature_store.bucket,
            key=s3_key,
            body=io.BytesIO(body.encode("utf-8")),
            content_type="application/json",
        )
        s3_uri = f"s3://{self.feature_store.bucket}/{s3_key}"
        logger.info("evaluation_report_persisted", path=s3_uri)
        return s3_uri
