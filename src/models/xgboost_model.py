"""
XGBoostWrapper — binary classification model for DTI prediction.

Implements the ModelWrapper interface for training, prediction,
serialization, and reproducibility metadata.
"""
import json
import pickle
import structlog
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Any
from datetime import datetime, timezone

import xgboost as xgb

from src.models.interfaces import ModelWrapper

logger = structlog.get_logger(__name__)

# Default hyperparameters for the baseline model
DEFAULT_PARAMS = {
    "objective": "binary:logistic",
    "eval_metric": "logloss",
    "max_depth": 6,
    "learning_rate": 0.1,
    "n_estimators": 100,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "min_child_weight": 1,
    "random_state": 42,
    "n_jobs": -1,
    "verbosity": 0,
}


class XGBoostWrapper(ModelWrapper):
    """XGBoost binary classifier for drug-target interaction prediction."""

    def __init__(
        self,
        params: Dict[str, Any] | None = None,
        feature_version: str = "v1",
    ):
        self.params = {**DEFAULT_PARAMS, **(params or {})}
        self.feature_version = feature_version
        self._model: xgb.XGBClassifier | None = None

    def train(self, df: pd.DataFrame) -> None:
        """
        Train XGBoost on a DataFrame with 'label' column + feature columns.

        Feature columns are identified as columns starting with 'fp_' or 'emb_'.
        """
        feature_cols = self._get_feature_cols(df)
        X = df[feature_cols].values
        y = df["label"].values

        self._model = xgb.XGBClassifier(**self.params)
        self._model.fit(X, y)

        logger.info(
            "xgboost_trained",
            samples=len(df),
            features=len(feature_cols),
            params=self.params,
        )

    def predict(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Run inference. Returns DataFrame with 'prediction' and 'probability' columns.
        """
        if self._model is None:
            raise RuntimeError("Model not trained. Call train() or load() first.")

        feature_cols = self._get_feature_cols(df)
        X = df[feature_cols].values

        probabilities = self._model.predict_proba(X)[:, 1]
        predictions = (probabilities >= 0.5).astype(int)

        result = df.copy()
        result["prediction"] = predictions
        result["probability"] = probabilities
        return result

    def save(self, path: Path) -> None:
        """Serialize model + metadata to a directory."""
        if self._model is None:
            raise RuntimeError("No model to save. Call train() first.")

        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)

        # Save model
        model_path = path / "model.pkl"
        with open(model_path, "wb") as f:
            pickle.dump(self._model, f)

        # Save signature alongside model
        sig_path = path / "signature.json"
        with open(sig_path, "w") as f:
            json.dump(self.model_signature(), f, indent=2)

        logger.info("xgboost_saved", path=str(path))

    def load(self, path: Path) -> "XGBoostWrapper":
        """Deserialize model from a directory."""
        path = Path(path)
        model_path = path / "model.pkl"

        with open(model_path, "rb") as f:
            self._model = pickle.load(f)

        # Load signature to restore metadata
        sig_path = path / "signature.json"
        if sig_path.exists():
            with open(sig_path) as f:
                sig = json.load(f)
                self.params = sig.get("hyperparameters", self.params)
                self.feature_version = sig.get("feature_version", self.feature_version)

        logger.info("xgboost_loaded", path=str(path))
        return self

    def model_signature(self) -> Dict[str, Any]:
        """Reproducibility metadata — pins hyperparameters and feature requirements."""
        return {
            "model_type": "xgboost",
            "feature_version": self.feature_version,
            "hyperparameters": self.params,
            "feature_columns": {
                "compound": "fp_0..fp_2047",
                "protein": "emb_0..emb_319",
            },
            "saved_at": datetime.now(timezone.utc).isoformat(),
        }

    @staticmethod
    def _get_feature_cols(df: pd.DataFrame) -> list:
        """Extract feature column names (fp_* and emb_*)."""
        return [c for c in df.columns if c.startswith(("fp_", "emb_"))]
