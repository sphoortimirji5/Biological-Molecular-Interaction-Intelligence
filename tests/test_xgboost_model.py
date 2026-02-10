"""Unit tests for XGBoostWrapper — trains on synthetic data, no DB needed."""
import pytest
import pandas as pd
import numpy as np
from pathlib import Path

from src.models.xgboost_model import XGBoostWrapper


@pytest.fixture
def synthetic_train_data():
    """50 samples with 4 fp cols + 3 emb cols + label."""
    rng = np.random.RandomState(42)
    n = 50
    return pd.DataFrame({
        "label": rng.randint(0, 2, n),
        "fp_0": rng.randint(0, 2, n),
        "fp_1": rng.randint(0, 2, n),
        "fp_2": rng.randint(0, 2, n),
        "fp_3": rng.randint(0, 2, n),
        "emb_0": rng.randn(n),
        "emb_1": rng.randn(n),
        "emb_2": rng.randn(n),
    })


@pytest.fixture
def trained_model(synthetic_train_data):
    wrapper = XGBoostWrapper(feature_version="test_v1")
    wrapper.train(synthetic_train_data)
    return wrapper


def test_train_and_predict(trained_model, synthetic_train_data):
    """Train → predict produces probability and prediction columns."""
    result = trained_model.predict(synthetic_train_data)

    assert "prediction" in result.columns
    assert "probability" in result.columns
    assert len(result) == len(synthetic_train_data)

    # Predictions are binary
    assert set(result["prediction"].unique()).issubset({0, 1})

    # Probabilities are between 0 and 1
    assert result["probability"].min() >= 0.0
    assert result["probability"].max() <= 1.0


def test_predict_without_training_raises():
    """Calling predict before train should raise RuntimeError."""
    wrapper = XGBoostWrapper()
    df = pd.DataFrame({"fp_0": [1], "emb_0": [0.5]})

    with pytest.raises(RuntimeError, match="not trained"):
        wrapper.predict(df)


def test_save_and_load(trained_model, synthetic_train_data, tmp_path):
    """Save → load → predict produces same results."""
    model_dir = tmp_path / "xgboost_test"

    # Save
    trained_model.save(model_dir)
    assert (model_dir / "model.pkl").exists()
    assert (model_dir / "signature.json").exists()

    # Load
    loaded = XGBoostWrapper()
    loaded.load(model_dir)

    # Predictions should match
    original_preds = trained_model.predict(synthetic_train_data)
    loaded_preds = loaded.predict(synthetic_train_data)

    np.testing.assert_allclose(
        original_preds["probability"].values,
        loaded_preds["probability"].values,
        atol=1e-6,
    )


def test_model_signature(trained_model):
    """Signature must contain model type, hyperparameters, and feature info."""
    sig = trained_model.model_signature()

    assert sig["model_type"] == "xgboost"
    assert sig["feature_version"] == "test_v1"
    assert "hyperparameters" in sig
    assert sig["hyperparameters"]["objective"] == "binary:logistic"
    assert "saved_at" in sig
    assert sig["feature_columns"]["compound"] == "fp_0..fp_2047"


def test_custom_params():
    """Custom hyperparameters should override defaults."""
    wrapper = XGBoostWrapper(params={"max_depth": 3, "n_estimators": 50})

    assert wrapper.params["max_depth"] == 3
    assert wrapper.params["n_estimators"] == 50
    assert wrapper.params["objective"] == "binary:logistic"  # default preserved
