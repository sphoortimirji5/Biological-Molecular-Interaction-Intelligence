"""Unit tests for TrainingManager — mocked pipeline, no DB or S3."""
import pytest
import pandas as pd
import numpy as np
from unittest.mock import MagicMock, AsyncMock, patch

from src.models.trainer import TrainingManager


@pytest.fixture
def mock_feature_store():
    store = MagicMock()

    # Simulate feature loading for TrainingDataBuilder
    def load_side_effect(entity_type, version):
        if entity_type == "compound":
            return pd.DataFrame({
                "external_id": [f"C{i}" for i in range(10)],
                "source": ["test"] * 10,
                "fp_0": np.random.randint(0, 2, 10),
                "fp_1": np.random.randint(0, 2, 10),
            })
        else:
            return pd.DataFrame({
                "external_id": [f"P{i}" for i in range(10)],
                "source": ["test"] * 10,
                "emb_0": np.random.randn(10),
                "emb_1": np.random.randn(10),
            })

    store.load.side_effect = load_side_effect
    store.bucket = "test-bucket"
    store.storage = MagicMock()
    store.storage.put_object = MagicMock()
    return store


@pytest.fixture
def interactions_df():
    """40 interactions: 10 compounds × 4 proteins."""
    rows = []
    for i in range(10):
        for j in range(4):
            rows.append({
                "label": 1 if (i + j) % 2 == 0 else 0,
                "compound_external_id": f"C{i}",
                "compound_source": "test",
                "protein_external_id": f"P{j}",
                "protein_source": "test",
            })
    return pd.DataFrame(rows)


@pytest.mark.asyncio
async def test_run_produces_result(mock_feature_store, interactions_df):
    """Full pipeline run returns model path, report, and evaluation."""
    manager = TrainingManager(
        feature_store=mock_feature_store,
        model_version="test_v1",
    )

    with patch(
        "src.models.training_data.TrainingDataBuilder._query_interactions",
        new_callable=AsyncMock,
        return_value=interactions_df,
    ):
        result = await manager.run()

    assert "model_path" in result
    assert "report_path" in result
    assert "evaluation" in result
    assert result["model_version"] == "test_v1"
    assert result["training_metadata"]["total_samples"] > 0


@pytest.mark.asyncio
async def test_model_and_report_persisted(mock_feature_store, interactions_df):
    """Model pickle and evaluation JSON must both be uploaded to S3."""
    manager = TrainingManager(
        feature_store=mock_feature_store,
        model_version="test_v1",
    )

    with patch(
        "src.models.training_data.TrainingDataBuilder._query_interactions",
        new_callable=AsyncMock,
        return_value=interactions_df,
    ):
        await manager.run()

    # Check that put_object was called for model.pkl, signature.json, and evaluation.json
    put_calls = mock_feature_store.storage.put_object.call_args_list
    keys = [call.kwargs.get("key", call.args[0] if call.args else "") for call in put_calls]

    # Flatten: check keyword args
    uploaded_keys = []
    for call in put_calls:
        _, kwargs = call
        uploaded_keys.append(kwargs.get("key", ""))

    assert any("model.pkl" in k for k in uploaded_keys), f"model.pkl not found in {uploaded_keys}"
    assert any("signature.json" in k for k in uploaded_keys), f"signature.json not found in {uploaded_keys}"
    assert any("evaluation.json" in k for k in uploaded_keys), f"evaluation.json not found in {uploaded_keys}"


@pytest.mark.asyncio
async def test_empty_data_raises(mock_feature_store):
    """No training data should raise ValueError."""
    manager = TrainingManager(feature_store=mock_feature_store)

    with patch(
        "src.models.training_data.TrainingDataBuilder._query_interactions",
        new_callable=AsyncMock,
        return_value=pd.DataFrame(),
    ):
        with pytest.raises(ValueError, match="No training data"):
            await manager.run()


@pytest.mark.asyncio
async def test_evaluation_has_metrics(mock_feature_store, interactions_df):
    """Evaluation report should contain AUROC and F1 for test split."""
    manager = TrainingManager(
        feature_store=mock_feature_store,
        model_version="test_v1",
    )

    with patch(
        "src.models.training_data.TrainingDataBuilder._query_interactions",
        new_callable=AsyncMock,
        return_value=interactions_df,
    ):
        result = await manager.run()

    comparison = result["evaluation"].get("comparison", {})
    assert "test" in comparison or "val" in comparison
    for split_data in comparison.values():
        assert "auroc" in split_data
        assert "f1" in split_data
