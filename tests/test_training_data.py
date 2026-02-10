"""Unit tests for TrainingDataBuilder — fully mocked, no DB or S3 needed."""
import pytest
import pandas as pd
import numpy as np
from unittest.mock import MagicMock, AsyncMock, patch

from src.models.training_data import TrainingDataBuilder


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def mock_feature_store():
    store = MagicMock()

    # Synthetic compound features: 4 compounds × 4 fingerprint cols
    store.load.side_effect = lambda entity_type, version: {
        ("compound", "morgan_v1"): pd.DataFrame({
            "external_id": ["C1", "C2", "C3", "C4"],
            "source": ["test"] * 4,
            "fp_0": [1, 0, 1, 0],
            "fp_1": [0, 1, 0, 1],
            "fp_2": [1, 1, 0, 0],
            "fp_3": [0, 0, 1, 1],
        }),
        ("protein", "esm2_v1"): pd.DataFrame({
            "external_id": ["P1", "P2", "P3", "P4"],
            "source": ["test"] * 4,
            "emb_0": [0.1, 0.2, 0.3, 0.4],
            "emb_1": [0.5, 0.6, 0.7, 0.8],
            "emb_2": [0.9, 1.0, 1.1, 1.2],
        }),
    }[(entity_type, version)]

    return store


@pytest.fixture
def interactions_df():
    """20 interactions with balanced labels (10 positive, 10 negative)."""
    rows = []
    compounds = [("C1", "test"), ("C2", "test"), ("C3", "test"), ("C4", "test")]
    proteins = [("P1", "test"), ("P2", "test"), ("P3", "test"), ("P4", "test")]

    idx = 0
    for c_id, c_src in compounds:
        for p_id, p_src in proteins:
            rows.append({
                "label": 1 if idx < 8 else 0,
                "compound_external_id": c_id,
                "compound_source": c_src,
                "protein_external_id": p_id,
                "protein_source": p_src,
            })
            idx += 1

    # 16 rows total: 8 positive, 8 negative
    return pd.DataFrame(rows)


@pytest.fixture
def builder(mock_feature_store):
    return TrainingDataBuilder(
        feature_store=mock_feature_store,
        test_size=0.25,
        val_size=0.25,
        random_state=42,
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.asyncio
async def test_build_produces_correct_splits(builder, interactions_df):
    """Build must return train, val, test DataFrames with correct proportions."""
    with patch.object(builder, "_query_interactions", new_callable=AsyncMock, return_value=interactions_df):
        result = await builder.build()

    total = len(result["train"]) + len(result["val"]) + len(result["test"])
    assert total == 16
    assert len(result["test"]) == 4   # 25% of 16
    assert len(result["val"]) == 4    # 25% of remaining 12 ≈ 3-4
    assert len(result["train"]) == 8  # remainder


@pytest.mark.asyncio
async def test_feature_columns_are_concatenated(builder, interactions_df):
    """Output must have both fp_ and emb_ columns."""
    with patch.object(builder, "_query_interactions", new_callable=AsyncMock, return_value=interactions_df):
        result = await builder.build()

    train = result["train"]
    fp_cols = [c for c in train.columns if c.startswith("fp_")]
    emb_cols = [c for c in train.columns if c.startswith("emb_")]

    assert len(fp_cols) == 4   # fp_0..fp_3
    assert len(emb_cols) == 3  # emb_0..emb_2
    assert "label" in train.columns
    assert train.shape[1] == 8  # label + 4 fp + 3 emb


@pytest.mark.asyncio
async def test_no_data_leakage(builder, interactions_df):
    """Train, val, and test must have no overlapping indices."""
    with patch.object(builder, "_query_interactions", new_callable=AsyncMock, return_value=interactions_df):
        result = await builder.build()

    train_idx = set(result["train"].index)
    val_idx = set(result["val"].index)
    test_idx = set(result["test"].index)

    # Indices are reset, so check by row content instead
    train_rows = set(result["train"].apply(tuple, axis=1))
    val_rows = set(result["val"].apply(tuple, axis=1))
    test_rows = set(result["test"].apply(tuple, axis=1))

    assert train_rows.isdisjoint(val_rows), "Train/val overlap"
    assert train_rows.isdisjoint(test_rows), "Train/test overlap"
    assert val_rows.isdisjoint(test_rows), "Val/test overlap"


@pytest.mark.asyncio
async def test_stratification_preserves_label_ratio(builder, interactions_df):
    """Each split should roughly preserve the label distribution."""
    with patch.object(builder, "_query_interactions", new_callable=AsyncMock, return_value=interactions_df):
        result = await builder.build()

    for split_name in ["train", "val", "test"]:
        split = result[split_name]
        if len(split) > 0:
            ratio = split["label"].mean()
            # Original ratio is 8/16 = 0.5, each split should be close
            assert 0.25 <= ratio <= 0.75, f"{split_name} label ratio {ratio} too skewed"


@pytest.mark.asyncio
async def test_metadata_is_populated(builder, interactions_df):
    """Metadata must contain sample counts and feature info."""
    with patch.object(builder, "_query_interactions", new_callable=AsyncMock, return_value=interactions_df):
        result = await builder.build()

    meta = result["metadata"]
    assert meta["total_samples"] == 16
    assert meta["train_samples"] + meta["val_samples"] + meta["test_samples"] == 16
    assert meta["feature_dim"] == 7  # 4 fp + 3 emb
    assert meta["compound_feature_version"] == "morgan_v1"
    assert meta["protein_feature_version"] == "esm2_v1"


@pytest.mark.asyncio
async def test_empty_interactions_returns_empty(builder):
    """No interactions → empty DataFrames, no crash."""
    empty_df = pd.DataFrame()

    with patch.object(builder, "_query_interactions", new_callable=AsyncMock, return_value=empty_df):
        result = await builder.build()

    assert result["train"].empty
    assert result["val"].empty
    assert result["test"].empty
