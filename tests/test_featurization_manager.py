"""Unit tests for FeaturizationManager — fully mocked, no DB or S3 needed."""
import pytest
import pandas as pd
from unittest.mock import MagicMock, AsyncMock, patch

from src.features.manager import FeaturizationManager


@pytest.fixture
def mock_feature_store():
    store = MagicMock()
    store.save.return_value = "s3://features/compound/morgan_v1/features.parquet"
    return store


@pytest.fixture
def mock_compound_generator():
    gen = MagicMock()
    gen.transform.return_value = pd.DataFrame({
        "external_id": ["ASPIRIN"],
        "source": ["chembl"],
        "fp_0": [1], "fp_1": [0],
    })
    gen.feature_manifest.return_value = {"type": "morgan", "radius": 2}
    return gen


@pytest.fixture
def mock_protein_generator():
    gen = MagicMock()
    gen.transform.return_value = pd.DataFrame({
        "external_id": ["P12345"],
        "source": ["uniprot"],
        "emb_0": [0.1], "emb_1": [0.2],
    })
    gen.feature_manifest.return_value = {"type": "esm2", "model": "esm2_t6_8M"}
    return gen


@pytest.fixture
def manager(mock_feature_store, mock_compound_generator, mock_protein_generator):
    return FeaturizationManager(
        feature_store=mock_feature_store,
        compound_generator=mock_compound_generator,
        protein_generator=mock_protein_generator,
    )


@pytest.mark.asyncio
async def test_featurize_compounds(manager, mock_feature_store, mock_compound_generator):
    """Compounds: query → transform → save to Feature Store."""
    compound_df = pd.DataFrame({
        "external_id": ["ASPIRIN"],
        "source": ["chembl"],
        "smiles": ["CC(=O)Oc1ccccc1C(=O)O"],
    })

    with patch.object(manager, "_query_compounds", new_callable=AsyncMock, return_value=compound_df):
        result = await manager.featurize_compounds("v1")

    mock_compound_generator.transform.assert_called_once()
    mock_feature_store.save.assert_called_once()
    assert result["count"] == 1
    assert result["uri"] is not None


@pytest.mark.asyncio
async def test_featurize_proteins(manager, mock_feature_store, mock_protein_generator):
    """Proteins: query → transform → save to Feature Store."""
    protein_df = pd.DataFrame({
        "external_id": ["P12345"],
        "source": ["uniprot"],
        "sequence": ["MALWMRLLPLL"],
    })

    with patch.object(manager, "_query_proteins", new_callable=AsyncMock, return_value=protein_df):
        result = await manager.featurize_proteins("v1")

    mock_protein_generator.transform.assert_called_once()
    mock_feature_store.save.assert_called_once()
    assert result["count"] == 1
    assert result["uri"] is not None


@pytest.mark.asyncio
async def test_full_run(manager, mock_feature_store):
    """Full run featurizes both compounds and proteins."""
    compound_df = pd.DataFrame({
        "external_id": ["ASPIRIN"],
        "source": ["chembl"],
        "smiles": ["CC(=O)Oc1ccccc1C(=O)O"],
    })
    protein_df = pd.DataFrame({
        "external_id": ["P12345"],
        "source": ["uniprot"],
        "sequence": ["MALWMRLLPLL"],
    })

    with patch.object(manager, "_query_compounds", new_callable=AsyncMock, return_value=compound_df), \
         patch.object(manager, "_query_proteins", new_callable=AsyncMock, return_value=protein_df):
        summary = await manager.run("v1")

    assert summary["version"] == "v1"
    assert summary["compounds"]["count"] == 1
    assert summary["proteins"]["count"] == 1
    assert summary["elapsed_seconds"] >= 0
    assert mock_feature_store.save.call_count == 2  # compounds + proteins


@pytest.mark.asyncio
async def test_empty_compounds_returns_zero(manager, mock_feature_store, mock_compound_generator):
    """No compounds in DB → no transform, no save."""
    empty_df = pd.DataFrame(columns=["external_id", "source", "smiles"])

    with patch.object(manager, "_query_compounds", new_callable=AsyncMock, return_value=empty_df):
        result = await manager.featurize_compounds("v1")

    mock_compound_generator.transform.assert_not_called()
    mock_feature_store.save.assert_not_called()
    assert result["count"] == 0
    assert result["uri"] is None


@pytest.mark.asyncio
async def test_empty_proteins_returns_zero(manager, mock_feature_store, mock_protein_generator):
    """No proteins in DB → no transform, no save."""
    empty_df = pd.DataFrame(columns=["external_id", "source", "sequence"])

    with patch.object(manager, "_query_proteins", new_callable=AsyncMock, return_value=empty_df):
        result = await manager.featurize_proteins("v1")

    mock_protein_generator.transform.assert_not_called()
    mock_feature_store.save.assert_not_called()
    assert result["count"] == 0
    assert result["uri"] is None
