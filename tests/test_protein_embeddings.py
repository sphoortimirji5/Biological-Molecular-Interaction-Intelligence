"""Unit tests for ProteinEmbeddingGenerator — real ESM-2 model, tiny sequences."""
import pytest
import pandas as pd
import numpy as np

from src.features.protein_embeddings import ProteinEmbeddingGenerator


# Short real protein sequences (fast inference)
INSULIN_FRAGMENT = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"
UBIQUITIN_FRAGMENT = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGR"


@pytest.fixture(scope="module")
def generator():
    """Module-scoped: loads model once for all tests."""
    gen = ProteinEmbeddingGenerator(
        model_name="facebook/esm2_t6_8M_UR50D",
        batch_size=2,
        max_length=128,
    )
    gen.fit(pd.DataFrame())  # Pre-load model
    return gen


@pytest.fixture
def sample_df():
    return pd.DataFrame({
        "external_id": ["P01308", "P0CG48"],
        "source": ["uniprot", "uniprot"],
        "sequence": [INSULIN_FRAGMENT, UBIQUITIN_FRAGMENT],
    })


def test_transform_produces_correct_shape(generator, sample_df):
    """2 sequences → 2 rows × (2 id cols + 320 embedding cols)."""
    result = generator.transform(sample_df)

    assert result.shape[0] == 2
    assert result.shape[1] == 322  # external_id + source + 320 dims
    assert list(result.columns[:2]) == ["external_id", "source"]
    assert result.columns[2] == "emb_0"
    assert result.columns[-1] == "emb_319"


def test_embeddings_are_deterministic(generator, sample_df):
    """Same sequences must produce identical embeddings."""
    result1 = generator.transform(sample_df)
    result2 = generator.transform(sample_df)

    emb_cols = [c for c in result1.columns if c.startswith("emb_")]
    np.testing.assert_allclose(
        result1[emb_cols].values,
        result2[emb_cols].values,
        atol=1e-6,
    )


def test_embeddings_are_float(generator, sample_df):
    """Embedding values must be floats, not binary."""
    result = generator.transform(sample_df)
    emb_cols = [c for c in result.columns if c.startswith("emb_")]

    # Should have non-integer values (unlike Morgan fingerprints)
    values = result[emb_cols].values
    assert not np.all(np.isin(values, [0, 1])), "Embeddings should not be binary"


def test_different_proteins_produce_different_embeddings(generator, sample_df):
    """Insulin and Ubiquitin must have different embeddings."""
    result = generator.transform(sample_df)
    emb1 = result.iloc[0, 2:].values.astype(float)
    emb2 = result.iloc[1, 2:].values.astype(float)

    assert not np.allclose(emb1, emb2, atol=1e-4)


def test_feature_manifest(generator):
    """Manifest must pin model name, dimension, and pooling strategy."""
    manifest = generator.feature_manifest()

    assert manifest["type"] == "esm2"
    assert manifest["model"] == "facebook/esm2_t6_8M_UR50D"
    assert manifest["dim"] == 320
    assert manifest["pooling"] == "mean"
    assert "generated_at" in manifest
