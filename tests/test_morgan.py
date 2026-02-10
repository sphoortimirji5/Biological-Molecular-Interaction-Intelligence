"""Unit tests for MorganFingerprintGenerator — deterministic fingerprint verification."""
import pytest
import pandas as pd
import numpy as np
from unittest.mock import patch

from src.features.morgan import MorganFingerprintGenerator


# Well-known SMILES for reproducibility
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"
ETHANOL = "CCO"


@pytest.fixture
def generator():
    return MorganFingerprintGenerator(radius=2, n_bits=2048)


@pytest.fixture
def sample_df():
    return pd.DataFrame({
        "external_id": ["ASPIRIN", "CAFFEINE", "ETHANOL"],
        "source": ["chembl", "chembl", "chembl"],
        "smiles": [ASPIRIN, CAFFEINE, ETHANOL],
    })


def test_transform_produces_correct_shape(generator, sample_df):
    """3 valid SMILES → 3 rows × (2 id cols + 2048 fp cols)."""
    result = generator.transform(sample_df)

    assert result.shape == (3, 2050)  # external_id + source + 2048 bits
    assert list(result.columns[:2]) == ["external_id", "source"]
    assert list(result.columns[2:5]) == ["fp_0", "fp_1", "fp_2"]


def test_fingerprints_are_deterministic(generator, sample_df):
    """Same SMILES must produce identical fingerprints across runs."""
    result1 = generator.transform(sample_df)
    result2 = generator.transform(sample_df)

    pd.testing.assert_frame_equal(result1, result2)


def test_fingerprints_are_binary(generator, sample_df):
    """All fingerprint values must be 0 or 1."""
    result = generator.transform(sample_df)
    fp_cols = [c for c in result.columns if c.startswith("fp_")]

    for col in fp_cols:
        assert result[col].isin([0, 1]).all(), f"Non-binary value in {col}"


def test_different_molecules_produce_different_fingerprints(generator, sample_df):
    """Aspirin and Caffeine must have different fingerprints."""
    result = generator.transform(sample_df)
    aspirin_fp = result.iloc[0, 2:].values
    caffeine_fp = result.iloc[1, 2:].values

    assert not np.array_equal(aspirin_fp, caffeine_fp)


def test_invalid_smiles_are_skipped(generator):
    """Invalid SMILES are skipped, logged, and counted."""
    df = pd.DataFrame({
        "external_id": ["GOOD", "BAD", "ALSO_GOOD"],
        "source": ["test", "test", "test"],
        "smiles": [ASPIRIN, "NOT_A_REAL_SMILES_XYZ", ETHANOL],
    })

    result = generator.transform(df)

    assert len(result) == 2  # Only valid ones
    assert generator.invalid_count == 1
    assert list(result["external_id"]) == ["GOOD", "ALSO_GOOD"]


def test_feature_manifest(generator):
    """Manifest must pin radius, n_bits, and RDKit version."""
    manifest = generator.feature_manifest()

    assert manifest["type"] == "morgan"
    assert manifest["radius"] == 2
    assert manifest["n_bits"] == 2048
    assert "rdkit_version" in manifest
    assert "generated_at" in manifest
