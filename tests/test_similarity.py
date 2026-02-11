"""
Unit tests for SimilarityService — cosine similarity for proteins and compounds.
"""
import pytest
import numpy as np
import pandas as pd
from unittest.mock import MagicMock, patch

from src.inference.similarity import (
    SimilarityService,
    _cosine_similarity,
)


# ---------------------------------------------------------------------------
# Low-level similarity helpers
# ---------------------------------------------------------------------------


class TestCosineSimilarity:
    """Verify cosine similarity edge cases."""

    def test_identical_vectors(self):
        q = np.array([1.0, 0.0, 1.0])
        c = np.array([[1.0, 0.0, 1.0]])
        scores = _cosine_similarity(q, c)
        assert scores[0] == pytest.approx(1.0, abs=1e-6)

    def test_orthogonal_vectors(self):
        q = np.array([1.0, 0.0])
        c = np.array([[0.0, 1.0]])
        scores = _cosine_similarity(q, c)
        assert scores[0] == pytest.approx(0.0, abs=1e-6)

    def test_zero_query(self):
        q = np.zeros(3)
        c = np.array([[1.0, 2.0, 3.0]])
        scores = _cosine_similarity(q, c)
        assert scores[0] == 0.0

    def test_multiple_candidates(self):
        q = np.array([1.0, 0.0])
        c = np.array([[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]])
        scores = _cosine_similarity(q, c)
        assert len(scores) == 3
        assert scores[0] > scores[2] > scores[1]



# ---------------------------------------------------------------------------
# SimilarityService — compound similarity (no heavy model dependency)
# ---------------------------------------------------------------------------


class TestSimilarityServiceCompounds:
    """Test compound similarity with mocked FeatureStore."""

    def _make_service(self, compound_df):
        store = MagicMock()
        store.load.return_value = compound_df
        return SimilarityService(feature_store=store)

    def test_similar_compounds_returns_ranked(self):
        """Mocked fingerprints: Aspirin should rank highest against itself."""
        compound_df = pd.DataFrame([
            {"external_id": "C1", "source": "chembl", "fp_0": 1, "fp_1": 1, "fp_2": 0, "fp_3": 0},
            {"external_id": "C2", "source": "chembl", "fp_0": 0, "fp_1": 0, "fp_2": 1, "fp_3": 1},
        ])
        service = self._make_service(compound_df)

        # Patch the fingerprint generator to return a known vector
        mock_fp_df = pd.DataFrame([
            {"external_id": "query", "source": "similarity_search", "fp_0": 1, "fp_1": 1, "fp_2": 0, "fp_3": 0},
        ])
        service._fingerprint_gen = MagicMock()
        service._fingerprint_gen.transform.return_value = mock_fp_df

        results = service.find_similar_compounds(smiles="CC(=O)Oc1ccccc1C(=O)O", top_k=2)

        assert len(results) == 2
        assert results[0]["external_id"] == "C1"
        assert results[0]["score"] == pytest.approx(1.0, abs=1e-6)
        assert results[0]["rank"] == 1
        assert results[1]["external_id"] == "C2"
        assert results[1]["score"] == pytest.approx(0.0, abs=1e-6)

    def test_similar_compounds_invalid_smiles(self):
        compound_df = pd.DataFrame([
            {"external_id": "C1", "source": "chembl", "fp_0": 1, "fp_1": 0},
        ])
        service = self._make_service(compound_df)

        with pytest.raises(ValueError, match="Invalid SMILES"):
            service.find_similar_compounds(smiles="NOT_A_MOLECULE", top_k=5)

    def test_similar_compounds_empty_store(self):
        service = self._make_service(pd.DataFrame())

        with pytest.raises(ValueError, match="No compound fingerprints"):
            service.find_similar_compounds(smiles="CC(=O)Oc1ccccc1C(=O)O", top_k=5)

    def test_scores_are_descending(self):
        compound_df = pd.DataFrame([
            {"external_id": f"C{i}", "source": "chembl", "fp_0": i % 2, "fp_1": (i + 1) % 2}
            for i in range(5)
        ])
        service = self._make_service(compound_df)

        mock_fp_df = pd.DataFrame([
            {"external_id": "query", "source": "similarity_search", "fp_0": 1, "fp_1": 0},
        ])
        service._fingerprint_gen = MagicMock()
        service._fingerprint_gen.transform.return_value = mock_fp_df

        results = service.find_similar_compounds(smiles="C", top_k=5)

        scores = [r["score"] for r in results]
        assert scores == sorted(scores, reverse=True)


# ---------------------------------------------------------------------------
# SimilarityService — protein similarity (mock ESM-2)
# ---------------------------------------------------------------------------


class TestSimilarityServiceProteins:
    """Test protein similarity with mocked FeatureStore and ESM-2."""

    def test_similar_proteins_returns_ranked(self):
        protein_df = pd.DataFrame([
            {"external_id": "P1", "source": "uniprot", "emb_0": 1.0, "emb_1": 0.0},
            {"external_id": "P2", "source": "uniprot", "emb_0": 0.0, "emb_1": 1.0},
        ])
        store = MagicMock()
        store.load.return_value = protein_df
        service = SimilarityService(feature_store=store)

        query_emb_df = pd.DataFrame([
            {"external_id": "query", "source": "similarity_search", "emb_0": 1.0, "emb_1": 0.0},
        ])

        with patch("src.features.protein_embeddings.ProteinEmbeddingGenerator") as MockEmbedder:
            mock_instance = MagicMock()
            mock_instance.transform.return_value = query_emb_df
            MockEmbedder.return_value = mock_instance

            results = service.find_similar_proteins(sequence="MKTAYIAKQRQISFVKSH", top_k=2)

        assert len(results) == 2
        assert results[0]["external_id"] == "P1"
        assert results[0]["score"] == pytest.approx(1.0, abs=1e-6)
        assert results[1]["external_id"] == "P2"
        assert results[1]["score"] == pytest.approx(0.0, abs=1e-6)

    def test_similar_proteins_empty_sequence(self):
        store = MagicMock()
        service = SimilarityService(feature_store=store)

        with pytest.raises(ValueError, match="must not be empty"):
            service.find_similar_proteins(sequence="", top_k=5)

    def test_similar_proteins_empty_store(self):
        store = MagicMock()
        store.load.return_value = pd.DataFrame()
        service = SimilarityService(feature_store=store)

        with pytest.raises(ValueError, match="No protein embeddings"):
            service.find_similar_proteins(sequence="MKTAYIAK", top_k=5)
