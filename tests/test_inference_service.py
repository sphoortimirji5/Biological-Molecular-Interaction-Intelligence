"""Unit tests for InferenceService — mock model + mock feature store."""
import pytest
import numpy as np
import pandas as pd
from unittest.mock import MagicMock

from src.inference.service import InferenceService
from src.models.scoring import ScoringStrategy


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_protein_df(n: int = 5, emb_dim: int = 4) -> pd.DataFrame:
    """Build a fake protein embeddings DataFrame."""
    rows = []
    for i in range(n):
        row = {
            "external_id": f"P{i:05d}",
            "source": "uniprot",
        }
        emb = np.random.rand(emb_dim)
        for j in range(emb_dim):
            row[f"emb_{j}"] = emb[j]
        rows.append(row)
    return pd.DataFrame(rows)


def _make_service(
    protein_df: pd.DataFrame | None = None,
    score_fn=None,
) -> InferenceService:
    """Factory: InferenceService with mocked feature store + scoring strategy."""
    if protein_df is None:
        protein_df = _make_protein_df()

    mock_store = MagicMock()
    mock_store.load.return_value = protein_df

    mock_scorer = MagicMock(spec=ScoringStrategy)
    if score_fn is None:
        # Default: random scores
        mock_scorer.score.side_effect = lambda q, c: np.random.rand(c.shape[0])
    else:
        mock_scorer.score.side_effect = score_fn

    return InferenceService(
        feature_store=mock_store,
        scoring_strategy=mock_scorer,
        model_version="test_v1",
    )


# Valid aspirin SMILES for all tests
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_rank_returns_top_k():
    """Happy path: returns exactly top_k results."""
    service = _make_service(protein_df=_make_protein_df(n=10))
    results = service.rank(ASPIRIN, top_k=3)

    assert len(results) == 3
    assert all("rank" in r and "external_id" in r and "score" in r for r in results)


def test_rank_invalid_smiles():
    """Invalid SMILES raises ValueError."""
    service = _make_service()
    with pytest.raises(ValueError, match="Invalid SMILES"):
        service.rank("NOT_A_MOLECULE", top_k=5)


def test_rank_respects_top_k():
    """Requesting top_k=2 from 5 candidates returns exactly 2."""
    service = _make_service(protein_df=_make_protein_df(n=5))
    results = service.rank(ASPIRIN, top_k=2)
    assert len(results) == 2


def test_rank_top_k_larger_than_candidates():
    """If top_k > candidates, return all candidates."""
    service = _make_service(protein_df=_make_protein_df(n=3))
    results = service.rank(ASPIRIN, top_k=100)
    assert len(results) == 3


def test_rank_scores_are_descending():
    """Results must be sorted highest score first."""
    # Deterministic scores: [0.1, 0.5, 0.9, 0.3, 0.7]
    fixed_scores = np.array([0.1, 0.5, 0.9, 0.3, 0.7])

    service = _make_service(
        protein_df=_make_protein_df(n=5),
        score_fn=lambda q, c: fixed_scores,
    )

    results = service.rank(ASPIRIN, top_k=5)
    scores = [r["score"] for r in results]

    assert scores == sorted(scores, reverse=True)
    assert scores[0] == pytest.approx(0.9)
    assert scores[-1] == pytest.approx(0.1)


def test_rank_empty_feature_store():
    """Empty protein embeddings raises ValueError."""
    empty_df = pd.DataFrame(columns=["external_id", "source", "emb_0"])
    service = _make_service(protein_df=empty_df)

    with pytest.raises(ValueError, match="No protein embeddings"):
        service.rank(ASPIRIN, top_k=5)


def test_rank_preserves_protein_identity():
    """Returned external_ids must match actual proteins in the feature store."""
    protein_df = _make_protein_df(n=3)
    expected_ids = set(protein_df["external_id"].tolist())

    service = _make_service(protein_df=protein_df)
    results = service.rank(ASPIRIN, top_k=3)
    returned_ids = {r["external_id"] for r in results}

    assert returned_ids == expected_ids


def test_rank_passes_correct_vectors_to_scorer():
    """Scorer receives query fingerprint + protein embedding matrix."""
    captured = {}

    def capture_score(query_vec, candidate_mat):
        captured["query_shape"] = query_vec.shape
        captured["candidate_shape"] = candidate_mat.shape
        return np.ones(candidate_mat.shape[0])

    protein_df = _make_protein_df(n=4, emb_dim=6)
    service = _make_service(protein_df=protein_df, score_fn=capture_score)
    service.rank(ASPIRIN, top_k=2)

    # Query should be 2048-dim Morgan fingerprint
    assert captured["query_shape"] == (2048,)
    # Candidates should be 4 proteins × 6 embedding dims
    assert captured["candidate_shape"] == (4, 6)
