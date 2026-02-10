"""Unit tests for scoring strategies — mock model, verify scoring output."""
import pytest
import numpy as np
import pandas as pd
from unittest.mock import MagicMock, patch
from pathlib import Path

from src.models.scoring import SimilarityScorer, LearnedScorer, ScoringStrategy


def test_similarity_scorer():
    """SimilarityScorer produces dot-product scores."""
    scorer = SimilarityScorer()
    query = np.array([1.0, 0.0, 0.0])
    candidates = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.5, 0.5, 0.0],
    ])

    scores = scorer.score(query, candidates)

    assert scores[0] == pytest.approx(1.0)
    assert scores[1] == pytest.approx(0.0)
    assert scores[2] == pytest.approx(0.5)


def _make_mock_scorer(predict_fn=None):
    """Factory: builds a LearnedScorer with a mocked internal wrapper."""
    mock_wrapper = MagicMock()

    if predict_fn is None:
        def predict_fn(df):
            result = df.copy()
            result["probability"] = np.random.rand(len(df))
            result["prediction"] = (result["probability"] >= 0.5).astype(int)
            return result

    mock_wrapper.predict.side_effect = predict_fn
    mock_wrapper.load.return_value = mock_wrapper

    scorer = LearnedScorer.__new__(LearnedScorer)
    scorer._wrapper = mock_wrapper
    return scorer, mock_wrapper


def test_learned_scorer_produces_probabilities():
    """LearnedScorer must return probability scores between 0 and 1."""
    scorer, _ = _make_mock_scorer()

    query = np.random.rand(4)
    candidates = np.random.rand(5, 3)

    scores = scorer.score(query, candidates)

    assert len(scores) == 5
    assert all(0 <= s <= 1 for s in scores)


def test_learned_scorer_concatenates_features():
    """LearnedScorer must pass [query | candidate] to the model."""
    received_dfs = []

    def capture_predict(df):
        received_dfs.append(df.copy())
        result = df.copy()
        result["probability"] = np.full(len(df), 0.5)
        result["prediction"] = np.zeros(len(df), dtype=int)
        return result

    scorer, _ = _make_mock_scorer(predict_fn=capture_predict)

    query = np.array([1.0, 2.0])
    candidates = np.array([[3.0, 4.0, 5.0], [6.0, 7.0, 8.0]])

    scorer.score(query, candidates)

    df = received_dfs[0]
    assert df.shape == (2, 5)  # 2 candidates × (2 fp + 3 emb)
    assert list(df.columns) == ["fp_0", "fp_1", "emb_0", "emb_1", "emb_2"]
    np.testing.assert_array_equal(df.iloc[0].values, [1.0, 2.0, 3.0, 4.0, 5.0])


def test_scoring_strategy_is_abstract():
    """Can't instantiate ScoringStrategy directly."""
    with pytest.raises(TypeError):
        ScoringStrategy()
