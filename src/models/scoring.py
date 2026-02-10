from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional
from pathlib import Path
import numpy as np
import structlog

logger = structlog.get_logger(__name__)


class ScoringStrategy(ABC):
    """
    Abstract Base Class for interaction scoring.
    Allows swapping between deterministic similarity (Phase 1)
    and learned models (Phase 3) without changing the pipeline.
    """

    @abstractmethod
    def score(self, query_vector: np.ndarray, candidate_vectors: np.ndarray) -> np.ndarray:
        """
        Computes scores for a query against a batch of candidates.

        Args:
            query_vector: Feature vector of the query (e.g. drug embedding)
            candidate_vectors: Matrix of candidate vectors (e.g. protein embeddings)

        Returns:
            Array of scores (higher is better/more similar)
        """
        pass


class SimilarityScorer(ScoringStrategy):
    """
    Default Phase 1 Scorer.
    Uses Cosine Similarity (deterministic).
    """
    def score(self, query_vector: np.ndarray, candidate_vectors: np.ndarray) -> np.ndarray:
        return np.dot(candidate_vectors, query_vector)


class LearnedScorer(ScoringStrategy):
    """
    Phase 3 Scorer — loads a trained XGBoost model for ranking.

    Concatenates query_vector (compound) with each candidate_vector (protein)
    and predicts interaction probability. Higher probability = better match.
    """

    def __init__(self, model_path: str | Path):
        from src.models.xgboost_model import XGBoostWrapper
        self._wrapper = XGBoostWrapper()
        self._wrapper.load(Path(model_path))
        logger.info("learned_scorer_loaded", model_path=str(model_path))

    def score(self, query_vector: np.ndarray, candidate_vectors: np.ndarray) -> np.ndarray:
        """
        Score by concatenating [query | candidate] and predicting probability.

        Args:
            query_vector: 1D compound feature vector (e.g. 2048-dim Morgan fp).
            candidate_vectors: 2D matrix of protein vectors (N × emb_dim).

        Returns:
            Array of interaction probabilities (N,).
        """
        import pandas as pd

        n_candidates = candidate_vectors.shape[0]
        query_repeated = np.tile(query_vector, (n_candidates, 1))
        combined = np.hstack([query_repeated, candidate_vectors])

        # Build DataFrame with feature column names matching training
        fp_cols = [f"fp_{i}" for i in range(query_vector.shape[0])]
        emb_cols = [f"emb_{i}" for i in range(candidate_vectors.shape[1])]
        df = pd.DataFrame(combined, columns=fp_cols + emb_cols)

        result = self._wrapper.predict(df)
        return result["probability"].values


def get_scoring_strategy() -> ScoringStrategy:
    """
    Factory to load the configured scoring strategy.
    """
    from src.config import settings

    if settings.scoring_strategy == "similarity":
        return SimilarityScorer()

    if settings.scoring_strategy == "learned":
        return LearnedScorer(model_path=settings.model_path)

    raise ValueError(f"Unknown scoring strategy: {settings.scoring_strategy}")
