from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional
import numpy as np

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
        # Simple dot product for normalized vectors (Cosine Similarity)
        # Ensure vectors are normalized before calling potential optimization
        # For now, raw dot product
        return np.dot(candidate_vectors, query_vector)

def get_scoring_strategy() -> ScoringStrategy:
    """
    Factory to load the configured scoring strategy.
    """
    from src.config import settings
    
    if settings.scoring_strategy == "similarity":
        return SimilarityScorer()
    
    # Future extension:
    # if settings.scoring_strategy == "learned":
    #     return LearnedScorer(model_path=settings.model_path)
        
    raise ValueError(f"Unknown scoring strategy: {settings.scoring_strategy}")
