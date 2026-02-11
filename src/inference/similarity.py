"""
SimilarityService — find similar proteins or compounds in the feature store.

Uses cosine similarity over ESM-2 protein embeddings and Morgan
compound fingerprints.
"""
import numpy as np
import pandas as pd
import structlog
from typing import List, Dict, Any

from rdkit import Chem

from src.features.store import FeatureStore
from src.features.morgan import MorganFingerprintGenerator


logger = structlog.get_logger(__name__)


def _cosine_similarity(query: np.ndarray, candidates: np.ndarray) -> np.ndarray:
    """Cosine similarity between a single query vector and many candidates."""
    query_norm = np.linalg.norm(query)
    if query_norm == 0:
        return np.zeros(len(candidates))

    candidate_norms = np.linalg.norm(candidates, axis=1)
    candidate_norms = np.where(candidate_norms == 0, 1.0, candidate_norms)

    return (candidates @ query) / (candidate_norms * query_norm)





class SimilarityService:
    """Find similar entities using embedding/fingerprint similarity."""

    def __init__(
        self,
        feature_store: FeatureStore,
        protein_feature_version: str = "esm2_v1",
        compound_feature_version: str = "morgan_v1",
    ):
        self.feature_store = feature_store
        self.protein_feature_version = protein_feature_version
        self.compound_feature_version = compound_feature_version
        self._fingerprint_gen = MorganFingerprintGenerator()

    def find_similar_proteins(
        self, sequence: str, top_k: int = 10
    ) -> List[Dict[str, Any]]:
        """
        Find proteins most similar to a query sequence.

        Args:
            sequence: Amino acid sequence of the query protein.
            top_k: Number of results to return.

        Returns:
            Ranked list of similar proteins with cosine similarity scores.

        Raises:
            ValueError: If sequence is empty or no proteins are in the store.
        """
        if not sequence or not sequence.strip():
            raise ValueError("Protein sequence must not be empty")

        # Load stored protein embeddings
        protein_df = self.feature_store.load(
            entity_type="protein",
            version=self.protein_feature_version,
        )

        if protein_df.empty:
            raise ValueError("No protein embeddings available in the feature store.")

        emb_cols = [c for c in protein_df.columns if c.startswith("emb_")]
        candidate_vectors = protein_df[emb_cols].values

        # Featurize query protein via ESM-2
        from src.features.protein_embeddings import ProteinEmbeddingGenerator

        embedder = ProteinEmbeddingGenerator()
        query_df = pd.DataFrame([{
            "external_id": "query",
            "source": "similarity_search",
            "sequence": sequence,
        }])
        query_emb_df = embedder.transform(query_df)

        if query_emb_df.empty:
            raise ValueError("Failed to generate embedding for query sequence.")

        query_emb_cols = [c for c in query_emb_df.columns if c.startswith("emb_")]
        query_vector = query_emb_df[query_emb_cols].values[0]

        # Cosine similarity
        scores = _cosine_similarity(query_vector, candidate_vectors)
        return self._build_results(protein_df, scores, top_k)

    def find_similar_compounds(
        self, smiles: str, top_k: int = 10
    ) -> List[Dict[str, Any]]:
        """
        Find compounds most similar to a query SMILES.

        Args:
            smiles: SMILES string of the query compound.
            top_k: Number of results to return.

        Returns:
            Ranked list of similar compounds with cosine similarity scores.

        Raises:
            ValueError: If SMILES is invalid or no compounds are in the store.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")

        # Load stored compound fingerprints
        compound_df = self.feature_store.load(
            entity_type="compound",
            version=self.compound_feature_version,
        )

        if compound_df.empty:
            raise ValueError("No compound fingerprints available in the feature store.")

        fp_cols = [c for c in compound_df.columns if c.startswith("fp_")]
        candidate_vectors = compound_df[fp_cols].values

        # Featurize query compound
        query_df = pd.DataFrame([{
            "external_id": "query",
            "source": "similarity_search",
            "smiles": smiles,
        }])
        query_fp_df = self._fingerprint_gen.transform(query_df)

        if query_fp_df.empty:
            raise ValueError(f"Failed to featurize SMILES: {smiles}")

        query_fp_cols = [c for c in query_fp_df.columns if c.startswith("fp_")]
        query_vector = query_fp_df[query_fp_cols].values[0]

        # Cosine similarity
        scores = _cosine_similarity(query_vector, candidate_vectors)
        return self._build_results(compound_df, scores, top_k)

    @staticmethod
    def _build_results(
        df: pd.DataFrame, scores: np.ndarray, top_k: int
    ) -> List[Dict[str, Any]]:
        """Sort by score descending, return top-K."""
        top_indices = np.argsort(scores)[::-1][:top_k]

        results = []
        for rank_pos, idx in enumerate(top_indices, start=1):
            results.append({
                "rank": rank_pos,
                "external_id": str(df.iloc[idx]["external_id"]),
                "source": str(df.iloc[idx]["source"]),
                "score": float(scores[idx]),
            })

        logger.info(
            "similarity_search_completed",
            candidates=len(df),
            top_k=top_k,
            top_score=results[0]["score"] if results else None,
        )

        return results
