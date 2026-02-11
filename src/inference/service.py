"""
InferenceService — stateless orchestrator for drug-target ranking.

Accepts a SMILES string, featurizes it on the fly, scores against
all pre-computed protein embeddings, and returns the top-K ranked targets.
"""
import structlog
import numpy as np
import pandas as pd
from typing import List, Dict, Any

from rdkit import Chem

from src.features.morgan import MorganFingerprintGenerator
from src.features.store import FeatureStore
from src.models.scoring import ScoringStrategy

logger = structlog.get_logger(__name__)


class InferenceService:
    """Orchestrate featurization → scoring → ranking for a query compound."""

    def __init__(
        self,
        feature_store: FeatureStore,
        scoring_strategy: ScoringStrategy,
        protein_feature_version: str = "esm2_v1",
        model_version: str = "v1",
    ):
        self.feature_store = feature_store
        self.scoring_strategy = scoring_strategy
        self.protein_feature_version = protein_feature_version
        self.model_version = model_version
        self._fingerprint_gen = MorganFingerprintGenerator()

    def rank(self, smiles: str, top_k: int = 10) -> List[Dict[str, Any]]:
        """
        Rank protein targets for a query compound.

        Args:
            smiles: SMILES string of the query compound.
            top_k: Number of top-ranked targets to return.

        Returns:
            List of dicts with rank, external_id, source, and score.

        Raises:
            ValueError: If the SMILES string is invalid or no proteins are available.
        """
        # 1. Validate SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")

        # 2. Featurize query compound (Morgan fingerprint)
        query_df = pd.DataFrame([{
            "external_id": "query",
            "source": "inference",
            "smiles": smiles,
        }])
        fp_df = self._fingerprint_gen.transform(query_df)

        if fp_df.empty:
            raise ValueError(f"Failed to featurize SMILES: {smiles}")

        fp_cols = [c for c in fp_df.columns if c.startswith("fp_")]
        query_vector = fp_df[fp_cols].values[0]

        # 3. Load pre-computed protein embeddings
        protein_df = self.feature_store.load(
            entity_type="protein",
            version=self.protein_feature_version,
        )

        if protein_df.empty:
            raise ValueError("No protein embeddings available in the feature store.")

        emb_cols = [c for c in protein_df.columns if c.startswith("emb_")]
        candidate_vectors = protein_df[emb_cols].values

        # 4. Score all candidates
        scores = self.scoring_strategy.score(query_vector, candidate_vectors)

        # 5. Rank and return top-K
        top_indices = np.argsort(scores)[::-1][:top_k]

        results = []
        for rank_pos, idx in enumerate(top_indices, start=1):
            results.append({
                "rank": rank_pos,
                "external_id": str(protein_df.iloc[idx]["external_id"]),
                "source": str(protein_df.iloc[idx]["source"]),
                "score": float(scores[idx]),
            })

        logger.info(
            "inference_completed",
            smiles=smiles[:50],
            candidates=len(protein_df),
            top_k=top_k,
            top_score=results[0]["score"] if results else None,
        )

        return results
