"""
TrainingDataBuilder — assembles feature matrices for model training.

Joins compound fingerprints and protein embeddings on the interactions
table to produce a single feature matrix with train/val/test splits.
"""
import structlog
import pandas as pd
import numpy as np
from typing import Dict, Any, Tuple, Optional
from sklearn.model_selection import train_test_split

from sqlalchemy import select, text

from src.db.database import AsyncSessionLocal
from src.db.models import Interaction, Compound, Protein
from src.features.store import FeatureStore

logger = structlog.get_logger(__name__)


class TrainingDataBuilder:
    """Build labeled feature matrices by joining features on interactions."""

    def __init__(
        self,
        feature_store: FeatureStore,
        test_size: float = 0.15,
        val_size: float = 0.15,
        random_state: int = 42,
    ):
        self.feature_store = feature_store
        self.test_size = test_size
        self.val_size = val_size
        self.random_state = random_state

    async def build(
        self,
        compound_feature_version: str = "morgan_v1",
        protein_feature_version: str = "esm2_v1",
    ) -> Dict[str, pd.DataFrame]:
        """
        Build train/val/test DataFrames.

        Returns dict with keys: 'train', 'val', 'test', 'metadata'.
        Each DataFrame has columns: [label, fp_0..fp_N, emb_0..emb_M]
        """
        logger.info(
            "building_training_data",
            compound_version=compound_feature_version,
            protein_version=protein_feature_version,
        )

        # 1. Load pre-computed features from S3
        compound_features = self.feature_store.load("compound", compound_feature_version)
        protein_features = self.feature_store.load("protein", protein_feature_version)

        # 2. Query interactions with compound/protein external_ids
        interactions_df = await self._query_interactions()

        if interactions_df.empty:
            logger.warning("no_interactions_found")
            return {"train": pd.DataFrame(), "val": pd.DataFrame(), "test": pd.DataFrame(), "metadata": {}}

        # 3. Merge: interactions → compound features → protein features
        merged = self._merge_features(interactions_df, compound_features, protein_features)

        if merged.empty:
            logger.warning("no_matched_features", interactions=len(interactions_df))
            return {"train": pd.DataFrame(), "val": pd.DataFrame(), "test": pd.DataFrame(), "metadata": {}}

        # 4. Split into train/val/test
        train_df, val_df, test_df = self._stratified_split(merged)

        metadata = {
            "total_samples": len(merged),
            "train_samples": len(train_df),
            "val_samples": len(val_df),
            "test_samples": len(test_df),
            "compound_feature_version": compound_feature_version,
            "protein_feature_version": protein_feature_version,
            "feature_dim": len([c for c in merged.columns if c.startswith(("fp_", "emb_"))]),
            "label_distribution": merged["label"].value_counts().to_dict(),
        }

        logger.info("training_data_built", **metadata)
        return {"train": train_df, "val": val_df, "test": test_df, "metadata": metadata}

    async def _query_interactions(self) -> pd.DataFrame:
        """Fetch interactions with joined compound/protein external_ids."""
        async with AsyncSessionLocal() as session:
            result = await session.execute(
                text("""
                    SELECT
                        i.label,
                        c.external_id AS compound_external_id,
                        c.source AS compound_source,
                        p.external_id AS protein_external_id,
                        p.source AS protein_source
                    FROM interactions i
                    JOIN compounds c ON i.compound_id = c.id
                    JOIN proteins p ON i.protein_id = p.id
                """)
            )
            rows = result.all()

        if not rows:
            return pd.DataFrame()

        return pd.DataFrame(
            rows,
            columns=["label", "compound_external_id", "compound_source",
                      "protein_external_id", "protein_source"],
        )

    def _merge_features(
        self,
        interactions: pd.DataFrame,
        compound_features: pd.DataFrame,
        protein_features: pd.DataFrame,
    ) -> pd.DataFrame:
        """Join interaction labels with compound + protein feature vectors."""
        # Merge compound features
        merged = interactions.merge(
            compound_features,
            left_on=["compound_external_id", "compound_source"],
            right_on=["external_id", "source"],
            how="inner",
        )

        # Merge protein features
        merged = merged.merge(
            protein_features,
            left_on=["protein_external_id", "protein_source"],
            right_on=["external_id", "source"],
            how="inner",
            suffixes=("_compound", "_protein"),
        )

        # Keep only label + feature columns
        feature_cols = [c for c in merged.columns if c.startswith(("fp_", "emb_"))]
        result = merged[["label"] + feature_cols].copy()

        logger.info(
            "features_merged",
            interactions=len(interactions),
            matched=len(result),
            feature_cols=len(feature_cols),
        )
        return result

    def _stratified_split(
        self, df: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Stratified train/val/test split preserving label distribution."""
        labels = df["label"]

        # First split: train+val vs test
        train_val, test = train_test_split(
            df,
            test_size=self.test_size,
            stratify=labels,
            random_state=self.random_state,
        )

        # Second split: train vs val (from train+val portion)
        val_fraction = self.val_size / (1 - self.test_size)
        train, val = train_test_split(
            train_val,
            test_size=val_fraction,
            stratify=train_val["label"],
            random_state=self.random_state,
        )

        return train.reset_index(drop=True), val.reset_index(drop=True), test.reset_index(drop=True)
