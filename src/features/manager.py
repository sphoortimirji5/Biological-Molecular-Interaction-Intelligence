"""
FeaturizationManager — batch orchestrator for molecular featurization.

Queries compounds/proteins from PostgreSQL, routes to the correct
FeatureGenerator, and persists results via FeatureStore.

Usage:
    python -m src.features.manager
"""
import asyncio
import structlog
import pandas as pd
from datetime import datetime, timezone
from typing import Optional

from sqlalchemy import select, func

from src.db.database import AsyncSessionLocal
from src.db.models import Compound, Protein
from src.features.interfaces import FeatureGenerator
from src.features.morgan import MorganFingerprintGenerator
from src.features.protein_embeddings import ProteinEmbeddingGenerator
from src.features.store import FeatureStore
from src.ingestion.storage import get_storage_provider
from src.config import settings

logger = structlog.get_logger(__name__)


class FeaturizationManager:
    """Orchestrate batch featurization of compounds and proteins."""

    def __init__(
        self,
        feature_store: Optional[FeatureStore] = None,
        compound_generator: Optional[FeatureGenerator] = None,
        protein_generator: Optional[FeatureGenerator] = None,
        batch_size: int = 1000,
    ):
        self.feature_store = feature_store or FeatureStore(
            storage=get_storage_provider(),
            bucket=settings.feature_store_bucket,
        )
        self.compound_generator = compound_generator or MorganFingerprintGenerator()
        self.protein_generator = protein_generator or ProteinEmbeddingGenerator()
        self.batch_size = batch_size

    async def run(self, version: str = "v1") -> dict:
        """
        Run full featurization pipeline.

        Returns summary stats dict.
        """
        logger.info("featurization_started", version=version)
        start = datetime.now(timezone.utc)

        compound_stats = await self.featurize_compounds(version)
        protein_stats = await self.featurize_proteins(version)

        elapsed = (datetime.now(timezone.utc) - start).total_seconds()
        summary = {
            "version": version,
            "compounds": compound_stats,
            "proteins": protein_stats,
            "elapsed_seconds": round(elapsed, 2),
        }

        logger.info("featurization_completed", **summary)
        return summary

    async def featurize_compounds(self, version: str) -> dict:
        """Query all compounds and generate Morgan fingerprints."""
        logger.info("featurizing_compounds", version=version)

        df = await self._query_compounds()
        if df.empty:
            logger.warning("no_compounds_found")
            return {"count": 0, "uri": None}

        result_df = self.compound_generator.transform(df)
        manifest = self.compound_generator.feature_manifest()
        manifest["row_count"] = len(result_df)

        uri = self.feature_store.save(
            entity_type="compound",
            version=f"morgan_{version}",
            df=result_df,
            manifest=manifest,
        )

        logger.info(
            "compounds_featurized",
            total=len(df),
            valid=len(result_df),
            uri=uri,
        )
        return {"count": len(result_df), "uri": uri}

    async def featurize_proteins(self, version: str) -> dict:
        """Query all proteins and generate ESM-2 embeddings."""
        logger.info("featurizing_proteins", version=version)

        df = await self._query_proteins()
        if df.empty:
            logger.warning("no_proteins_found")
            return {"count": 0, "uri": None}

        result_df = self.protein_generator.transform(df)
        manifest = self.protein_generator.feature_manifest()
        manifest["row_count"] = len(result_df)

        uri = self.feature_store.save(
            entity_type="protein",
            version=f"esm2_{version}",
            df=result_df,
            manifest=manifest,
        )

        logger.info(
            "proteins_featurized",
            total=len(df),
            valid=len(result_df),
            uri=uri,
        )
        return {"count": len(result_df), "uri": uri}

    async def _query_compounds(self) -> pd.DataFrame:
        """Fetch all compounds from Postgres as a DataFrame."""
        async with AsyncSessionLocal() as session:
            result = await session.execute(
                select(
                    Compound.external_id,
                    Compound.source,
                    Compound.smiles,
                )
            )
            rows = result.all()

        if not rows:
            return pd.DataFrame(columns=["external_id", "source", "smiles"])

        return pd.DataFrame(rows, columns=["external_id", "source", "smiles"])

    async def _query_proteins(self) -> pd.DataFrame:
        """Fetch all proteins from Postgres as a DataFrame."""
        async with AsyncSessionLocal() as session:
            result = await session.execute(
                select(
                    Protein.external_id,
                    Protein.source,
                    Protein.sequence,
                )
            )
            rows = result.all()

        if not rows:
            return pd.DataFrame(columns=["external_id", "source", "sequence"])

        return pd.DataFrame(rows, columns=["external_id", "source", "sequence"])


async def main():
    """CLI entrypoint for batch featurization."""
    manager = FeaturizationManager()
    summary = await manager.run()
    print(f"Featurization complete: {summary}")


if __name__ == "__main__":
    asyncio.run(main())
