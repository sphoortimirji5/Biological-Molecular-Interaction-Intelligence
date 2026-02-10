"""
Integration test: Seed DB → Featurize → Verify Parquet in S3.

Requires running services (Postgres, MinIO) and migrations applied.
Run with:
    docker-compose run --rm app poetry run alembic upgrade head
    docker-compose run --rm app poetry run pytest tests/test_featurization_integration.py -v -m integration
"""
import pytest
import pandas as pd
from sqlalchemy import text

from src.db.database import AsyncSessionLocal
from src.features.manager import FeaturizationManager
from src.features.morgan import MorganFingerprintGenerator
from src.features.store import FeatureStore
from src.ingestion.storage import get_storage_provider
from src.config import settings


# Use real Morgan (fast, CPU-only). Skip ESM-2 in integration tests
# to avoid model download in CI. Protein embeddings are covered by unit tests.
SEED_COMPOUNDS = [
    ("ASPIRIN", "integration_test", "CC(=O)Oc1ccccc1C(=O)O"),
    ("CAFFEINE", "integration_test", "Cn1c(=O)c2c(ncn2C)n(C)c1=O"),
    ("ETHANOL", "integration_test", "CCO"),
]


@pytest.fixture
def storage():
    return get_storage_provider()


@pytest.fixture
def feature_store(storage):
    return FeatureStore(storage=storage, bucket=settings.feature_store_bucket)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_compound_featurization_end_to_end(feature_store, storage, db_cleanup):
    """Seed compounds → featurize → verify Parquet exists in S3."""

    # 1. Seed compounds into Postgres
    async with AsyncSessionLocal() as session:
        for ext_id, source, smiles in SEED_COMPOUNDS:
            await session.execute(
                text(
                    "INSERT INTO compounds (source, external_id, smiles) "
                    "VALUES (:source, :ext_id, :smiles) "
                    "ON CONFLICT (source, external_id) DO NOTHING"
                ),
                {"source": source, "ext_id": ext_id, "smiles": smiles},
            )
        await session.commit()

    # 2. Verify seed
    async with AsyncSessionLocal() as session:
        result = await session.execute(
            text("SELECT COUNT(*) FROM compounds WHERE source = 'integration_test'")
        )
        count = result.scalar()
        assert count == 3, f"Expected 3 seeded compounds, got {count}"

    # 3. Run featurization (compounds only, skip proteins)
    manager = FeaturizationManager(
        feature_store=feature_store,
        compound_generator=MorganFingerprintGenerator(),
        protein_generator=None,  # Skip protein embedding
    )
    result = await manager.featurize_compounds("integration_test")

    # 4. Verify results
    assert result["count"] == 3
    assert result["uri"] is not None
    assert "integration_test" in result["uri"]

    # 5. Verify Parquet is loadable from S3
    loaded_df = feature_store.load("compound", "morgan_integration_test")
    assert len(loaded_df) == 3
    assert "fp_0" in loaded_df.columns
    assert "external_id" in loaded_df.columns

    # 6. Verify manifest
    manifest = feature_store.load_manifest("compound", "morgan_integration_test")
    assert manifest["type"] == "morgan"
    assert manifest["radius"] == 2
    assert manifest["n_bits"] == 2048
    assert manifest["row_count"] == 3
