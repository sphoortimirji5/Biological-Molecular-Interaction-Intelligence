"""
CLI entrypoint for model training.

Usage:
    docker-compose run --rm app python -m src.models.trainer
"""
import asyncio
import json
import structlog

from src.config import settings
from src.ingestion.storage import create_storage_provider
from src.features.store import FeatureStore
from src.models.trainer import TrainingManager

logger = structlog.get_logger(__name__)


async def main():
    storage = create_storage_provider()
    feature_store = FeatureStore(storage=storage, bucket=settings.s3_bucket)

    manager = TrainingManager(
        feature_store=feature_store,
        model_version="v1",
        compound_feature_version="morgan_v1",
        protein_feature_version="esm2_v1",
    )

    result = await manager.run()
    print(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    asyncio.run(main())
