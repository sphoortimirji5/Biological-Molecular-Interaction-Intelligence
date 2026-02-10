"""
FeatureStore — Parquet-based feature persistence on S3.

Saves and loads feature DataFrames with versioned manifests
for reproducibility tracking.
"""
import io
import json
import structlog
import pandas as pd
from typing import Dict, Any

from src.ingestion.storage_interface import ObjectStorageProvider

logger = structlog.get_logger(__name__)


class FeatureStore:
    """Persist and load feature sets as Parquet files in S3."""

    def __init__(self, storage: ObjectStorageProvider, bucket: str):
        self.storage = storage
        self.bucket = bucket

    def save(
        self,
        entity_type: str,
        version: str,
        df: pd.DataFrame,
        manifest: Dict[str, Any],
    ) -> str:
        """Save features + manifest. Returns the S3 URI of the Parquet file."""
        prefix = f"{entity_type}/{version}"

        # 1. Save Parquet
        parquet_key = f"{prefix}/features.parquet"
        buffer = io.BytesIO()
        df.to_parquet(buffer, index=False)
        buffer.seek(0)
        self.storage.put_object(
            bucket=self.bucket,
            key=parquet_key,
            body=buffer,
            content_type="application/vnd.apache.parquet",
        )

        # 2. Save Manifest
        manifest_key = f"{prefix}/manifest.json"
        manifest_json = json.dumps(manifest, indent=2)
        self.storage.put_object(
            bucket=self.bucket,
            key=manifest_key,
            body=io.BytesIO(manifest_json.encode("utf-8")),
            content_type="application/json",
        )

        logger.info(
            "features_saved",
            entity_type=entity_type,
            version=version,
            rows=len(df),
            path=f"s3://{self.bucket}/{parquet_key}",
        )
        return f"s3://{self.bucket}/{parquet_key}"

    def load(self, entity_type: str, version: str) -> pd.DataFrame:
        """Load features from S3 as a DataFrame."""
        key = f"{entity_type}/{version}/features.parquet"
        try:
            body = self.storage.get_object(bucket=self.bucket, key=key)
            return pd.read_parquet(body)
        except Exception as e:
            logger.error(
                "feature_load_failed",
                entity_type=entity_type,
                version=version,
                error=str(e),
            )
            raise

    def load_manifest(self, entity_type: str, version: str) -> Dict[str, Any]:
        """Load the manifest JSON for a feature version."""
        key = f"{entity_type}/{version}/manifest.json"
        body = self.storage.get_object(bucket=self.bucket, key=key)
        return json.loads(body.read().decode("utf-8"))
