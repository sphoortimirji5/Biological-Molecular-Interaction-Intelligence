"""
S3StorageProvider — boto3-backed implementation of ObjectStorageProvider.

Works with any S3-compatible backend (AWS S3, GCS, MinIO) by pointing
``OBJECT_STORE_ENDPOINT`` to the appropriate service URL.

Also exposes ``get_storage_provider()`` factory, which resolves the active
provider from ``settings.object_store_type``.
"""
import logging
from src.logging_config import get_logger
from pathlib import Path
from typing import Dict, Any

import boto3
from botocore.exceptions import ClientError

from src.config import settings
from src.ingestion.storage_interface import ObjectStorageProvider

logger = get_logger(__name__)


class S3StorageProvider(ObjectStorageProvider):
    """
    S3-compatible object storage provider.

    Enforces a content-addressed layout for reproducibility::

        raw/<source>/<version>/<sha256>/<filename>

    Backed by any S3-compatible service — set ``OBJECT_STORE_ENDPOINT``
    to ``http://minio:9000`` locally or omit for AWS S3.
    """

    def __init__(self) -> None:
        """Initialize boto3 client from environment configuration."""
        endpoint_url = self._resolve_endpoint()
        self.client = boto3.client(
            "s3",
            endpoint_url=endpoint_url,
            aws_access_key_id=settings.object_store_access_key,
            aws_secret_access_key=settings.object_store_secret_key,
            region_name=settings.object_store_region,
        )
        self._ensure_buckets()

    @staticmethod
    def _resolve_endpoint() -> str:
        """Build the full endpoint URL from configuration."""
        endpoint = settings.object_store_endpoint
        if not endpoint.startswith(("http://", "https://")):
            scheme = "https" if settings.object_store_use_ssl else "http"
            endpoint = f"{scheme}://{endpoint}"
        return endpoint

    def _ensure_buckets(self) -> None:
        """Create required buckets if they do not already exist."""
        buckets = [
            settings.object_store_bucket_raw,
            settings.object_store_bucket_processed,
            settings.object_store_bucket_artifacts,
        ]
        for bucket in buckets:
            try:
                self.client.head_bucket(Bucket=bucket)
            except ClientError:
                try:
                    self.client.create_bucket(Bucket=bucket)
                except ClientError as e:
                    logger.warning("Could not create bucket %s: %s", bucket, e)

    def upload_raw(self, source: str, version: str, file_path: Path) -> Dict[str, Any]:
        """
        Upload a raw artifact to the content-addressed layout.

        Returns:
            dict with ``bucket``, ``key``, and ``sha256`` fields.
        """
        sha256 = self.calculate_checksum(file_path)
        filename = file_path.name
        object_key = f"{source}/{version}/{sha256}/{filename}"

        self.client.upload_file(
            Filename=str(file_path),
            Bucket=settings.object_store_bucket_raw,
            Key=object_key,
            ExtraArgs={
                "Metadata": {
                    "source": source,
                    "version": version,
                    "sha256": sha256,
                }
            },
        )

        return {
            "bucket": settings.object_store_bucket_raw,
            "key": object_key,
            "sha256": sha256,
        }

    def download(self, bucket: str, key: str, dest: Path) -> Path:
        """Download an object to *dest*."""
        dest.parent.mkdir(parents=True, exist_ok=True)
        self.client.download_file(Bucket=bucket, Key=key, Filename=str(dest))
        return dest

    def exists(self, bucket: str, key: str) -> bool:
        """Check whether an object exists in the given bucket."""
        try:
            self.client.head_object(Bucket=bucket, Key=key)
            return True
        except ClientError:
            return False


def get_storage_provider() -> ObjectStorageProvider:
    """
    Factory: resolve the active storage provider from configuration.

    Currently supports ``s3`` (default). Extend with additional
    ``elif`` branches for GCS, Azure Blob, or local filesystem providers.
    """
    provider_type = settings.object_store_type.lower()
    if provider_type == "s3":
        return S3StorageProvider()
    raise ValueError(
        f"Unknown OBJECT_STORE_TYPE: '{provider_type}'. Supported: s3"
    )
