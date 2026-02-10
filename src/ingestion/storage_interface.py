"""
ObjectStorageProvider — vendor-neutral interface for object storage.

Implementations must handle:
    - Content-addressed uploads (raw/<source>/<version>/<sha256>/<file>)
    - Downloads to local ephemeral paths
    - Existence checks for idempotent re-ingestion

The active provider is resolved at runtime from ``settings.object_store_type``.
"""
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Any, IO
import hashlib


class ObjectStorageProvider(ABC):
    """Abstract base for S3-compatible object storage backends."""

    @abstractmethod
    def upload_raw(self, source: str, version: str, file_path: Path) -> Dict[str, Any]:
        """
        Upload a raw artifact to the content-addressed layout.

        Expected key format: ``<source>/<version>/<sha256>/<filename>``

        Returns:
            dict with ``bucket``, ``key``, and ``sha256`` fields.
        """

    @abstractmethod
    def download(self, bucket: str, key: str, dest: Path) -> Path:
        """
        Download an object to *dest*.

        Returns:
            The local path of the downloaded file.
        """

    @abstractmethod
    def exists(self, bucket: str, key: str) -> bool:
        """Check whether an object exists in the given bucket."""

    @abstractmethod
    def put_object(self, bucket: str, key: str, body: IO, content_type: str = "application/octet-stream") -> None:
        """Upload a file-like object."""

    @abstractmethod
    def get_object(self, bucket: str, key: str) -> IO:
        """Download an object as a file-like stream."""

    @staticmethod
    def calculate_checksum(file_path: Path) -> str:
        """Compute SHA-256 digest for integrity verification."""
        sha256_hash = hashlib.sha256()
        with open(file_path, "rb") as f:
            for block in iter(lambda: f.read(8192), b""):
                sha256_hash.update(block)
        return sha256_hash.hexdigest()
