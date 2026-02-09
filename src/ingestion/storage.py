import logging

# Configure logging (or rely on root logger configured in main)
logger = logging.getLogger(__name__)

class MinIOStorage:
    """
    MinIO storage wrapper.
    Enforces layout: raw/<source>/<version>/<sha256>/<file>
    """

    def __init__(self):
        # Assume env is correct. Fail fast if not.
        self.client = Minio(
            settings.object_store_endpoint,
            access_key=settings.object_store_access_key,
            secret_key=settings.object_store_secret_key,
            secure=settings.object_store_use_ssl,
            region=settings.object_store_region
        )
        self._ensure_buckets()

    def _ensure_buckets(self):
        """Ensure all required buckets exist."""
        buckets = [
            settings.object_store_bucket_raw,
            settings.object_store_bucket_processed,
            settings.object_store_bucket_artifacts
        ]
        for bucket in buckets:
            try:
                if not self.client.bucket_exists(bucket):
                    self.client.make_bucket(bucket)
            except Exception as e:
                logger.warning(f"Could not verify bucket {bucket}: {e}")

    def upload_raw(self, source: str, version: str, file_path: Path) -> Dict[str, Any]:
        """
        Uploads a raw file to: raw/<source>/<version>/<sha256>/<filename>
        Returns structured metadata (bucket, key, sha256).
        """
        sha256 = self._calculate_checksum(file_path)
        filename = file_path.name
        object_name = f"{source}/{version}/{sha256}/{filename}"
        
        metadata = {
            "source": source, 
            "version": version, 
            "sha256": sha256
        }

        self.client.fput_object(
            settings.object_store_bucket_raw,
            object_name,
            str(file_path),
            metadata=metadata
        )
        
        return {
            "bucket": settings.object_store_bucket_raw,
            "key": object_name,
            "sha256": sha256
        }

    def _calculate_checksum(self, file_path: Path) -> str:
        """Calculates SHA256 checksum of a file."""
        sha256_hash = hashlib.sha256()
        with open(file_path, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
        return sha256_hash.hexdigest()
