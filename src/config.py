from pydantic_settings import BaseSettings, SettingsConfigDict
from pydantic import computed_field

class Settings(BaseSettings):
    model_config = SettingsConfigDict(env_file=".env", env_ignore_empty=True, extra="ignore")

    # Application
    app_env: str
    app_port: int
    log_level: str

    # Database (Sync URL for Alembic, computed Async for App)
    database_url: str

    # Scoring
    scoring_strategy: str = "similarity"  # "similarity" or "learned"

    # Object Storage (S3-compatible)
    object_store_type: str = "s3"  # Provider key resolved by get_storage_provider()
    object_store_endpoint: str
    object_store_access_key: str
    object_store_secret_key: str
    object_store_bucket_raw: str
    object_store_bucket_processed: str
    object_store_bucket_artifacts: str
    object_store_region: str = "us-east-1"
    object_store_use_ssl: bool = False
    object_store_versioning: bool = False

    # Observability
    otel_exporter: str = "console"
    otel_service_name: str = "bio-interaction-intelligence"
    otel_exporter_otlp_endpoint: str = "http://localhost:4317"

    # Feature Generation
    feature_store_bucket: str = "features"
    esm2_model_name: str = "facebook/esm2_t6_8M_UR50D"
    morgan_radius: int = 2
    morgan_n_bits: int = 2048

    # Default source URLs (override via environment variables)
    chembl_sdf_url: str = ""
    uniprot_fasta_url: str = ""

    @computed_field
    @property
    def async_database_url(self) -> str:
        """
        Converts the standard sync DATABASE_URL (postgresql+psycopg2) 
        to an async driver URL (postgresql+asyncpg) for the app.
        """
        return self.database_url.replace("postgresql+psycopg2://", "postgresql+asyncpg://")

settings = Settings()
