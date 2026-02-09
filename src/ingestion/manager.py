"""
IngestionManager — orchestrates the data ingestion pipeline.

Stages:
    1. Discover — resolve source metadata (version, format, entity type)
    2. Fetch    — download raw data to ephemeral storage
    3. Store    — upload to S3-compatible object store (content-addressed)
    4. Parse    — stream records through format-driven parser
    5. Upsert   — idempotent write to PostgreSQL via ON CONFLICT

Routing is fully data-driven: the ``format`` field selects the parser,
the ``entity_type`` field selects the upsert target. No source-specific
branches exist in this module.

Registering a proprietary source
---------------------------------

.. code-block:: python

    from src.ingestion.sources import HTTPSource
    from src.ingestion.manager import register_source, IngestionManager
    import os

    # 1. Register — URL from environment, format drives the pipeline.
    register_source("pharma_x", HTTPSource(
        source="pharma_x",
        url=os.environ["PHARMA_X_SDF_URL"],
        fmt="sdf",               # drives parser selection
        entity_type="compound",  # drives upsert target
        version="2024Q1",
    ))

    # 2. Run — identical pipeline, zero code changes.
    manager = IngestionManager()
    await manager.run_pipeline("pharma_x")
"""
import time
import tempfile
from pathlib import Path
from typing import Dict, Any, List
from datetime import datetime

from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.dialects.postgresql import insert

from src.db.database import AsyncSessionLocal
from src.db.models import IngestionRun, Protein, Compound
from src.ingestion.interfaces import DataSource
from src.ingestion.storage_interface import ObjectStorageProvider
from src.ingestion.storage import get_storage_provider
from src.ingestion.parsers import get_parser
from src.logging_config import get_logger
from src.metrics import (
    RUNS_TOTAL,
    RECORDS_UPSERTED,
    STAGE_DURATION,
    BATCH_SIZE,
)
from src.tracing import get_tracer

logger = get_logger(__name__)
tracer = get_tracer(__name__)

# ---------------------------------------------------------------------------
# Source Registry — maps source labels to DataSource instances.
# Register built-in and proprietary sources here.
# ---------------------------------------------------------------------------
_SOURCE_REGISTRY: Dict[str, DataSource] = {}


def register_source(name: str, source: DataSource) -> None:
    """Register a DataSource under *name*. Overwrites existing entries."""
    _SOURCE_REGISTRY[name] = source


def get_registered_source(name: str) -> DataSource:
    """Retrieve a registered source by name. Raises ValueError if unknown."""
    source = _SOURCE_REGISTRY.get(name)
    if source is None:
        raise ValueError(
            f"Unknown source: '{name}'. Registered: {list(_SOURCE_REGISTRY.keys())}"
        )
    return source


# ---------------------------------------------------------------------------
# Register built-in public sources (imported lazily to avoid circular deps)
# ---------------------------------------------------------------------------
def _register_defaults() -> None:
    """Register ChEMBL and UniProt with default configuration."""
    from src.ingestion.sources import chembl_source, uniprot_source

    if "chembl" not in _SOURCE_REGISTRY:
        register_source("chembl", chembl_source())
    if "uniprot" not in _SOURCE_REGISTRY:
        register_source("uniprot", uniprot_source())


_register_defaults()


class IngestionManager:
    """
    Orchestrates the data ingestion pipeline.

    Routing is driven by ``DataSource.discover()`` metadata:
    - ``format`` → selects parser (sdf, fasta, csv)
    - ``entity_type`` → selects upsert target (compound, protein)
    """

    def __init__(self, storage: ObjectStorageProvider | None = None) -> None:
        self.storage = storage or get_storage_provider()

    async def run_pipeline(self, source_name: str) -> None:
        """
        Execute the full ingestion pipeline for a registered data source.

        Creates an ``ingestion_runs`` audit record, then progresses through
        Discover → Fetch → Store → Parse → Upsert. On failure the run is
        marked FAILED with error details before re-raising.
        """
        with tracer.start_as_current_span("ingestion.run_pipeline") as root_span:
            root_span.set_attribute("source", source_name)
            pipeline_start = time.monotonic()

            async with AsyncSessionLocal() as session:
                run = IngestionRun(
                    status="STARTED",
                    source=source_name,
                    started_at=datetime.utcnow(),
                )
                session.add(run)
                await session.commit()
                await session.refresh(run)

                logger.info(
                    "pipeline_started",
                    source=source_name,
                    run_id=str(run.run_id),
                )

                try:
                    # -- Discover --
                    with tracer.start_as_current_span("ingestion.discover") as span:
                        t0 = time.monotonic()
                        source = get_registered_source(source_name)
                        discovery = source.discover()
                        span.set_attribute("format", discovery["format"])
                        span.set_attribute("entity_type", discovery["entity_type"])
                        STAGE_DURATION.labels(source=source_name, stage="discover").observe(
                            time.monotonic() - t0
                        )

                    run.stats = discovery.get("stats")
                    session.add(run)
                    await session.commit()

                    fmt = discovery["format"]
                    entity_type = discovery["entity_type"]

                    logger.info(
                        "discovery_complete",
                        source=source_name,
                        format=fmt,
                        entity_type=entity_type,
                        version=discovery.get("version"),
                    )

                    with tempfile.TemporaryDirectory() as temp_dir:
                        temp_path = Path(temp_dir) / f"{source_name}_raw.gz"

                        # -- Fetch --
                        with tracer.start_as_current_span("ingestion.fetch") as span:
                            t0 = time.monotonic()
                            downloaded_path = source.fetch(temp_path)
                            fetch_duration = time.monotonic() - t0
                            file_size = downloaded_path.stat().st_size
                            span.set_attribute("file_size_bytes", file_size)
                            STAGE_DURATION.labels(source=source_name, stage="fetch").observe(
                                fetch_duration
                            )

                        logger.info(
                            "fetch_complete",
                            source=source_name,
                            bytes=file_size,
                            duration_s=round(fetch_duration, 2),
                        )

                        # -- Store --
                        with tracer.start_as_current_span("ingestion.store") as span:
                            t0 = time.monotonic()
                            meta = self.storage.upload_raw(
                                source_name,
                                discovery.get("version", "unknown"),
                                downloaded_path,
                            )
                            STAGE_DURATION.labels(source=source_name, stage="store").observe(
                                time.monotonic() - t0
                            )
                            span.set_attribute("sha256", meta["sha256"])

                        run.checksums = {"sha256": meta["sha256"]}
                        session.add(run)
                        await session.commit()

                        logger.info(
                            "store_complete",
                            source=source_name,
                            sha256=meta["sha256"],
                        )

                        # -- Parse + Upsert --
                        with tracer.start_as_current_span("ingestion.parse_and_upsert") as span:
                            t0 = time.monotonic()
                            parser = get_parser(fmt)
                            total_records = await self._process(
                                session, parser, downloaded_path, source_name, entity_type
                            )
                            STAGE_DURATION.labels(
                                source=source_name, stage="parse_and_upsert"
                            ).observe(time.monotonic() - t0)
                            span.set_attribute("total_records", total_records)

                        logger.info(
                            "parse_upsert_complete",
                            source=source_name,
                            records=total_records,
                            entity_type=entity_type,
                        )

                    run.status = "COMPLETED"
                    run.ended_at = datetime.utcnow()
                    session.add(run)
                    await session.commit()

                    total_duration = time.monotonic() - pipeline_start
                    RUNS_TOTAL.labels(source=source_name, status="COMPLETED").inc()

                    logger.info(
                        "pipeline_completed",
                        source=source_name,
                        run_id=str(run.run_id),
                        duration_s=round(total_duration, 2),
                        records=total_records,
                    )

                except Exception as e:
                    RUNS_TOTAL.labels(source=source_name, status="FAILED").inc()
                    root_span.set_attribute("error", True)
                    root_span.set_attribute("error.message", str(e))

                    logger.error(
                        "pipeline_failed",
                        source=source_name,
                        run_id=str(run.run_id),
                        error=str(e),
                        exc_info=True,
                    )

                    run.status = "FAILED"
                    run.ended_at = datetime.utcnow()
                    run.stats = {"error": str(e)}
                    session.add(run)
                    await session.commit()
                    raise

    # ------------------------------------------------------------------
    # Unified processing — driven by entity_type, not source name
    # ------------------------------------------------------------------

    async def _process(
        self,
        session: AsyncSession,
        parser,
        file_path: Path,
        source_name: str,
        entity_type: str,
        batch_size: int = 1000,
    ) -> int:
        """Parse file and batch-upsert to the table matching *entity_type*.

        Returns the total number of records processed.
        """
        batch: List[Dict[str, Any]] = []
        total = 0

        for record in parser.parse(file_path):
            row = self._to_row(record, source_name, entity_type)
            batch.append(row)
            if len(batch) >= batch_size:
                await self._upsert(session, entity_type, batch, source_name)
                total += len(batch)
                batch = []

        if batch:
            await self._upsert(session, entity_type, batch, source_name)
            total += len(batch)

        return total

    @staticmethod
    def _to_row(
        record: Dict[str, Any], source_name: str, entity_type: str
    ) -> Dict[str, Any]:
        """Map a parser record to a database row dict based on entity type."""
        if entity_type == "compound":
            return {
                "source": source_name,
                "external_id": record.get("external_id") or record.get("source_id"),
                "smiles": record.get("smiles"),
                "metadata": record.get("props"),
            }
        if entity_type == "protein":
            return {
                "source": source_name,
                "external_id": record.get("external_id"),
                "sequence": record.get("sequence"),
                "metadata": {"description": record.get("description")},
            }
        raise ValueError(f"Unknown entity_type: '{entity_type}'")

    @staticmethod
    async def _upsert(
        session: AsyncSession,
        entity_type: str,
        batch: List[Dict[str, Any]],
        source_name: str,
    ) -> None:
        """Idempotent bulk upsert via PostgreSQL ON CONFLICT."""
        with tracer.start_as_current_span("ingestion.upsert_batch") as span:
            span.set_attribute("entity_type", entity_type)
            span.set_attribute("batch_size", len(batch))

            BATCH_SIZE.labels(source=source_name, entity_type=entity_type).observe(
                len(batch)
            )

            if entity_type == "compound":
                table = Compound.__table__
                constraint = "uq_compound_source_external_id"
                update_cols = {"smiles": "smiles", "metadata": "metadata"}
            elif entity_type == "protein":
                table = Protein.__table__
                constraint = "uq_protein_source_external_id"
                update_cols = {"sequence": "sequence", "metadata": "metadata"}
            else:
                raise ValueError(f"Unknown entity_type: '{entity_type}'")

            stmt = insert(table).values(batch)
            stmt = stmt.on_conflict_do_update(
                constraint=constraint,
                set_={col: getattr(stmt.excluded, col) for col in update_cols},
            )
            await session.execute(stmt)
            await session.commit()

            RECORDS_UPSERTED.labels(
                source=source_name, entity_type=entity_type
            ).inc(len(batch))
