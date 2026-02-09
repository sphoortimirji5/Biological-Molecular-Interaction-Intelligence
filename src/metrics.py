"""
Prometheus metrics for the ingestion pipeline.

All metrics are prefixed with ``ingestion_`` and labeled by ``source``
and/or ``entity_type`` for per-source drill-down.

Metrics are exposed via ``/metrics`` in the FastAPI app.
"""
from prometheus_client import Counter, Histogram

# ---------------------------------------------------------------------------
# Counters
# ---------------------------------------------------------------------------

RUNS_TOTAL = Counter(
    "ingestion_runs_total",
    "Total ingestion runs by source and final status",
    ["source", "status"],
)

RECORDS_UPSERTED = Counter(
    "ingestion_records_upserted_total",
    "Total records upserted by source and entity type",
    ["source", "entity_type"],
)

PARSE_ERRORS = Counter(
    "ingestion_parse_errors_total",
    "Parse-level errors by source and format",
    ["source", "format"],
)

# ---------------------------------------------------------------------------
# Histograms (stage latency)
# ---------------------------------------------------------------------------

STAGE_DURATION = Histogram(
    "ingestion_stage_duration_seconds",
    "Duration of each pipeline stage in seconds",
    ["source", "stage"],
    buckets=(0.1, 0.5, 1, 5, 10, 30, 60, 120, 300, 600),
)

BATCH_SIZE = Histogram(
    "ingestion_batch_size",
    "Number of records per upsert batch",
    ["source", "entity_type"],
    buckets=(10, 50, 100, 250, 500, 1000, 2500, 5000),
)
