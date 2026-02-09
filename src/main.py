from fastapi import FastAPI, Response
from prometheus_client import generate_latest, CONTENT_TYPE_LATEST
from src.logging_config import get_logger

logger = get_logger(__name__)

app = FastAPI(title="Biological Molecular Interaction Intelligence", version="0.1.0")


@app.get("/")
async def root():
    return {"message": "Welcome to the Biological Molecular Interaction Intelligence API"}


@app.get("/health")
async def health():
    return {"status": "ok"}


@app.get("/metrics")
async def metrics():
    """Prometheus-compatible metrics endpoint."""
    return Response(
        content=generate_latest(),
        media_type=CONTENT_TYPE_LATEST,
    )
