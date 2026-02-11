"""
API integration tests — verifies all FastAPI endpoints via TestClient.

These tests do NOT require a running database or S3. The /health endpoint
checks are patched to isolate from infrastructure.
"""
import pytest
from unittest.mock import patch, AsyncMock, MagicMock
from fastapi.testclient import TestClient


@pytest.fixture
def client():
    """Fresh TestClient for each test (resets any app-level state)."""
    from src.main import app
    return TestClient(app)


# ---------------------------------------------------------------------------
# GET /
# ---------------------------------------------------------------------------


def test_root_returns_welcome(client):
    """GET / → 200 with welcome message."""
    resp = client.get("/")
    assert resp.status_code == 200
    body = resp.json()
    assert "message" in body
    assert "Biological Molecular Interaction Intelligence" in body["message"]


# ---------------------------------------------------------------------------
# GET /health
# ---------------------------------------------------------------------------


def test_health_returns_ok_when_all_components_up(client):
    """GET /health → 200 when DB and S3 are reachable."""
    mock_session = AsyncMock()
    mock_session.__aenter__ = AsyncMock(return_value=mock_session)
    mock_session.__aexit__ = AsyncMock(return_value=False)
    mock_session.execute = AsyncMock()

    mock_storage_instance = MagicMock()
    mock_storage_instance.client.head_bucket = MagicMock()

    with patch("src.db.database.AsyncSessionLocal", return_value=mock_session), \
         patch("src.ingestion.storage.S3StorageProvider", return_value=mock_storage_instance), \
         patch("src.ingestion.storage.settings"):
        resp = client.get("/health")

    assert resp.status_code == 200
    body = resp.json()
    assert body["status"] == "healthy"
    assert body["checks"]["database"] == "ok"
    assert body["checks"]["object_store"] == "ok"
    assert isinstance(body["checks"]["model_loaded"], bool)


def test_health_returns_503_when_db_down(client):
    """GET /health → 503 when database is unreachable."""
    mock_session = AsyncMock()
    mock_session.__aenter__ = AsyncMock(return_value=mock_session)
    mock_session.__aexit__ = AsyncMock(return_value=False)
    mock_session.execute = AsyncMock(side_effect=ConnectionRefusedError("connection refused"))

    mock_storage_instance = MagicMock()
    mock_storage_instance.client.head_bucket = MagicMock()

    with patch("src.db.database.AsyncSessionLocal", return_value=mock_session), \
         patch("src.ingestion.storage.S3StorageProvider", return_value=mock_storage_instance), \
         patch("src.ingestion.storage.settings"):
        resp = client.get("/health")

    assert resp.status_code == 503
    body = resp.json()
    assert body["status"] == "degraded"
    assert "refused" in body["checks"]["database"].lower()


# ---------------------------------------------------------------------------
# GET /metrics
# ---------------------------------------------------------------------------


def test_metrics_returns_prometheus(client):
    """GET /metrics → 200 with Prometheus text format."""
    resp = client.get("/metrics")
    assert resp.status_code == 200
    assert "text/plain" in resp.headers["content-type"] or "text/plain" in resp.headers.get("content-type", "")
    # Prometheus output always contains HELP or TYPE lines
    assert "# " in resp.text or "python_" in resp.text


# ---------------------------------------------------------------------------
# POST /rank
# ---------------------------------------------------------------------------


def test_rank_invalid_smiles_returns_422(client):
    """POST /rank with garbage SMILES → 422."""
    # Patch _get_inference_service so it actually loads InferenceService
    # but InferenceService.rank will raise ValueError for bad SMILES
    mock_service = MagicMock()
    mock_service.rank.side_effect = ValueError("Invalid SMILES string: NOT_A_MOLECULE")

    with patch("src.main._get_inference_service", return_value=mock_service):
        resp = client.post("/rank", json={"smiles": "NOT_A_MOLECULE", "top_k": 5})

    assert resp.status_code == 422


def test_rank_success_with_mock(client):
    """POST /rank with mocked inference → 200 with ranked results."""
    mock_service = MagicMock()
    mock_service.model_version = "v1"
    mock_service.rank.return_value = [
        {"rank": 1, "external_id": "P31749", "source": "uniprot", "score": 0.943},
        {"rank": 2, "external_id": "P00533", "source": "uniprot", "score": 0.871},
    ]

    with patch("src.main._get_inference_service", return_value=mock_service):
        resp = client.post("/rank", json={"smiles": "CC(=O)Oc1ccccc1C(=O)O", "top_k": 2})

    assert resp.status_code == 200
    body = resp.json()
    assert body["smiles"] == "CC(=O)Oc1ccccc1C(=O)O"
    assert body["top_k"] == 2
    assert body["model_version"] == "v1"
    assert len(body["results"]) == 2
    assert body["results"][0]["rank"] == 1
    assert body["results"][0]["score"] == 0.943


# ---------------------------------------------------------------------------
# OpenAPI Schema
# ---------------------------------------------------------------------------


def test_openapi_schema_accessible(client):
    """GET /openapi.json → 200 with valid schema."""
    resp = client.get("/openapi.json")
    assert resp.status_code == 200
    schema = resp.json()
    assert schema["info"]["title"] == "Biological Molecular Interaction Intelligence"
    assert "/rank" in schema["paths"]
    assert "/health" in schema["paths"]
