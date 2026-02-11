"""
Unit tests for production middleware — API key auth and request logging.
"""
import pytest
from unittest.mock import patch, MagicMock
from fastapi import FastAPI
from fastapi.testclient import TestClient
from starlette.responses import JSONResponse

from src.middleware import APIKeyMiddleware, RequestLoggingMiddleware


# ---------------------------------------------------------------------------
# Helpers — build minimal apps with middleware for isolated testing
# ---------------------------------------------------------------------------


def _make_app(api_key: str = "") -> FastAPI:
    """Build a tiny FastAPI app with our middleware stack."""
    app = FastAPI()

    app.add_middleware(RequestLoggingMiddleware)
    app.add_middleware(APIKeyMiddleware, api_key=api_key)

    @app.get("/")
    async def root():
        return {"message": "ok"}

    @app.get("/health")
    async def health():
        return {"status": "healthy"}

    @app.post("/rank")
    async def rank():
        return {"results": []}

    return app


# ---------------------------------------------------------------------------
# APIKeyMiddleware
# ---------------------------------------------------------------------------


class TestAPIKeyMiddleware:
    """Validates X-API-Key header behaviour."""

    def test_no_api_key_required_when_disabled(self):
        """Empty api_key config → all requests pass without a header."""
        client = TestClient(_make_app(api_key=""))
        resp = client.get("/")
        assert resp.status_code == 200

        resp = client.post("/rank")
        assert resp.status_code == 200

    def test_valid_api_key_passes(self):
        """Correct X-API-Key header → 200."""
        client = TestClient(_make_app(api_key="secret-123"))
        resp = client.post("/rank", headers={"X-API-Key": "secret-123"})
        assert resp.status_code == 200

    def test_invalid_api_key_rejected(self):
        """Wrong X-API-Key header → 401."""
        client = TestClient(_make_app(api_key="secret-123"))
        resp = client.post("/rank", headers={"X-API-Key": "wrong-key"})
        assert resp.status_code == 401
        assert "Invalid" in resp.json()["detail"]

    def test_missing_api_key_rejected(self):
        """No X-API-Key header → 401."""
        client = TestClient(_make_app(api_key="secret-123"))
        resp = client.post("/rank")
        assert resp.status_code == 401

    def test_public_paths_skip_auth(self):
        """Health and root bypass authentication."""
        client = TestClient(_make_app(api_key="secret-123"))

        resp = client.get("/")
        assert resp.status_code == 200

        resp = client.get("/health")
        assert resp.status_code == 200


# ---------------------------------------------------------------------------
# RequestLoggingMiddleware
# ---------------------------------------------------------------------------


class TestRequestLoggingMiddleware:
    """Validates request logging and Prometheus histogram recording."""

    def test_request_logging_emits_log(self):
        """Verify structlog emits an http_request event."""
        client = TestClient(_make_app())

        with patch("src.middleware.logger") as mock_logger:
            client.get("/")

        mock_logger.info.assert_called()
        call_args = mock_logger.info.call_args
        assert call_args[0][0] == "http_request"

    def test_prometheus_histogram_recorded(self):
        """Verify the histogram observe() is called."""
        client = TestClient(_make_app())

        with patch("src.middleware.HTTP_REQUEST_DURATION") as mock_hist:
            mock_labels = MagicMock()
            mock_hist.labels.return_value = mock_labels
            client.get("/")

        mock_hist.labels.assert_called()
        mock_labels.observe.assert_called_once()
