"""
Production middleware for the Inference API.

- **APIKeyMiddleware**: Validates ``X-API-Key`` header.  Disabled when
  ``api_key`` setting is empty (local dev).
- **RequestLoggingMiddleware**: Emits structured logs and records a
  Prometheus histogram for every request.
"""
import time
import structlog
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request
from starlette.responses import JSONResponse
from prometheus_client import Histogram

logger = structlog.get_logger(__name__)

# Prometheus latency histogram
HTTP_REQUEST_DURATION = Histogram(
    "http_request_duration_seconds",
    "HTTP request latency in seconds",
    labelnames=["method", "path", "status"],
    buckets=[0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0],
)

# Paths that bypass API key authentication
PUBLIC_PATHS = frozenset({
    "/", "/health", "/metrics", "/docs", "/redoc", "/openapi.json",
})


class APIKeyMiddleware(BaseHTTPMiddleware):
    """
    Reject requests without a valid ``X-API-Key`` header.

    Behaviour:
    - If ``api_key`` is falsy (empty string / None), authentication is
      **disabled** — every request passes through.
    - Public paths (``/health``, ``/docs``, etc.) always bypass auth.
    """

    def __init__(self, app, api_key: str = ""):
        super().__init__(app)
        self.api_key = api_key

    async def dispatch(self, request: Request, call_next):
        # Auth disabled when no key configured
        if not self.api_key:
            return await call_next(request)

        # Public paths skip auth
        if request.url.path in PUBLIC_PATHS:
            return await call_next(request)

        # Validate header
        provided = request.headers.get("X-API-Key", "")
        if provided != self.api_key:
            logger.warning(
                "auth_rejected",
                path=request.url.path,
                client=request.client.host if request.client else "unknown",
            )
            return JSONResponse(
                status_code=401,
                content={"detail": "Invalid or missing API key"},
            )

        return await call_next(request)


class RequestLoggingMiddleware(BaseHTTPMiddleware):
    """
    Log every request/response and record latency in Prometheus.

    Emits a structlog event with method, path, status code, and duration.
    """

    async def dispatch(self, request: Request, call_next):
        start = time.perf_counter()

        response = await call_next(request)

        duration = time.perf_counter() - start
        path = request.url.path
        method = request.method
        status = response.status_code

        # Prometheus histogram
        HTTP_REQUEST_DURATION.labels(
            method=method,
            path=path,
            status=str(status),
        ).observe(duration)

        # Structured log (skip noisy /metrics polling)
        if path != "/metrics":
            logger.info(
                "http_request",
                method=method,
                path=path,
                status=status,
                duration_ms=round(duration * 1000, 2),
                client=request.client.host if request.client else "unknown",
            )

        return response
