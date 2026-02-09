"""
OpenTelemetry tracing configuration.

Exporter is selected via the ``OTEL_EXPORTER`` environment variable:

- ``console``  — prints spans to stdout (default for development)
- ``otlp``     — sends to any OTLP-compatible backend (Jaeger, Datadog, Grafana Cloud)
- ``none``     — tracing disabled

Additional env vars for OTLP:
- ``OTEL_EXPORTER_OTLP_ENDPOINT`` — e.g. ``http://jaeger:4317``
- ``OTEL_SERVICE_NAME``           — defaults to ``bio-interaction-intelligence``

Usage::

    from src.tracing import get_tracer
    tracer = get_tracer(__name__)

    with tracer.start_as_current_span("my_operation") as span:
        span.set_attribute("source", "chembl")
        ...
"""
import os
from opentelemetry import trace
from opentelemetry.sdk.trace import TracerProvider
from opentelemetry.sdk.resources import Resource


def _init_tracing() -> None:
    """Initialize the tracer provider with the configured exporter."""
    exporter_type = os.environ.get("OTEL_EXPORTER", "console").lower()
    service_name = os.environ.get("OTEL_SERVICE_NAME", "bio-interaction-intelligence")

    resource = Resource.create({"service.name": service_name})
    provider = TracerProvider(resource=resource)

    if exporter_type == "none":
        pass  # No exporter — tracing disabled
    elif exporter_type == "otlp":
        from opentelemetry.exporter.otlp.proto.grpc.trace_exporter import (
            OTLPSpanExporter,
        )
        from opentelemetry.sdk.trace.export import BatchSpanProcessor

        endpoint = os.environ.get("OTEL_EXPORTER_OTLP_ENDPOINT", "http://localhost:4317")
        exporter = OTLPSpanExporter(endpoint=endpoint)
        provider.add_span_processor(BatchSpanProcessor(exporter))
    else:
        # Default: console exporter for local dev
        from opentelemetry.sdk.trace.export import (
            SimpleSpanProcessor,
            ConsoleSpanExporter,
        )

        provider.add_span_processor(SimpleSpanProcessor(ConsoleSpanExporter()))

    trace.set_tracer_provider(provider)


def get_tracer(name: str) -> trace.Tracer:
    """Return a tracer for *name*."""
    return trace.get_tracer(name)


# Auto-initialize on import
_init_tracing()
