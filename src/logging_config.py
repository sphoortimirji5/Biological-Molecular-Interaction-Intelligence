"""
Structured logging configuration.

- JSON output in production (APP_ENV=production)
- Colored console output in development
- Automatic context binding (timestamp, level, logger name)

Usage::

    from src.logging_config import get_logger
    logger = get_logger(__name__)
    logger.info("fetch_started", source="chembl", url="...")
"""
import logging
import os
import sys
import structlog


def configure_logging() -> None:
    """Configure structlog processors and stdlib integration."""
    shared_processors: list = [
        structlog.contextvars.merge_contextvars,
        structlog.stdlib.add_log_level,
        structlog.stdlib.add_logger_name,
        structlog.processors.TimeStamper(fmt="iso"),
        structlog.processors.StackInfoRenderer(),
        structlog.processors.format_exc_info,
    ]

    env = os.environ.get("APP_ENV", "development").lower()

    if env == "production":
        renderer = structlog.processors.JSONRenderer()
    else:
        renderer = structlog.dev.ConsoleRenderer()

    structlog.configure(
        processors=[
            *shared_processors,
            structlog.stdlib.ProcessorFormatter.wrap_for_formatter,
        ],
        logger_factory=structlog.stdlib.LoggerFactory(),
        wrapper_class=structlog.stdlib.BoundLogger,
        cache_logger_on_first_use=True,
    )

    formatter = structlog.stdlib.ProcessorFormatter(
        processors=[
            structlog.stdlib.ProcessorFormatter.remove_processors_meta,
            renderer,
        ],
    )

    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(formatter)

    root_logger = logging.getLogger()
    root_logger.handlers.clear()
    root_logger.addHandler(handler)
    root_logger.setLevel(os.environ.get("LOG_LEVEL", "INFO").upper())


def get_logger(name: str) -> structlog.stdlib.BoundLogger:
    """Return a bound structlog logger for *name*."""
    return structlog.get_logger(name)


# Auto-configure on import
configure_logging()
