import pytest
import pytest_asyncio
import asyncio
from typing import Generator
from sqlalchemy import text
from src.db.database import AsyncSessionLocal


@pytest.fixture(scope="session")
def event_loop() -> Generator[asyncio.AbstractEventLoop, None, None]:
    """
    Session-scoped event loop.
    Required because SQLAlchemy's async engine binds to the loop at creation time.
    """
    loop = asyncio.get_event_loop_policy().new_event_loop()
    yield loop
    loop.close()


@pytest_asyncio.fixture(autouse=False)
async def db_cleanup(request):
    """
    Per-test DB isolation via TRUNCATE.

    Only runs on tests marked with @pytest.mark.integration.
    Unit tests using mocks skip this entirely.
    """
    marker = request.node.get_closest_marker("integration")
    if marker is None:
        yield
        return

    async with AsyncSessionLocal() as session:
        await session.execute(
            text("TRUNCATE TABLE interactions, compounds, proteins, ingestion_runs RESTART IDENTITY CASCADE")
        )
        await session.commit()
    yield
