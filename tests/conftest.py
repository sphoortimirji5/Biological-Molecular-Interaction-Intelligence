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


@pytest_asyncio.fixture(autouse=True)
async def db_cleanup():
    """
    Per-test DB isolation via TRUNCATE before each test.
    
    Uses @pytest_asyncio.fixture (not @pytest.fixture) so the async
    fixture is actually awaited by pytest-asyncio 0.23+.
    """
    async with AsyncSessionLocal() as session:
        await session.execute(
            text("TRUNCATE TABLE interactions, compounds, proteins, ingestion_runs RESTART IDENTITY CASCADE")
        )
        await session.commit()
    yield
