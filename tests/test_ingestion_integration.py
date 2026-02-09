"""
Integration tests for the Ingestion Pipeline.

Covers:
- Happy path (SDF compounds, FASTA proteins)
- Idempotency (ingest twice, assert no duplicates)
- Failure path (fetch throws, run marked FAILED)

Sources and storage are injected via constructor (DI), not patched.
"""
import pytest
import gzip
from unittest.mock import MagicMock, patch
from pathlib import Path
from sqlalchemy import select, func
from src.ingestion.manager import IngestionManager, register_source
from src.ingestion.sources import HTTPSource
from src.ingestion.storage_interface import ObjectStorageProvider
from src.db.models import Compound, Protein, IngestionRun
from src.db.database import AsyncSessionLocal

from rdkit import Chem
from io import StringIO

# ---------------------------------------------------------------------------
# Sample Data
# ---------------------------------------------------------------------------

_mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
_mol.SetProp("_Name", "ASPIRIN")
_mol.SetProp("chembl_id", "CHEMBL113")
_sio = StringIO()
_w = Chem.SDWriter(_sio)
_w.write(_mol)
_w.close()
SAMPLE_SDF_BYTES = _sio.getvalue().encode("utf-8")

SAMPLE_FASTA_BYTES = b""">P12345 Test Protein
MKWVTFISLLFLFSSAYSRGVFRR
"""

# ---------------------------------------------------------------------------
# Discovery metadata — now includes format and entity_type
# ---------------------------------------------------------------------------

COMPOUND_DISCOVERY = {
    "source": "test_compound_src",
    "version": "vTEST",
    "format": "sdf",
    "entity_type": "compound",
    "stats": {"estimated": 1},
}

PROTEIN_DISCOVERY = {
    "source": "test_protein_src",
    "version": "vTEST",
    "format": "fasta",
    "entity_type": "protein",
    "stats": {"estimated": 1},
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_gz(dest: Path, data: bytes) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(dest, "wb") as f:
        f.write(data)
    return dest


def _make_mock_storage(sha256: str = "fakehash") -> MagicMock:
    """Create a mock ObjectStorageProvider with preconfigured upload_raw."""
    mock = MagicMock(spec=ObjectStorageProvider)
    mock.upload_raw.return_value = {"sha256": sha256, "bucket": "raw", "key": "test"}
    return mock


def _register_test_source(
    name: str, discovery: dict, data: bytes
) -> None:
    """Register a test source that returns canned data."""
    src = MagicMock(spec=HTTPSource)
    src.discover.return_value = discovery

    def _mock_fetch(dest: Path) -> Path:
        return _write_gz(dest, data)

    src.fetch.side_effect = _mock_fetch
    register_source(name, src)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_compound_ingestion_flow():
    """Happy path: SDF format → compounds table."""
    _register_test_source("test_compound_src", COMPOUND_DISCOVERY, SAMPLE_SDF_BYTES)
    mock_storage = _make_mock_storage("fakehash")

    manager = IngestionManager(storage=mock_storage)
    await manager.run_pipeline("test_compound_src")

    mock_storage.upload_raw.assert_called_once()

    async with AsyncSessionLocal() as session:
        result = await session.execute(
            select(IngestionRun).where(IngestionRun.source == "test_compound_src")
        )
        run = result.scalars().first()
        assert run is not None, "No ingestion_run row"
        assert run.status == "COMPLETED"
        assert run.checksums["sha256"] == "fakehash"

        result = await session.execute(
            select(Compound).where(Compound.source == "test_compound_src")
        )
        compound = result.scalars().first()
        assert compound is not None, "Compound not found"
        assert compound.smiles == "CC(=O)Oc1ccccc1C(=O)O"
        assert compound.external_id == "ASPIRIN"


@pytest.mark.asyncio
async def test_protein_ingestion_flow():
    """Happy path: FASTA format → proteins table."""
    _register_test_source("test_protein_src", PROTEIN_DISCOVERY, SAMPLE_FASTA_BYTES)
    mock_storage = _make_mock_storage("fakehash_prot")

    manager = IngestionManager(storage=mock_storage)
    await manager.run_pipeline("test_protein_src")

    async with AsyncSessionLocal() as session:
        result = await session.execute(
            select(Protein).where(Protein.source == "test_protein_src")
        )
        protein = result.scalars().first()
        assert protein is not None, "Protein not found"
        assert protein.sequence == "MKWVTFISLLFLFSSAYSRGVFRR"
        assert protein.external_id == "P12345"


@pytest.mark.asyncio
async def test_idempotency_no_duplicates():
    """Ingest the same source twice → exactly 1 compound row, not 2."""
    _register_test_source("test_compound_src", COMPOUND_DISCOVERY, SAMPLE_SDF_BYTES)

    for run_num in range(2):
        mock_storage = _make_mock_storage(f"hash_{run_num}")
        manager = IngestionManager(storage=mock_storage)
        await manager.run_pipeline("test_compound_src")

    async with AsyncSessionLocal() as session:
        result = await session.execute(
            select(func.count()).select_from(Compound).where(
                Compound.source == "test_compound_src"
            )
        )
        count = result.scalar()
        assert count == 1, f"Expected 1 compound, got {count} — upsert is broken"

        result = await session.execute(
            select(func.count()).select_from(IngestionRun).where(
                IngestionRun.source == "test_compound_src"
            )
        )
        run_count = result.scalar()
        assert run_count == 2, f"Expected 2 runs, got {run_count}"


@pytest.mark.asyncio
async def test_failure_path_fetch_throws():
    """When fetch raises, the run should be marked FAILED with error details."""
    failing_src = MagicMock(spec=HTTPSource)
    failing_src.discover.return_value = COMPOUND_DISCOVERY

    def _exploding_fetch(dest: Path) -> Path:
        raise ConnectionError("Simulated network failure")

    failing_src.fetch.side_effect = _exploding_fetch
    register_source("test_fail_src", failing_src)

    mock_storage = _make_mock_storage()
    manager = IngestionManager(storage=mock_storage)

    with pytest.raises(ConnectionError):
        await manager.run_pipeline("test_fail_src")

    async with AsyncSessionLocal() as session:
        result = await session.execute(
            select(IngestionRun).where(IngestionRun.source == "test_fail_src")
        )
        run = result.scalars().first()
        assert run is not None, "No ingestion_run row after failure"
        assert run.status == "FAILED"
        assert run.ended_at is not None
        assert "Simulated network failure" in str(run.stats.get("error", ""))
