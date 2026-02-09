"""
Preconfigured DataSource implementations for public databases.

Each source is fully configurable via constructor parameters.
URLs default to environment variables; override via constructor for
proprietary sources.
"""
import logging
import urllib.request
import shutil
from pathlib import Path
from typing import Dict, Any, Optional
from src.ingestion.interfaces import DataSource
from src.config import settings
from src.logging_config import get_logger

logger = get_logger(__name__)


class HTTPSource(DataSource):
    """
    Generic HTTP-based data source.

    Downloads a file from a URL and declares its format and entity type.
    All parameters are set via constructor — no global config dependency.
    """

    def __init__(
        self,
        source: str,
        url: str,
        fmt: str,
        entity_type: str,
        version: str = "latest",
        stats: Optional[Dict[str, Any]] = None,
    ) -> None:
        self._source = source
        self._url = url
        self._format = fmt
        self._entity_type = entity_type
        self._version = version
        self._stats = stats or {}

    def discover(self) -> Dict[str, Any]:
        return {
            "source": self._source,
            "version": self._version,
            "format": self._format,
            "entity_type": self._entity_type,
            "stats": self._stats,
            "url": self._url,
        }

    def fetch(self, destination: Path) -> Path:
        """Stream-download from the configured URL."""
        logger.info("Fetching %s from %s", self._source, self._url)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(self._url) as response, \
             open(destination, "wb") as out_file:
            shutil.copyfileobj(response, out_file)
        logger.info("Downloaded %s to %s", self._source, destination)
        return destination


def chembl_source(
    url: str | None = None,
    version: str = "v33",
) -> HTTPSource:
    """Factory: preconfigured ChEMBL compound source. URL from CHEMBL_SDF_URL env var."""
    return HTTPSource(
        source="chembl",
        url=url or settings.chembl_sdf_url,
        fmt="sdf",
        entity_type="compound",
        version=version,
        stats={"estimated_compounds": 2_400_000},
    )


def uniprot_source(
    url: str | None = None,
    version: str = "2024_01",
) -> HTTPSource:
    """Factory: preconfigured UniProt protein source. URL from UNIPROT_FASTA_URL env var."""
    return HTTPSource(
        source="uniprot",
        url=url or settings.uniprot_fasta_url,
        fmt="fasta",
        entity_type="protein",
        version=version,
        stats={"estimated_proteins": 570_000},
    )
