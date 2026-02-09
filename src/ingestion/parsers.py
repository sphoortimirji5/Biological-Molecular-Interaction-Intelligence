"""
Format-specific parsers for ingested data files.

Each parser yields normalised dicts. The manager selects the parser
based on the ``format`` field returned by ``DataSource.discover()``.

Supported formats: ``sdf``, ``fasta``, ``csv``.
"""
import gzip
from pathlib import Path
from typing import Generator, Dict, Any
import pandas as pd
from Bio import SeqIO
from rdkit import Chem
from src.logging_config import get_logger

logger = get_logger(__name__)


class FileParser:
    """Base class providing transparent gzip handling for format-specific parsers."""

    @staticmethod
    def _open(file_path: Path):
        """Open *file_path*, decompressing gzip on the fly when detected."""
        if file_path.suffix == ".gz":
            return gzip.open(file_path, "rt")
        return open(file_path, "r")


class SDFParser(FileParser):
    """
    Stream-parse SDF/SD files, yielding normalised compound dicts.

    The ``external_id_prop`` parameter controls which SDF property
    is extracted as the external identifier. Defaults to ``_Name``
    (universal SDF convention). Override for sources that use a
    custom property (e.g. ``chembl_id``).
    """

    def __init__(self, external_id_prop: str = "_Name") -> None:
        self._id_prop = external_id_prop

    def parse(self, file_path: Path) -> Generator[Dict[str, Any], None, None]:
        if file_path.suffix == ".gz":
            suppl = Chem.ForwardSDMolSupplier(gzip.open(file_path))
        else:
            suppl = Chem.ForwardSDMolSupplier(str(file_path))

        for mol in suppl:
            if mol is None:
                continue
            try:
                yield {
                    "external_id": (
                        mol.GetProp(self._id_prop)
                        if mol.HasProp(self._id_prop)
                        else None
                    ),
                    "smiles": Chem.MolToSmiles(mol),
                    "props": mol.GetPropsAsDict(),
                }
            except Exception as e:
                logger.warning("Failed to parse molecule in %s: %s", file_path, e)


class FastaParser(FileParser):
    """Stream-parse FASTA files, yielding normalised protein dicts."""

    @staticmethod
    def parse(file_path: Path) -> Generator[Dict[str, Any], None, None]:
        with FileParser._open(file_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                yield {
                    "external_id": record.id,
                    "description": record.description,
                    "sequence": str(record.seq),
                }


class CSVParser(FileParser):
    """Chunked CSV reader for tabular interaction data."""

    @staticmethod
    def parse(
        file_path: Path, chunk_size: int = 10000
    ) -> Generator[pd.DataFrame, None, None]:
        """Yield fixed-size DataFrame chunks for memory-bounded processing."""
        for chunk in pd.read_csv(file_path, chunksize=chunk_size):
            yield chunk


# ---------------------------------------------------------------------------
# Format → Parser registry
# ---------------------------------------------------------------------------

PARSER_REGISTRY: Dict[str, type] = {
    "sdf": SDFParser,
    "fasta": FastaParser,
    "csv": CSVParser,
}


def get_parser(fmt: str, **kwargs) -> FileParser:
    """Resolve a parser instance by format string. Raises ValueError if unknown."""
    cls = PARSER_REGISTRY.get(fmt.lower())
    if cls is None:
        raise ValueError(
            f"Unknown format: '{fmt}'. Supported: {list(PARSER_REGISTRY.keys())}"
        )
    return cls(**kwargs)
