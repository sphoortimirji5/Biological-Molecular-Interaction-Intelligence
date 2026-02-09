import gzip
from pathlib import Path
from typing import Generator, Dict, Any, Union
import pandas as pd
from Bio import SeqIO
from rdkit import Chem
import logging

logger = logging.getLogger(__name__)

class FileParser:
    """Base class for file parsers."""
    
    @staticmethod
    def _open(file_path: Path):
        """Opens file, handling gzip if necessary."""
        if file_path.suffix == '.gz':
            return gzip.open(file_path, 'rt')
        return open(file_path, 'r')

class SDFParser(FileParser):
    """Parses .sdf files for chemical compounds."""
    
    @staticmethod
    def parse(file_path: Path) -> Generator[Dict[str, Any], None, None]:
        if file_path.suffix == '.gz':
             suppl = Chem.ForwardSDMolSupplier(gzip.open(file_path))
        else:
             suppl = Chem.ForwardSDMolSupplier(str(file_path))
             
        for mol in suppl:
            if mol is None:
                continue
            
            try:
                yield {
                    "source_id": mol.GetProp("_Name") if mol.HasProp("_Name") else None,
                    "external_id": mol.GetProp("chembl_id") if mol.HasProp("chembl_id") else None,
                    "smiles": Chem.MolToSmiles(mol),
                    "props": mol.GetPropsAsDict()
                }
            except Exception as e:
                logger.warning(f"Failed to parse molecule in {file_path}: {e}")

class FastaParser(FileParser):
    """Parses .fasta files for protein sequences."""
    
    @staticmethod
    def parse(file_path: Path) -> Generator[Dict[str, Any], None, None]:
        with FileParser._open(file_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                yield {
                    "external_id": record.id,
                    "description": record.description,
                    "sequence": str(record.seq)
                }

class CSVParser(FileParser):
    """Parses .csv files using Pandas (chunked)."""
    
    @staticmethod
    def parse(file_path: Path, chunk_size: int = 10000) -> Generator[pd.DataFrame, None, None]:
        """Yields chunks of DataFrames."""
        for chunk in pd.read_csv(file_path, chunksize=chunk_size):
            yield chunk
