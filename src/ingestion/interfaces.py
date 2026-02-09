from abc import ABC, abstractmethod
from typing import Dict, Any
from pathlib import Path

class DataSource(ABC):
    """
    Abstract contract for fetching and parsing data from various sources (ChEMBL, UniProt, Local Files).
    """

    @abstractmethod
    def discover(self) -> Dict[str, Any]:
        """
        Returns dataset version, estimated entry stats, and checksums if available.
        Ensures we know what we are about to fetch.
        """
        pass

    @abstractmethod
    def fetch(self, destination: Path) -> Path:
        """
        Fetches the raw data artifact to the specified location.
        MUST be idempotent and content-addressed.
        
        Args:
            destination: The local path to save the artifact.
            
        Returns:
            Publisher path to the downloaded file.
        """
        pass
