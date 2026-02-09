from abc import ABC, abstractmethod
from typing import Dict, Any
from pathlib import Path


class DataSource(ABC):
    """
    Abstract contract for data sources.

    Each implementation provides:
        - ``discover()`` — metadata including ``format`` and ``entity_type``
        - ``fetch()`` — download raw artifact to a local path

    The manager uses ``format`` to select the parser and ``entity_type``
    to route records to the correct database table.
    """

    @abstractmethod
    def discover(self) -> Dict[str, Any]:
        """
        Return source metadata. Must include at minimum:

        - ``source``: unique label (e.g. ``"chembl"``, ``"pharma_x"``)
        - ``version``: dataset version string
        - ``format``: file format — ``"sdf"`` | ``"fasta"`` | ``"csv"``
        - ``entity_type``: target entity — ``"compound"`` | ``"protein"`` | ``"interaction"``
        """

    @abstractmethod
    def fetch(self, destination: Path) -> Path:
        """
        Download the raw data artifact to *destination*.

        Returns:
            Path to the downloaded file.
        """
