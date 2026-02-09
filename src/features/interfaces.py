from abc import ABC, abstractmethod
from typing import Dict, Any
import pandas as pd

class FeatureGenerator(ABC):
    """
    Abstract strategy for featurizing raw molecular or protein data.
    """

    @abstractmethod
    def fit(self, df: pd.DataFrame) -> None:
        """
        Fit the featurizer (e.g., learn vocabulary) if stateful.
        Optional for stateless featurizers (like Morgan Fingerprints).
        """
        pass

    @abstractmethod
    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Generate features from the raw input columns.
        """
        pass
    
    @abstractmethod
    def feature_manifest(self) -> Dict[str, Any]:
        """
        Returns metadata about the feature generation strategy.
        MUST pin dependency versions (e.g., RDKit), parameters (radius, nBits), and embedding model/version.
        """
        pass
