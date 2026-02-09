from abc import ABC, abstractmethod
from typing import Dict, Any
import pandas as pd
from pathlib import Path

class ModelWrapper(ABC):
    """
    Standard contract for ML model interaction.
    Enforces Pandas DataFrame I/O for consistency across researchers.
    """

    @abstractmethod
    def train(self, df: pd.DataFrame) -> None:
        """
        Train the underlying model using the provided DataFrame.
        """
        pass

    @abstractmethod
    def predict(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Run inference on the provided DataFrame.
        Returns a DataFrame with predictions and confidence scores.
        """
        pass

    @abstractmethod
    def save(self, path: Path) -> None:
        """
        Serialize the model to disk.
        """
        pass

    @abstractmethod
    def load(self, path: Path) -> 'ModelWrapper':
        """
        Deserialize the model from disk.
        """
        pass
    
    @abstractmethod
    def model_signature(self) -> Dict[str, Any]:
        """
        Defines required feature names, dtypes, and feature versions.
        Used to validate inference/training inputs before execution.
        """
        pass
