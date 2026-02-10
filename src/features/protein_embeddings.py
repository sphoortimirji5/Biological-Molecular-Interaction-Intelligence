"""
ProteinEmbeddingGenerator — ESM-2 mean-pooled embeddings via HuggingFace.

Converts protein amino-acid sequences into fixed-length dense vectors
using Facebook's ESM-2 language model. Supports batched inference
and configurable model selection.
"""
import structlog
import torch
import numpy as np
import pandas as pd
from typing import Dict, Any, List
from datetime import datetime, timezone

from transformers import AutoTokenizer, AutoModel

from src.features.interfaces import FeatureGenerator
from src.config import settings

logger = structlog.get_logger(__name__)


class ProteinEmbeddingGenerator(FeatureGenerator):
    """Generate ESM-2 embeddings for protein sequences."""

    def __init__(
        self,
        model_name: str | None = None,
        batch_size: int = 8,
        max_length: int = 1024,
    ):
        self.model_name = model_name or settings.esm2_model_name
        self.batch_size = batch_size
        self.max_length = max_length
        self._dim: int | None = None

        # Lazy-load model (expensive)
        self._tokenizer = None
        self._model = None

    def _load_model(self) -> None:
        """Load tokenizer and model on first use."""
        if self._model is not None:
            return

        logger.info("loading_esm2_model", model=self.model_name)
        self._tokenizer = AutoTokenizer.from_pretrained(self.model_name)
        self._model = AutoModel.from_pretrained(self.model_name)
        self._model.eval()

        # Detect embedding dimension from model config
        self._dim = self._model.config.hidden_size
        logger.info("esm2_model_loaded", model=self.model_name, dim=self._dim)

    def fit(self, df: pd.DataFrame) -> None:
        """Pre-load the model (optional — transform auto-loads)."""
        self._load_model()

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Generate embeddings from a DataFrame with `external_id`, `source`,
        and `sequence` columns.

        Returns a DataFrame with `external_id`, `source`, and `emb_0`..`emb_{dim-1}`.
        """
        self._load_model()
        rows: List[Dict[str, Any]] = []

        sequences = df.to_dict("records")

        for batch_start in range(0, len(sequences), self.batch_size):
            batch = sequences[batch_start : batch_start + self.batch_size]
            batch_seqs = [r["sequence"][:self.max_length] for r in batch]

            embeddings = self._embed_batch(batch_seqs)

            for record, emb in zip(batch, embeddings):
                rows.append(
                    {
                        "external_id": record["external_id"],
                        "source": record["source"],
                        **{f"emb_{i}": float(emb[i]) for i in range(len(emb))},
                    }
                )

        logger.info(
            "protein_embeddings_generated",
            count=len(rows),
            model=self.model_name,
            dim=self._dim,
        )
        return pd.DataFrame(rows)

    def _embed_batch(self, sequences: List[str]) -> np.ndarray:
        """Tokenize, run through ESM-2, and mean-pool per-residue outputs."""
        inputs = self._tokenizer(
            sequences,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=self.max_length,
        )

        with torch.no_grad():
            outputs = self._model(**inputs)

        # Mean-pool over sequence length (dim=1), ignoring padding
        attention_mask = inputs["attention_mask"].unsqueeze(-1)
        token_embeddings = outputs.last_hidden_state
        masked = token_embeddings * attention_mask
        summed = masked.sum(dim=1)
        counts = attention_mask.sum(dim=1)
        mean_pooled = summed / counts

        return mean_pooled.numpy()

    def feature_manifest(self) -> Dict[str, Any]:
        """Reproducibility metadata — pins model, dimension, and pooling."""
        self._load_model()
        return {
            "type": "esm2",
            "model": self.model_name,
            "dim": self._dim,
            "pooling": "mean",
            "max_seq_length": self.max_length,
            "generated_at": datetime.now(timezone.utc).isoformat(),
        }
