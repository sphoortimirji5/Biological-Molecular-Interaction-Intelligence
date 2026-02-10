"""
MorganFingerprintGenerator — ECFP4-style circular fingerprints via RDKit.

Stateless featurizer: takes a DataFrame with a `smiles` column and
produces 2048-bit (configurable) binary fingerprint vectors.
"""
import structlog
import numpy as np
import pandas as pd
from typing import Dict, Any
from datetime import datetime, timezone

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

from src.features.interfaces import FeatureGenerator
from src.config import settings

logger = structlog.get_logger(__name__)


class MorganFingerprintGenerator(FeatureGenerator):
    """Generate Morgan (ECFP) fingerprints for small-molecule compounds."""

    def __init__(
        self,
        radius: int | None = None,
        n_bits: int | None = None,
    ):
        self.radius = radius or settings.morgan_radius
        self.n_bits = n_bits or settings.morgan_n_bits
        self._invalid_count = 0

    def fit(self, df: pd.DataFrame) -> None:
        """No-op — Morgan fingerprints are stateless."""
        pass

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Generate fingerprints from a DataFrame with `external_id`, `source`,
        and `smiles` columns.

        Invalid SMILES are skipped and logged.
        Returns a DataFrame with `external_id`, `source`, and `fp_0`..`fp_{n_bits-1}`.
        """
        rows = []
        self._invalid_count = 0

        for _, row in df.iterrows():
            smiles = row.get("smiles", "")
            mol = Chem.MolFromSmiles(smiles) if smiles else None

            if mol is None:
                self._invalid_count += 1
                logger.warning(
                    "invalid_smiles_skipped",
                    external_id=row.get("external_id", "unknown"),
                    smiles=smiles[:50] if smiles else "",
                )
                continue

            fp = AllChem.GetMorganFingerprintAsBitVect(
                mol, radius=self.radius, nBits=self.n_bits
            )
            bits = np.array(fp, dtype=np.uint8)

            rows.append(
                {
                    "external_id": row["external_id"],
                    "source": row["source"],
                    **{f"fp_{i}": bits[i] for i in range(self.n_bits)},
                }
            )

        logger.info(
            "morgan_fingerprints_generated",
            valid=len(rows),
            invalid=self._invalid_count,
            total=len(df),
            radius=self.radius,
            n_bits=self.n_bits,
        )

        return pd.DataFrame(rows)

    def feature_manifest(self) -> Dict[str, Any]:
        """Reproducibility metadata — pins RDKit version and parameters."""
        return {
            "type": "morgan",
            "radius": self.radius,
            "n_bits": self.n_bits,
            "rdkit_version": rdkit.__version__,
            "generated_at": datetime.now(timezone.utc).isoformat(),
        }

    @property
    def invalid_count(self) -> int:
        """Number of invalid SMILES skipped in the last transform call."""
        return self._invalid_count
