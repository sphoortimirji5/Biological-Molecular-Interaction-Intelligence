"""
Pydantic schemas for the Inference API.

Provides request/response models with OpenAPI examples for Swagger UI.
"""
from pydantic import BaseModel, Field
from typing import List


class RankRequest(BaseModel):
    """Input: a SMILES string and desired number of top targets."""

    smiles: str = Field(
        ...,
        description="SMILES string of the query compound",
        json_schema_extra={"examples": ["CC(=O)Oc1ccccc1C(=O)O"]},
    )
    top_k: int = Field(
        default=10,
        ge=1,
        le=500,
        description="Number of top-ranked targets to return",
    )

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                    "top_k": 5,
                },
                {
                    "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                    "top_k": 10,
                },
            ]
        }
    }


class RankedTarget(BaseModel):
    """A single ranked protein target with its interaction probability."""

    rank: int = Field(..., description="1-indexed rank position", examples=[1])
    external_id: str = Field(
        ...,
        description="UniProt accession or compound external ID",
        examples=["P12345"],
    )
    source: str = Field(
        ...,
        description="Data source the protein was ingested from",
        examples=["uniprot"],
    )
    score: float = Field(
        ...,
        description="Interaction probability (0–1 for learned, cosine for similarity)",
        examples=[0.943],
    )


class RankResponse(BaseModel):
    """Output: ranked list of predicted protein targets."""

    smiles: str = Field(
        ...,
        description="Echo of the input SMILES query",
        examples=["CC(=O)Oc1ccccc1C(=O)O"],
    )
    top_k: int = Field(..., description="Number of results requested", examples=[5])
    model_version: str = Field(
        ...,
        description="Version of the scoring model used",
        examples=["v1"],
    )
    results: List[RankedTarget] = Field(
        ..., description="Ranked protein targets, highest score first"
    )

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                    "top_k": 3,
                    "model_version": "v1",
                    "results": [
                        {"rank": 1, "external_id": "P31749", "source": "uniprot", "score": 0.943},
                        {"rank": 2, "external_id": "P00533", "source": "uniprot", "score": 0.871},
                        {"rank": 3, "external_id": "P04637", "source": "uniprot", "score": 0.812},
                    ],
                }
            ]
        }
    }
