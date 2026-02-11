from fastapi import FastAPI, Response, HTTPException
from fastapi.openapi.utils import get_openapi
from prometheus_client import generate_latest, CONTENT_TYPE_LATEST
from src.logging_config import get_logger
from src.inference.schemas import RankRequest, RankResponse, RankedTarget

logger = get_logger(__name__)

app = FastAPI(
    title="Biological Molecular Interaction Intelligence",
    version="0.1.0",
    description=(
        "Predict drug–protein interactions using molecular fingerprints, "
        "protein embeddings, and machine learning. Given a compound SMILES string, "
        "the API ranks all known protein targets by predicted interaction probability."
    ),
    contact={
        "name": "BMII Team",
    },
    license_info={
        "name": "MIT",
    },
)

# Lazy-initialized inference service (loaded on first /rank request)
_inference_service = None


def _get_inference_service():
    """Factory: build InferenceService once from config, cache globally."""
    global _inference_service
    if _inference_service is not None:
        return _inference_service

    from src.config import settings
    from src.features.store import FeatureStore
    from src.ingestion.s3_provider import S3StorageProvider
    from src.models.scoring import get_scoring_strategy
    from src.inference.service import InferenceService

    storage = S3StorageProvider()
    feature_store = FeatureStore(storage=storage, bucket=settings.feature_store_bucket)
    scoring_strategy = get_scoring_strategy()

    _inference_service = InferenceService(
        feature_store=feature_store,
        scoring_strategy=scoring_strategy,
        model_version=settings.model_path.split("/")[-1],
    )
    return _inference_service


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------


@app.get(
    "/",
    tags=["General"],
    summary="API Landing",
    response_description="Welcome message",
    responses={
        200: {
            "content": {
                "application/json": {
                    "example": {
                        "message": "Welcome to the Biological Molecular Interaction Intelligence API"
                    }
                }
            }
        }
    },
)
async def root():
    """Returns a welcome message confirming the API is reachable."""
    return {"message": "Welcome to the Biological Molecular Interaction Intelligence API"}


@app.get(
    "/health",
    tags=["Operations"],
    summary="Health Check",
    response_description="Service health status",
    responses={
        200: {
            "content": {
                "application/json": {
                    "example": {"status": "ok"}
                }
            }
        }
    },
)
async def health():
    """Lightweight liveness probe. Returns `ok` if the service is running."""
    return {"status": "ok"}


@app.get(
    "/metrics",
    tags=["Operations"],
    summary="Prometheus Metrics",
    response_description="Prometheus text-format metrics",
    responses={
        200: {
            "content": {
                "text/plain": {
                    "example": '# HELP python_gc_objects_collected_total ...\npython_gc_objects_collected_total{generation="0"} 350.0'
                }
            }
        }
    },
)
async def metrics():
    """Prometheus-compatible metrics endpoint. Scrape this with your Prometheus server."""
    return Response(
        content=generate_latest(),
        media_type=CONTENT_TYPE_LATEST,
    )


@app.post(
    "/rank",
    response_model=RankResponse,
    tags=["Inference"],
    summary="Rank Protein Targets",
    response_description="Ranked list of predicted protein targets",
    responses={
        200: {
            "description": "Successful ranking",
            "content": {
                "application/json": {
                    "example": {
                        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                        "top_k": 3,
                        "model_version": "v1",
                        "results": [
                            {"rank": 1, "external_id": "P31749", "source": "uniprot", "score": 0.943},
                            {"rank": 2, "external_id": "P00533", "source": "uniprot", "score": 0.871},
                            {"rank": 3, "external_id": "P04637", "source": "uniprot", "score": 0.812},
                        ],
                    }
                }
            },
        },
        422: {
            "description": "Invalid SMILES string or request parameters",
            "content": {
                "application/json": {
                    "example": {"detail": "Invalid SMILES: 'NOT_A_MOLECULE' could not be parsed by RDKit"}
                }
            },
        },
    },
)
async def rank(request: RankRequest):
    """
    Rank protein targets for a query compound.

    Accepts a **SMILES** string representing a small molecule and returns the
    top-K most likely interacting protein targets, ranked by predicted
    interaction probability.

    **How it works:**
    1. The SMILES string is converted to a 2048-bit Morgan fingerprint (ECFP4).
    2. Protein embeddings are loaded from the feature store (ESM-2, 320-dim).
    3. The scoring strategy (cosine similarity or trained XGBoost model) ranks
       all proteins against the compound.
    4. The top-K results are returned, sorted by descending score.

    **Example SMILES:**
    - Aspirin: `CC(=O)Oc1ccccc1C(=O)O`
    - Caffeine: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`
    - Ibuprofen: `CC(C)Cc1ccc(cc1)C(C)C(=O)O`
    """
    try:
        service = _get_inference_service()
        results = service.rank(smiles=request.smiles, top_k=request.top_k)
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))

    return RankResponse(
        smiles=request.smiles,
        top_k=request.top_k,
        model_version=service.model_version,
        results=[RankedTarget(**r) for r in results],
    )
