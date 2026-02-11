from fastapi import FastAPI, Response, HTTPException
from fastapi.openapi.utils import get_openapi
from prometheus_client import generate_latest, CONTENT_TYPE_LATEST
from src.logging_config import get_logger
from src.inference.schemas import (
    RankRequest, RankResponse, RankedTarget,
    SimilarRequest, SimilarResponse,
)

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
    response_description="Service health with component readiness",
    responses={
        200: {
            "description": "All components healthy",
            "content": {
                "application/json": {
                    "example": {
                        "status": "healthy",
                        "checks": {
                            "database": "ok",
                            "object_store": "ok",
                            "model_loaded": True,
                        },
                    }
                }
            },
        },
        503: {
            "description": "One or more components degraded",
            "content": {
                "application/json": {
                    "example": {
                        "status": "degraded",
                        "checks": {
                            "database": "connection refused",
                            "object_store": "ok",
                            "model_loaded": False,
                        },
                    }
                }
            },
        },
    },
)
async def health():
    """
    Readiness probe with component-level checks.

    - **database**: Runs `SELECT 1` to verify Postgres connectivity.
    - **object_store**: Calls `head_bucket` on the feature store bucket.
    - **model_loaded**: Reports whether the inference service has been initialized.

    Returns HTTP 200 if all infrastructure checks pass, HTTP 503 otherwise.
    """
    checks = {}
    all_ok = True

    # --- Database ---
    try:
        from src.db.database import AsyncSessionLocal
        from sqlalchemy import text

        async with AsyncSessionLocal() as session:
            await session.execute(text("SELECT 1"))
        checks["database"] = "ok"
    except Exception as e:
        checks["database"] = str(e)[:120]
        all_ok = False

    # --- Object Store (S3 / MinIO) ---
    try:
        from src.config import settings
        from src.ingestion.storage import S3StorageProvider

        storage = S3StorageProvider()
        storage.client.head_bucket(Bucket=settings.feature_store_bucket)
        checks["object_store"] = "ok"
    except Exception as e:
        checks["object_store"] = str(e)[:120]
        all_ok = False

    # --- Model ---
    checks["model_loaded"] = _inference_service is not None

    status_code = 200 if all_ok else 503
    return Response(
        content=__import__("json").dumps({
            "status": "healthy" if all_ok else "degraded",
            "checks": checks,
        }),
        media_type="application/json",
        status_code=status_code,
    )


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


# Lazy-initialized similarity service
_similarity_service = None


def _get_similarity_service():
    """Factory: build SimilarityService once from config, cache globally."""
    global _similarity_service
    if _similarity_service is not None:
        return _similarity_service

    from src.config import settings
    from src.features.store import FeatureStore
    from src.ingestion.storage import S3StorageProvider
    from src.inference.similarity import SimilarityService

    storage = S3StorageProvider()
    feature_store = FeatureStore(storage=storage, bucket=settings.feature_store_bucket)

    _similarity_service = SimilarityService(feature_store=feature_store)
    return _similarity_service


@app.post(
    "/similar",
    response_model=SimilarResponse,
    tags=["Inference"],
    summary="Find Similar Entities",
    response_description="Ranked list of similar proteins or compounds",
    responses={
        200: {
            "description": "Successful similarity search",
            "content": {
                "application/json": {
                    "example": {
                        "entity_type": "protein",
                        "query": "MKTAYIAKQRQISFVKSH",
                        "top_k": 3,
                        "similarity_metric": "cosine",
                        "results": [
                            {"rank": 1, "external_id": "P31749", "source": "uniprot", "score": 0.982},
                            {"rank": 2, "external_id": "P00533", "source": "uniprot", "score": 0.941},
                            {"rank": 3, "external_id": "P04637", "source": "uniprot", "score": 0.897},
                        ],
                    }
                }
            },
        },
        422: {
            "description": "Invalid query or unknown entity type",
            "content": {
                "application/json": {
                    "example": {"detail": "entity_type must be 'protein' or 'compound', got 'gene'"}
                }
            },
        },
    },
)
async def similar(request: SimilarRequest):
    """
    Find similar proteins or compounds.

    **Protein similarity** (`entity_type: "protein"`):
    - Query is an amino-acid sequence
    - Generates ESM-2 embedding on the fly
    - Ranks all stored proteins by **cosine similarity**

    **Compound similarity** (`entity_type: "compound"`):
    - Query is a SMILES string
    - Generates Morgan fingerprint (ECFP4) on the fly
    - Ranks all stored compounds by **Tanimoto similarity**
    """
    entity_type = request.entity_type.lower()

    if entity_type not in ("protein", "compound"):
        raise HTTPException(
            status_code=422,
            detail=f"entity_type must be 'protein' or 'compound', got '{request.entity_type}'",
        )

    try:
        service = _get_similarity_service()

        if entity_type == "protein":
            results = service.find_similar_proteins(sequence=request.query, top_k=request.top_k)
            metric = "cosine"
        else:
            results = service.find_similar_compounds(smiles=request.query, top_k=request.top_k)
            metric = "cosine"
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))

    return SimilarResponse(
        entity_type=entity_type,
        query=request.query,
        top_k=request.top_k,
        similarity_metric=metric,
        results=[RankedTarget(**r) for r in results],
    )
