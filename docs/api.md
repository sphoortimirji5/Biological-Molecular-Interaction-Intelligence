# API Documentation

## Swagger UI

Interactive API documentation is auto-generated from the endpoint definitions and Pydantic schemas.

### Start the Server

```bash
# Kill any process on port 8000, then start
lsof -ti :8000 | xargs kill -9 2>/dev/null
uvicorn src.main:app --reload
```

### Open Swagger

Once the server is running, open:

- **Swagger UI**: [http://localhost:8000/docs](http://localhost:8000/docs)
- **ReDoc** (alternative viewer): [http://localhost:8000/redoc](http://localhost:8000/redoc)
- **OpenAPI JSON**: [http://localhost:8000/openapi.json](http://localhost:8000/openapi.json)

---

## Endpoints

| Method | Path | Tag | Description |
|--------|------|-----|-------------|
| `GET` | `/` | General | Welcome message — confirms API is reachable |
| `GET` | `/health` | Operations | Readiness probe — checks DB, S3, and model status |
| `GET` | `/metrics` | Operations | Prometheus-compatible metrics scrape endpoint |
| `POST` | `/rank` | Inference | Rank protein targets for a query compound |

---

## GET /health — Example

**200 (healthy):**
```json
{
  "status": "healthy",
  "checks": {
    "database": "ok",
    "object_store": "ok",
    "model_loaded": true
  }
}
```

**503 (degraded):**
```json
{
  "status": "degraded",
  "checks": {
    "database": "connection refused",
    "object_store": "ok",
    "model_loaded": false
  }
}
```

---

## POST /rank — Example

**Request:**
```bash
curl -X POST http://localhost:8000/rank \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CC(=O)Oc1ccccc1C(=O)O", "top_k": 3}'
```

**Response:**
```json
{
  "smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "top_k": 3,
  "model_version": "v1",
  "results": [
    { "rank": 1, "external_id": "P31749", "source": "uniprot", "score": 0.943 },
    { "rank": 2, "external_id": "P00533", "source": "uniprot", "score": 0.871 },
    { "rank": 3, "external_id": "P04637", "source": "uniprot", "score": 0.812 }
  ]
}
```

**Common SMILES for testing:**

| Drug | SMILES |
|------|--------|
| Aspirin | `CC(=O)Oc1ccccc1C(=O)O` |
| Caffeine | `CN1C=NC2=C1C(=O)N(C(=O)N2C)C` |
| Ibuprofen | `CC(C)Cc1ccc(cc1)C(C)C(=O)O` |

**Error Response (422):**
```json
{
  "detail": "Invalid SMILES: 'NOT_A_MOLECULE' could not be parsed by RDKit"
}
```
