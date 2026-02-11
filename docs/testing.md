# Testing Guide

> **57 tests** across 12 test modules — covering featurization, scoring, training, inference, and ingestion.

---

## Quick Reference

```bash
# All unit tests (no DB/S3 required)
python3 -m pytest tests/ -v -m "not integration"

# Only inference tests
python3 -m pytest tests/test_inference_service.py -v

# Integration tests (requires Postgres + MinIO via docker-compose)
python3 -m pytest tests/ -v -m integration
```

---

## Test Architecture

| Tier | Marker | Dependencies | Count |
|------|--------|-------------|-------|
| **Unit** | _(default)_ | None — all external deps mocked | 52 |
| **Integration** | `@pytest.mark.integration` | Postgres, MinIO, migrations applied | 5 |

---

## Unit Test Cases

### Inference Service — `test_inference_service.py`

| # | Test | Use Case | Validates |
|---|------|----------|-----------|
| 1 | `test_rank_returns_top_k` | Happy path: rank proteins for a valid compound | Returns exactly `top_k` results with rank, external_id, score |
| 2 | `test_rank_invalid_smiles` | Malformed SMILES input | Raises `ValueError` with descriptive message |
| 3 | `test_rank_respects_top_k` | Requesting fewer results than candidates | Returns exactly the requested number |
| 4 | `test_rank_top_k_larger_than_candidates` | `top_k` exceeds available proteins | Returns all candidates without error |
| 5 | `test_rank_scores_are_descending` | Score ordering | Results are sorted highest-score-first |
| 6 | `test_rank_empty_feature_store` | No protein embeddings loaded | Raises `ValueError` — not a silent empty list |
| 7 | `test_rank_preserves_protein_identity` | Returned IDs match feature store | External IDs in results ⊆ stored proteins |
| 8 | `test_rank_passes_correct_vectors_to_scorer` | Internal wiring | Query = 2048-dim fingerprint, candidates = N×emb_dim matrix |

---

### Morgan Fingerprints — `test_morgan.py`

| # | Test | Use Case | Validates |
|---|------|----------|-----------|
| 1 | `test_transform_produces_correct_shape` | Generate fingerprints from SMILES | Output has `fp_0`…`fp_2047` columns |
| 2 | `test_fingerprints_are_deterministic` | Same molecule → same fingerprint | Two runs produce identical values |
| 3 | `test_fingerprints_are_binary` | Fingerprint values are 0 or 1 | All values ∈ {0, 1} |
| 4 | `test_different_molecules_produce_different_fingerprints` | Distinct compounds yield distinct vectors | Aspirin ≠ Caffeine fingerprints |
| 5 | `test_invalid_smiles_are_skipped` | Graceful handling of bad SMILES | Invalid rows dropped, valid rows retained |
| 6 | `test_feature_manifest` | Manifest metadata | Contains radius, n_bits, generator name |

---

### Protein Embeddings — `test_protein_embeddings.py`

| # | Test | Use Case | Validates |
|---|------|----------|-----------|
| 1 | `test_transform_produces_correct_shape` | ESM-2 embedding generation | Output has `emb_0`…`emb_N` columns |
| 2 | `test_embeddings_are_deterministic` | Same sequence → same embedding | Two runs produce identical vectors |
| 3 | `test_embeddings_are_float` | Embedding dtype | All values are floats, not ints |
| 4 | `test_different_proteins_produce_different_embeddings` | Distinct sequences yield distinct vectors | Two different proteins ≠ same embedding |
| 5 | `test_feature_manifest` | Manifest metadata | Contains model name and embedding dim |

---

### Scoring Strategies — `test_scoring.py`

| # | Test | Use Case | Validates |
|---|------|----------|-----------|
| 1 | `test_similarity_scorer` | Cosine similarity scoring | Identical vectors → score ≈ 1.0 |
| 2 | `test_learned_scorer_produces_probabilities` | XGBoost-based scoring | All scores ∈ [0, 1] |
| 3 | `test_learned_scorer_concatenates_features` | Feature vector assembly | Query + candidate concatenated before prediction |
| 4 | `test_scoring_strategy_is_abstract` | ABC enforcement | Cannot instantiate `ScoringStrategy` directly |

---

### XGBoost Model — `test_xgboost_model.py`

| # | Test | Use Case | Validates |
|---|------|----------|-----------|
| 1 | `test_train_and_predict` | Train → predict cycle | Predictions have `probability` column ∈ [0, 1] |
| 2 | `test_predict_without_training_raises` | Predict before train | Raises error — no silent fallback |
| 3 | `test_save_and_load` | Model serialization | Loaded model produces identical predictions |
| 4 | `test_model_signature` | Artifact signing | Signature dict has `sha256`, `feature_version` |
| 5 | `test_custom_params` | Parameter override | Custom hyperparameters override defaults |

---

### Model Evaluator — `test_evaluator.py`

| # | Test | Use Case | Validates |
|---|------|----------|-----------|
| 1 | `test_perfect_predictions` | Oracle predictions | AUC-ROC = 1.0, accuracy = 1.0 |
| 2 | `test_random_predictions` | Noise predictions | Metrics are defined (not NaN) but not perfect |
| 3 | `test_ranking_metrics` | Ranking quality | Report includes `average_precision` |
| 4 | `test_report_structure` | Report schema | Report dict has expected keys |
| 5 | `test_single_class_returns_none` | All-positive / all-negative data | Gracefully returns `None` for undefined metrics |
| 6 | `test_compare_reports` | Model comparison | Comparison dict shows deltas between two reports |

---

### Training Data Builder — `test_training_data.py`

| # | Test | Use Case | Validates |
|---|------|----------|-----------|
| 1 | `test_build_produces_correct_splits` | Train/test split | Both splits are non-empty DataFrames |
| 2 | `test_feature_columns_are_concatenated` | Feature assembly | Compound fp + protein emb columns present |
| 3 | `test_no_data_leakage` | Split integrity | Train and test sets share no compound IDs |
| 4 | `test_stratification_preserves_label_ratio` | Balanced splits | Label ratio similar in train and test |
| 5 | `test_metadata_is_populated` | Build metadata | Metadata has feature versions and counts |
| 6 | `test_empty_interactions_returns_empty` | No interaction data | Returns empty DataFrames, not errors |

---

### Training Manager — `test_trainer.py`

| # | Test | Use Case | Validates |
|---|------|----------|-----------|
| 1 | `test_run_produces_result` | End-to-end training | Returns result with model + report |
| 2 | `test_model_and_report_persisted` | Artifact persistence | Model and report saved to expected paths |
| 3 | `test_empty_data_raises` | No training data | Raises error — not silent no-op |
| 4 | `test_evaluation_has_metrics` | Post-training evaluation | Evaluation report has AUC-ROC, accuracy |

---

### Feature Store — `test_feature_store.py`

| # | Test | Use Case | Validates |
|---|------|----------|-----------|
| 1 | `test_save_features_and_manifest` | Save Parquet + manifest to S3 | Correct keys, compressed content |
| 2 | `test_load_features` | Load Parquet from S3 | DataFrame matches original |
| 3 | `test_load_manifest` | Load manifest JSON from S3 | Manifest dict matches original |

---

### Featurization Manager — `test_featurization_manager.py`

| # | Test | Use Case | Validates |
|---|------|----------|-----------|
| 1 | `test_featurize_compounds` | Compound featurization pipeline | Generates and saves Morgan fingerprints |
| 2 | `test_featurize_proteins` | Protein featurization pipeline | Generates and saves ESM-2 embeddings |
| 3 | `test_full_run` | End-to-end featurization | Both compound + protein features generated |
| 4 | `test_empty_compounds_returns_zero` | No compounds in DB | Returns 0 count, no crash |
| 5 | `test_empty_proteins_returns_zero` | No proteins in DB | Returns 0 count, no crash |

---

## Integration Test Cases

Integration tests connect to **real** Postgres and MinIO services. They will **fail** if these services are not running.

### Prerequisites

| Service | Purpose | Default Port |
|---------|---------|-------------|
| **PostgreSQL** | Stores compounds, proteins, ingestion runs | `5432` |
| **MinIO** (S3-compatible) | Stores raw files, feature Parquet files | `9000` |

### How to Run

```bash
# Step 1: Free up ports used by the project (kill any conflicting processes)
#         Ports: 5432 (Postgres), 8000 (App), 9000-9001 (MinIO)
lsof -ti :5432,:8000,:9000,:9001 | xargs kill -9 2>/dev/null

# macOS Homebrew Postgres often auto-starts and shadows Docker on port 5432:
brew services stop postgresql@17 2>/dev/null  # or postgresql@14, postgresql@16, etc.

# Step 2: Start infrastructure
docker compose up -d

# Step 3: Wait for healthy containers
docker compose ps   # Verify postgres and minio show "healthy"

# Step 4: Create database tables
python3 -m alembic upgrade head

# Step 5: Run integration tests
python3 -m pytest tests/ -v -m integration

# Step 6: (Optional) Tear down
docker compose down
```

> [!IMPORTANT]
> If you skip Steps 1–4 and run `python3 -m pytest tests/ -v` (without the `-m "not integration"` filter),
> you will see failures like:
> ```
> asyncpg.exceptions.InvalidCatalogNameError: database "bio_bind_rank" does not exist
> ```
> This means Postgres is either not running or the database has not been created.

> [!WARNING]
> **macOS users**: Homebrew may auto-start a local PostgreSQL (`postgresql@14`, `@17`, etc.) on port 5432.
> This shadows the Docker Postgres container and causes the `database does not exist` error even when
> Docker is running. Always run `brew services stop postgresql@17` (or your version) before starting Docker.
>
> To check: `lsof -i :5432` — if you see a `postgres` process that is NOT Docker, stop it first.

### Troubleshooting

| Error | Cause | Fix |
|-------|-------|-----|
| `database "bio_bind_rank" does not exist` | Local Postgres shadowing Docker, or DB not created | `brew services stop postgresql@17` then `docker compose up -d` |
| `port is already allocated` | Another process using the same port | `lsof -ti :5432,:8000,:9000,:9001 \| xargs kill -9` |
| `connection refused` on port 5432 | Postgres container hasn't started | `docker compose ps` — wait for healthy status |
| `relation "compounds" does not exist` | Migrations not applied | `python3 -m alembic upgrade head` |
| `connection refused` on port 9000 | MinIO container hasn't started | `docker compose ps` — wait for healthy status |

### Test Cases

#### Ingestion Pipeline — `test_ingestion_integration.py`

| # | Test | Use Case | Validates |
|---|------|----------|-----------|
| 1 | `test_compound_ingestion_flow` | SDF file → compounds table | Parsed compound in DB, run marked COMPLETED |
| 2 | `test_protein_ingestion_flow` | FASTA file → proteins table | Parsed protein in DB with sequence + ID |
| 3 | `test_idempotency_no_duplicates` | Ingest same source twice | Exactly 1 compound row (upsert), 2 run rows |
| 4 | `test_failure_path_fetch_throws` | Network failure during fetch | Run marked FAILED with error details |

#### Featurization Pipeline — `test_featurization_integration.py`

| # | Test | Use Case | Validates |
|---|------|----------|-----------|
| 1 | `test_compound_featurization_end_to_end` | DB compounds → Parquet in S3 | Feature file written to MinIO |

---

## Coverage by Component

| Component | Unit Tests | Integration Tests | Total |
|-----------|-----------|------------------|-------|
| Inference Service | 8 | — | 8 |
| Morgan Fingerprints | 6 | — | 6 |
| Protein Embeddings | 5 | — | 5 |
| Scoring Strategies | 4 | — | 4 |
| XGBoost Model | 5 | — | 5 |
| Model Evaluator | 6 | — | 6 |
| Training Data | 6 | — | 6 |
| Training Manager | 4 | — | 4 |
| Feature Store | 3 | — | 3 |
| Featurization Manager | 5 | 1 | 6 |
| Ingestion Pipeline | — | 4 | 4 |
| **Total** | **52** | **5** | **57** |
