# Production Deployment Plan

## Local → Production Matrix

| Concern | Local (docker-compose) | Production (AWS) |
|---|---|---|
| **Compute** | Single container | ECS Fargate (API) + ECS GPU tasks (featurization) |
| **Database** | Postgres container | RDS PostgreSQL (Multi-AZ) |
| **Object Storage** | MinIO container | S3 |
| **PyTorch** | CPU-only (200 MB) | GPU/CUDA 12.1 (2 GB) |
| **Secrets** | `.env` file | AWS Secrets Manager |
| **Observability** | Console exporter | OTLP → Grafana Cloud / CloudWatch |
| **CI/CD** | `git push` | GitHub Actions → ECR → ECS |

---

## Infrastructure

```mermaid
graph TD
    subgraph "VPC (Private Subnets)"
        ALB["Application Load Balancer"] --> API["ECS Fargate<br/>FastAPI Service"]
        API --> RDS["RDS PostgreSQL<br/>(Multi-AZ)"]
        API --> S3["S3 Buckets<br/>raw / features / models"]
        
        GPU["ECS GPU Task<br/>(p3.2xlarge)<br/>ESM-2 Embeddings"] --> S3
        GPU --> RDS
    end
    
    GHA["GitHub Actions"] -->|Push Image| ECR["ECR"]
    ECR --> API
    ECR --> GPU
    
    style ALB fill:#f9f,stroke:#333,color:#000
    style RDS fill:#e1d5e7,stroke:#333,color:#000
    style S3 fill:#fff2cc,stroke:#333,color:#000
    style GPU fill:#d5e8d4,stroke:#333,color:#000
```

---

## CI/CD Pipeline (GitHub Actions)

```yaml
# .github/workflows/deploy.yml (planned)
name: Deploy
on:
  push:
    branches: [main]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: docker-compose build
      - run: docker-compose run --rm app poetry run alembic upgrade head
      - run: docker-compose run --rm app poetry run pytest

  build-and-push:
    needs: test
    runs-on: ubuntu-latest
    steps:
      - uses: aws-actions/configure-aws-credentials@v4
      - uses: aws-actions/amazon-ecr-login@v2
      # API image (CPU)
      - run: docker build -t $ECR_REPO:$SHA .
      # Featurization image (GPU) — only when src/features/ changes
      - run: |
          docker build \
            --build-arg TORCH_INDEX_URL=https://download.pytorch.org/whl/cu121 \
            -t $ECR_REPO:$SHA-gpu .

  deploy:
    needs: build-and-push
    runs-on: ubuntu-latest
    steps:
      - run: aws ecs update-service --cluster prod --service api --force-new-deployment
```

---

## Environment Promotion

### Step 1: Secrets Migration
```bash
# Move from .env to Secrets Manager
aws secretsmanager create-secret \
  --name bio-intel/prod \
  --secret-string file://prod-secrets.json
```

Secrets needed:
- `DATABASE_URL` → RDS connection string
- `OBJECT_STORE_*` → S3 credentials (or use IAM role, preferred)
- `OTEL_EXPORTER_OTLP_ENDPOINT` → Grafana Cloud / ADOT endpoint

### Step 2: Database
```bash
# Run migrations against RDS
docker run --rm \
  -e DATABASE_URL=$RDS_URL \
  $ECR_IMAGE poetry run alembic upgrade head
```

### Step 3: Object Storage
No migration needed — S3 uses the same key layout as MinIO:
```
raw/<source>/<version>/<sha256>/<file>
features/<entity_type>/<version>/features.parquet
```

### Step 4: PyTorch Switch
```bash
# GPU image for featurization tasks
docker build --build-arg TORCH_INDEX_URL=https://download.pytorch.org/whl/cu121 .
```
API service continues using CPU image (no PyTorch needed at inference time for Phase 2).

---

## Config Changes Required (Zero Code Changes)

| Variable | Local Value | Production Value |
|---|---|---|
| `APP_ENV` | `development` | `production` |
| `DATABASE_URL` | `postgresql+psycopg2://postgres:postgres@postgres:5432/bio_bind_rank` | `postgresql+psycopg2://<user>:<pass>@<rds-host>:5432/bio_bind_rank` |
| `OBJECT_STORE_ENDPOINT` | `minio:9000` | _(omit — uses AWS S3 default)_ |
| `OBJECT_STORE_ACCESS_KEY` | `minioadmin` | IAM role (no key needed) |
| `OTEL_EXPORTER` | `console` | `otlp` |
| `OTEL_EXPORTER_OTLP_ENDPOINT` | `http://localhost:4317` | `https://otlp.grafana.net:4317` |
| `LOG_LEVEL` | `DEBUG` | `INFO` |

> [!IMPORTANT]
> **Zero code changes between local and production.** Everything is driven by environment variables defined in `src/config.py`. The same Docker image runs everywhere — only the config changes.

---

## Deployment Checklist

- [ ] Provision RDS PostgreSQL (Multi-AZ, `db.r6g.large`)
- [ ] Create S3 buckets (`raw`, `features`, `models`)
- [ ] Create ECR repository
- [ ] Set up GitHub Actions workflow
- [ ] Configure Secrets Manager
- [ ] Run Alembic migrations against RDS
- [ ] Deploy API service (ECS Fargate, CPU)
- [ ] Deploy featurization task (ECS GPU, on-demand)
- [ ] Configure ALB + Route53
- [ ] Verify OTEL traces in Grafana/Jaeger
- [ ] Smoke test: ingest → featurize → verify S3 artifacts
