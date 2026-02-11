# Migration Guide: Local → Production (AWS)

Step-by-step plan for migrating the Biological Molecular Interaction Intelligence platform from a local Docker Compose stack to AWS production infrastructure.

---

## Current State (Local)

| Component | Local | Production Target |
|-----------|-------|-------------------|
| API | `uvicorn --reload` (Docker) | ECS Fargate (multi-worker) |
| Database | PostgreSQL 15 (Docker) | RDS PostgreSQL 15 (private subnet) |
| Object Storage | MinIO (Docker) | S3 |
| Secrets | `.env` file | SSM Parameter Store |
| Load Balancer | Direct port 8000 | Application Load Balancer |
| Auth | Disabled (`api_key=""`) | `X-API-Key` header enforced |
| Logging | Console (colorized) | JSON → CloudWatch Logs (ECS default) |
| Metrics | `/metrics` endpoint | Prometheus → Grafana |

---

## Migration Checklist

### Phase 1: Infrastructure Provisioning

```bash
cd infra
cp terraform.tfvars.example terraform.tfvars
```

- [ ] **Configure `terraform.tfvars`**
  - Set `api_key` to a strong secret
  - Choose `aws_region` and `azs`
  - Select `db_instance_class` (start with `db.t3.micro`, upgrade as needed)
  - Set `ecs_cpu` / `ecs_memory` (1024/2048 minimum for ML models)

- [ ] **Provision AWS resources**
  ```bash
  terraform init
  terraform plan    # Review changes
  terraform apply   # Create infrastructure
  ```

- [ ] **Verify outputs**
  ```bash
  terraform output alb_url             # Public API endpoint
  terraform output ecr_repository_url  # Docker push target
  terraform output rds_endpoint        # Database host
  terraform output s3_buckets          # Bucket names
  ```

### Phase 2: Database Migration

- [ ] **Run Alembic migrations against RDS**
  ```bash
  RDS_ENDPOINT=$(cd infra && terraform output -raw rds_endpoint)

  DATABASE_URL="postgresql+psycopg2://bmii_admin:<password>@${RDS_ENDPOINT}/bio_bind_rank" \
    poetry run alembic upgrade head
  ```

- [ ] **Seed reference data** (if needed)
  - Export local Postgres data: `pg_dump -Fc bio_bind_rank > local.dump`
  - Import to RDS: `pg_restore --host=<rds-endpoint> --dbname=bio_bind_rank local.dump`

- [ ] **Verify schema**
  ```bash
  psql "postgresql://bmii_admin:<password>@${RDS_ENDPOINT}/bio_bind_rank" \
    -c "\dt"
  ```

### Phase 3: Object Storage Migration

MinIO → S3 migration:

- [ ] **Sync local MinIO buckets to S3**
  ```bash
  # For each bucket: raw, processed, artifacts, features
  aws s3 sync s3://raw s3://bmii-raw --source-region us-east-1
  ```

  Or if using local MinIO with `mc`:
  ```bash
  mc alias set local http://localhost:9000 minioadmin minioadmin
  mc mirror local/raw s3/bmii-raw
  mc mirror local/processed s3/bmii-processed
  mc mirror local/artifacts s3/bmii-artifacts
  mc mirror local/features s3/bmii-features
  ```

- [ ] **Update config references**
  - `OBJECT_STORE_ENDPOINT` → no longer needed (native S3)
  - `OBJECT_STORE_ACCESS_KEY/SECRET_KEY` → IAM role (no static credentials)
  - `OBJECT_STORE_USE_SSL` → `true`

### Phase 4: Container Build & Push

- [ ] **Authenticate with ECR**
  ```bash
  ECR_URL=$(cd infra && terraform output -raw ecr_repository_url)
  REGION=$(cd infra && terraform output -raw aws_region 2>/dev/null || echo "us-east-1")

  aws ecr get-login-password --region $REGION | \
    docker login --username AWS --password-stdin $ECR_URL
  ```

- [ ] **Build production image**
  ```bash
  docker build -t bmii:latest .
  docker tag bmii:latest $ECR_URL:latest
  docker push $ECR_URL:latest
  ```

- [ ] **Verify image in ECR**
  ```bash
  aws ecr describe-images --repository-name bmii
  ```

### Phase 5: Deploy & Verify

- [ ] **Deploy ECS service**
  ```bash
  aws ecs update-service \
    --cluster bmii \
    --service bmii \
    --force-new-deployment
  ```

- [ ] **Monitor deployment**
  ```bash
  aws ecs wait services-stable --cluster bmii --services bmii
  ```

- [ ] **Smoke tests**
  ```bash
  ALB_URL=$(cd infra && terraform output -raw alb_url)

  # Health
  curl -s $ALB_URL/health | jq .

  # Docs
  curl -s $ALB_URL/openapi.json | jq .info

  # Authenticated request
  curl -s -X POST $ALB_URL/rank \
    -H "Content-Type: application/json" \
    -H "X-API-Key: <your-key>" \
    -d '{"smiles": "CC(=O)Oc1ccccc1C(=O)O", "top_k": 5}' | jq .
  ```

---

## Configuration Changes Summary

| Variable | Local (`.env`) | Production (SSM / ECS env) |
|----------|----------------|---------------------------|
| `APP_ENV` | `development` | `production` |
| `DATABASE_URL` | `postgresql+psycopg2://postgres:postgres@postgres:5432/bio_bind_rank` | `postgresql+psycopg2://bmii_admin:<pwd>@<rds-endpoint>/bio_bind_rank` (from SSM) |
| `OBJECT_STORE_ENDPOINT` | `minio:9000` | *(not needed — native S3)* |
| `OBJECT_STORE_ACCESS_KEY` | `minioadmin` | *(not needed — IAM role)* |
| `OBJECT_STORE_SECRET_KEY` | `minioadmin` | *(not needed — IAM role)* |
| `OBJECT_STORE_USE_SSL` | `false` | `true` |
| `API_KEY` | *(empty)* | Strong secret (from SSM) |
| `LOG_LEVEL` | `DEBUG` | `INFO` |
| `OTEL_EXPORTER` | `console` | `otlp` |

---

## Code Changes Required

### S3 Client (Zero-Credential Mode)

The `S3StorageProvider` currently uses explicit `access_key` / `secret_key`. In production with ECS task IAM roles, the boto3 client should fall back to the credential chain:

```python
# src/ingestion/s3_provider.py — production path
if settings.object_store_access_key:
    # Local dev: explicit MinIO credentials
    client = boto3.client("s3",
        endpoint_url=f"http://{settings.object_store_endpoint}",
        aws_access_key_id=settings.object_store_access_key,
        aws_secret_access_key=settings.object_store_secret_key,
    )
else:
    # Production: IAM role credentials (automatic)
    client = boto3.client("s3", region_name=settings.object_store_region)
```

### Logging Mode

Already handled — `APP_ENV=production` switches structlog to JSON output automatically.

---

## Rollback Plan

If production deployment fails:

1. **ECS**: Roll back to previous task definition revision
   ```bash
   aws ecs update-service --cluster bmii --service bmii \
     --task-definition bmii:<previous-revision>
   ```

2. **Database**: RDS has 7-day automated backups
   ```bash
   aws rds restore-db-instance-to-point-in-time \
     --source-db-instance-identifier bmii-postgres \
     --target-db-instance-identifier bmii-postgres-restored \
     --restore-time <timestamp>
   ```

3. **Infrastructure**: `terraform plan` shows current vs desired state; revert `.tfvars` and `apply`

---

## Post-Migration Tasks

- [ ] Set up Grafana alerting rules (5xx rate, latency p99, ECS CPU) via Prometheus datasource
- [ ] Enable RDS Performance Insights
- [ ] Configure S3 bucket policies for cross-account access (if needed)
- [ ] Add HTTPS (ACM certificate + ALB HTTPS listener)
- [ ] Set up CI/CD pipeline (GitHub Actions → ECR → ECS)
- [ ] Enable Terraform remote state (S3 + DynamoDB locking)
- [ ] Point OTel collector to Grafana Cloud or self-hosted Grafana stack
