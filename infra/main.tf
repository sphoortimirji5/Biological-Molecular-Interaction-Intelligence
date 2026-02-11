##############################################################################
# Biological Molecular Interaction Intelligence — AWS Infrastructure
#
# Usage:
#   cp terraform.tfvars.example terraform.tfvars
#   terraform init
#   terraform plan
#   terraform apply
##############################################################################

terraform {
  required_version = ">= 1.5"

  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
    random = {
      source  = "hashicorp/random"
      version = "~> 3.6"
    }
  }

  # Remote state (uncomment for team use)
  # backend "s3" {
  #   bucket         = "bmii-terraform-state"
  #   key            = "infra/terraform.tfstate"
  #   region         = "us-east-1"
  #   dynamodb_table = "bmii-terraform-lock"
  #   encrypt        = true
  # }
}

provider "aws" {
  region = var.aws_region

  default_tags {
    tags = {
      Project     = var.project_name
      Environment = var.environment
      ManagedBy   = "terraform"
    }
  }
}

# ── Modules ──

module "networking" {
  source       = "./modules/networking"
  project_name = var.project_name
  vpc_cidr     = var.vpc_cidr
  azs          = var.azs
}

module "database" {
  source             = "./modules/database"
  project_name       = var.project_name
  private_subnet_ids = module.networking.private_subnet_ids
  rds_sg_id          = module.networking.rds_sg_id
  db_instance_class  = var.db_instance_class
}

module "storage" {
  source       = "./modules/storage"
  project_name = var.project_name
}

module "ecr" {
  source       = "./modules/ecr"
  project_name = var.project_name
}

module "secrets" {
  source       = "./modules/secrets"
  project_name = var.project_name
  database_url = module.database.connection_url
  api_key      = var.api_key
}

module "ecs" {
  source              = "./modules/ecs"
  project_name        = var.project_name
  vpc_id              = module.networking.vpc_id
  public_subnet_ids   = module.networking.public_subnet_ids
  private_subnet_ids  = module.networking.private_subnet_ids
  alb_sg_id           = module.networking.alb_sg_id
  ecs_sg_id           = module.networking.ecs_sg_id
  ecr_repository_url  = module.ecr.repository_url
  image_tag           = var.image_tag
  cpu                 = var.ecs_cpu
  memory              = var.ecs_memory
  desired_count       = var.ecs_desired_count

  environment_variables = [
    { name = "APP_ENV", value = var.environment },
    { name = "APP_PORT", value = "8000" },
    { name = "LOG_LEVEL", value = "INFO" },
    { name = "OBJECT_STORE_TYPE", value = "s3" },
    { name = "OBJECT_STORE_REGION", value = var.aws_region },
    { name = "OBJECT_STORE_USE_SSL", value = "true" },
    { name = "OBJECT_STORE_BUCKET_RAW", value = "${var.project_name}-raw" },
    { name = "OBJECT_STORE_BUCKET_PROCESSED", value = "${var.project_name}-processed" },
    { name = "OBJECT_STORE_BUCKET_ARTIFACTS", value = "${var.project_name}-artifacts" },
    { name = "FEATURE_STORE_BUCKET", value = "${var.project_name}-features" },
    { name = "OTEL_EXPORTER", value = "otlp" },
    { name = "OTEL_SERVICE_NAME", value = var.project_name },
  ]

  secrets = [
    { name = "DATABASE_URL", valueFrom = module.secrets.database_url_arn },
    { name = "API_KEY", valueFrom = module.secrets.api_key_arn },
  ]
}
