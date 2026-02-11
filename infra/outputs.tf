output "alb_url" {
  description = "Public URL for the API (ALB DNS)"
  value       = module.ecs.alb_url
}

output "ecr_repository_url" {
  description = "ECR repository URL for docker push"
  value       = module.ecr.repository_url
}

output "rds_endpoint" {
  description = "RDS PostgreSQL endpoint"
  value       = module.database.endpoint
}

output "s3_buckets" {
  description = "S3 bucket names"
  value       = module.storage.bucket_names
}

output "ecs_cluster" {
  description = "ECS cluster name"
  value       = module.ecs.cluster_name
}

output "ecs_service" {
  description = "ECS service name"
  value       = module.ecs.service_name
}
