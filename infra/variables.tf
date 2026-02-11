variable "project_name" {
  type        = string
  description = "Project name used as prefix for all resources"
  default     = "bmii"
}

variable "environment" {
  type        = string
  description = "Deployment environment (development, staging, production)"
  default     = "production"
}

variable "aws_region" {
  type    = string
  default = "us-east-1"
}

variable "vpc_cidr" {
  type    = string
  default = "10.0.0.0/16"
}

variable "azs" {
  type    = list(string)
  default = ["us-east-1a", "us-east-1b"]
}

variable "db_instance_class" {
  type    = string
  default = "db.t3.micro"
}

variable "image_tag" {
  type    = string
  default = "latest"
}

variable "ecs_cpu" {
  type    = number
  default = 1024
}

variable "ecs_memory" {
  type    = number
  default = 2048
}

variable "ecs_desired_count" {
  type    = number
  default = 1
}

variable "api_key" {
  type        = string
  description = "API key for X-API-Key authentication"
  sensitive   = true
  default     = ""
}
