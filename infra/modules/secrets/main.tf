##############################################################################
# Secrets Module — SSM Parameter Store
##############################################################################

variable "project_name" {
  type = string
}

variable "database_url" {
  type      = string
  sensitive = true
}

variable "api_key" {
  type      = string
  sensitive = true
  default   = ""
}

# ── Parameters ──

resource "aws_ssm_parameter" "database_url" {
  name  = "/${var.project_name}/DATABASE_URL"
  type  = "SecureString"
  value = var.database_url

  tags = { Name = "${var.project_name}-db-url" }
}

resource "aws_ssm_parameter" "api_key" {
  name  = "/${var.project_name}/API_KEY"
  type  = "SecureString"
  value = var.api_key != "" ? var.api_key : "change-me-in-production"

  tags = { Name = "${var.project_name}-api-key" }
}

# ── Outputs ──

output "database_url_arn" {
  value = aws_ssm_parameter.database_url.arn
}

output "api_key_arn" {
  value = aws_ssm_parameter.api_key.arn
}

output "secret_arns" {
  value = [
    aws_ssm_parameter.database_url.arn,
    aws_ssm_parameter.api_key.arn,
  ]
}
