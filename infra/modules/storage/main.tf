##############################################################################
# Storage Module — S3 Buckets
##############################################################################

variable "project_name" {
  type = string
}

variable "bucket_names" {
  type    = list(string)
  default = ["raw", "processed", "artifacts", "features"]
}

resource "aws_s3_bucket" "buckets" {
  count  = length(var.bucket_names)
  bucket = "${var.project_name}-${var.bucket_names[count.index]}"

  tags = { Name = "${var.project_name}-${var.bucket_names[count.index]}" }
}

resource "aws_s3_bucket_versioning" "buckets" {
  count  = length(var.bucket_names)
  bucket = aws_s3_bucket.buckets[count.index].id

  versioning_configuration {
    status = "Enabled"
  }
}

resource "aws_s3_bucket_server_side_encryption_configuration" "buckets" {
  count  = length(var.bucket_names)
  bucket = aws_s3_bucket.buckets[count.index].id

  rule {
    apply_server_side_encryption_by_default {
      sse_algorithm = "AES256"
    }
  }
}

resource "aws_s3_bucket_public_access_block" "buckets" {
  count  = length(var.bucket_names)
  bucket = aws_s3_bucket.buckets[count.index].id

  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

resource "aws_s3_bucket_lifecycle_configuration" "buckets" {
  count  = length(var.bucket_names)
  bucket = aws_s3_bucket.buckets[count.index].id

  rule {
    id     = "transition-to-ia"
    status = "Enabled"

    transition {
      days          = 90
      storage_class = "STANDARD_IA"
    }
  }
}

# ── Outputs ──

output "bucket_names" {
  value = { for i, name in var.bucket_names : name => aws_s3_bucket.buckets[i].bucket }
}

output "bucket_arns" {
  value = { for i, name in var.bucket_names : name => aws_s3_bucket.buckets[i].arn }
}
