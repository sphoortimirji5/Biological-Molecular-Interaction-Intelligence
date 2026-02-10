FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libpq-dev \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install Poetry
RUN curl -sSL https://install.python-poetry.org | python3 -
ENV PATH="/root/.local/bin:$PATH"

# Configure Poetry to not use virtualenvs inside container
RUN poetry config virtualenvs.create false

# Copy dependency definition
COPY pyproject.toml poetry.lock* ./

# ── PyTorch: CPU-only by default (fast local builds) ──
# Override for GPU/CUDA in production:
#   docker build --build-arg TORCH_INDEX_URL=https://download.pytorch.org/whl/cu121 .
ARG TORCH_INDEX_URL=https://download.pytorch.org/whl/cpu
RUN pip install torch --index-url ${TORCH_INDEX_URL} --no-cache-dir

# Install remaining dependencies (torch already satisfied, skipped)
RUN poetry install --no-interaction --no-ansi --no-root

# Copy application code
COPY . .

# Run application
CMD ["uvicorn", "src.main:app", "--host", "0.0.0.0", "--port", "8000", "--reload"]
