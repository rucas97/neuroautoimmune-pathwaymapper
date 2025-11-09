FROM python:3.12-slim
WORKDIR /app

# --- System build dependencies ---
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        g++ \
        make \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        git \
        ca-certificates \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# ---------- Python dependencies ----------
COPY requirements.txt .
RUN pip install --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --force-reinstall -U six python-dateutil

# ---------- Project source ----------
COPY src/ ./src
COPY data/ ./data

CMD ["python", "src/pathway_mapper.py"]
