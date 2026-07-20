# ChiralFold — stereochemistry toolkit container image
# Usage:
#   docker run --rm -v "$PWD:/data" ghcr.io/tommaso-r-marena/chiralfold:3.4.0 audit /data/protein.pdb
FROM python:3.11-slim-bookworm

LABEL org.opencontainers.image.source="https://github.com/Tommaso-R-Marena/ChiralFold"
LABEL org.opencontainers.image.description="ChiralFold stereochemistry toolkit"
LABEL org.opencontainers.image.licenses="MIT"

RUN apt-get update && apt-get install -y --no-install-recommends \
    libxrender1 \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY pyproject.toml setup.py README.md LICENSE ./
COPY chiralfold/ chiralfold/
COPY web/ web/

RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir "rdkit>=2023.3,<2027" \
    && pip install --no-cache-dir .

WORKDIR /data

ENTRYPOINT ["chiralfold"]
CMD ["--help"]
