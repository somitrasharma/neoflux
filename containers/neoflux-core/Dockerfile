## NeoFlux runtime image (SITC 2026)
## -----------------------------------
## Reference data (GRCh38 FASTA, GTF, VEP cache) is volume-mounted at
## /opt/neoflux/data to keep the image small and avoid 30GB+ downloads
## during docker build. See scripts/download_refdata.sh to prepare data.
##
## Build:   docker build -t neoflux/core:latest .
## Run:     docker run --rm -v /path/to/refdata:/opt/neoflux/data neoflux/core neoflux ...
## -----------------------------------------------------------------------

FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
WORKDIR /opt/neoflux

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates && rm -rf /var/lib/apt/lists/*

COPY lib/pyproject.toml /opt/neoflux/
COPY lib/neoflux /opt/neoflux/neoflux

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir ".[ml]"

RUN mkdir -p /opt/neoflux/data/genome /opt/neoflux/data/vep_cache

ENTRYPOINT ["neoflux"]
