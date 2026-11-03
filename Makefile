PYTHON ?= python3
PIP ?= pip3
NF ?= nextflow

.PHONY: install test train-demo stub-run demo-run docker-build clean lint format help download-refdata build-salmon-index build-star-index

help:
	@echo "NeoFlux development targets:"
	@echo "  make install      Install Python package in editable mode"
	@echo "  make test         Run pytest suite"
	@echo "  make test-cli     Run CLI smoke tests"
	@echo "  make train-demo   Train a demo ML model to models/demo_model.pkl"
	@echo "  make stub-run     Run Nextflow stub validation (requires nextflow)"
	@echo "  make demo-run     Run Nextflow demo profile (stub, 2 samples, no QC, no MS)"
	@echo "  make docker-build Build neoflux/core container image"
	@echo "  make lint         Run ruff linter"
	@echo "  make format       Run ruff formatter"
	@echo "  make clean        Remove Python bytecode, build artifacts, and test output"

install:
	cd lib && $(PIP) install -e ".[dev,ml]"

test:
	cd lib && $(PYTHON) -m pytest ../tests -v

test-cli:
	$(PYTHON) -m pytest tests/test_cli.py -v

train-demo:
	mkdir -p models
	cd lib && $(PYTHON) -m neoflux.cli train-model --output ../models/demo_model.pkl --n-samples 400

stub-run:
	$(NF) run main.nf -profile test,stub -stub-run

demo-run:
	$(NF) run main.nf -profile lite -stub-run

docker-build:
	docker build -t neoflux/core:latest -f containers/neoflux-core/Dockerfile .

lint:
	cd lib && ruff check neoflux/ ../tests/

format:
	cd lib && ruff format neoflux/ ../tests/

download-refdata:
	bash scripts/download_refdata.sh

build-salmon-index:
	mkdir -p /Volumes/OMNYCS/neoflux_data/indexes
	docker run --rm --entrypoint gffread \
		-v /Volumes/OMNYCS/neoflux_data:/data \
		quay.io/biocontainers/gffread:0.12.9--hf426362_0 \
		-w /data/indexes/transcripts.fa \
		-g /data/genome/GRCh38.fa /data/genome/GRCh38.gtf
	docker run --rm \
		-v /Volumes/OMNYCS/neoflux_data:/data \
		combinelab/salmon:1.10.3 \
		salmon index -t /data/indexes/transcripts.fa \
		-i /data/indexes/salmon_index -p 4

build-star-index:
	docker run --rm \
		-v /Volumes/OMNYCS/neoflux_data:/data \
		trinityctat/starfusion:1.12.0 \
		STAR --runMode genomeGenerate \
		--genomeDir /data/indexes/star_index \
		--genomeFastaFiles /data/genome/GRCh38.fa \
		--sjdbGTFfile /data/genome/GRCh38.gtf \
		--runThreadN 8

clean:
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	rm -rf build/ dist/ models/ test_results/ demo_results/ .ruff_cache/
