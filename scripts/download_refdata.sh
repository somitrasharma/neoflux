#!/usr/bin/env bash
# NeoFlux reference data downloader
# Usage: bash scripts/download_refdata.sh /Volumes/OMNYCS/neoflux_data
#
# Downloads GRCh38 reference genome, GTF annotation, VEP cache,
# extracts transcriptome, and builds salmon index.
set -euo pipefail

DATA_DIR="${1:-/Volumes/OMNYCS/neoflux_data}"
ENSEMBL_RELEASE="${ENSEMBL_RELEASE:-111}"
mkdir -p "$DATA_DIR"/{genome,vep_cache,indexes}

if command -v curl &>/dev/null; then
    DL="curl -sL -o"
elif command -v wget &>/dev/null; then
    DL="wget -q --show-progress -O"
else
    echo "ERROR: neither curl nor wget is installed"
    exit 1
fi

echo "=== NeoFlux Reference Data Downloader ==="
echo "Target: $DATA_DIR"
echo "Ensembl release: $ENSEMBL_RELEASE"
echo ""

# -- GRCh38 primary assembly FASTA (~841 MB compressed, ~3.1 GB decompressed) --
FASTA="$DATA_DIR/genome/GRCh38.fa"

if [ -f "$FASTA" ] && [ -s "$FASTA" ]; then
    echo "[SKIP] $FASTA exists"
else
    echo "[DOWNLOAD] GRCh38 primary assembly FASTA (~841 MB) ..."
    $DL /tmp/GRCh38.fa.gz \
        "https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    echo "[DECOMPRESS] to /tmp ..."
    gunzip -c /tmp/GRCh38.fa.gz > /tmp/GRCh38.fa
    echo "[COPY] to $FASTA ..."
    cp /tmp/GRCh38.fa "$FASTA"
    rm -f /tmp/GRCh38.fa.gz /tmp/GRCh38.fa
    echo "[OK] $FASTA ($(du -h "$FASTA" | cut -f1))"
fi

# -- GRCh38 GTF annotation (~52 MB compressed, ~1.4 GB decompressed) --
GTF="$DATA_DIR/genome/GRCh38.gtf"

if [ -f "$GTF" ] && [ -s "$GTF" ]; then
    echo "[SKIP] $GTF exists"
else
    echo "[DOWNLOAD] GRCh38 GTF annotation (~52 MB) ..."
    $DL /tmp/GRCh38.gtf.gz \
        "https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gtf.gz"
    echo "[DECOMPRESS] to /tmp ..."
    gunzip -c /tmp/GRCh38.gtf.gz > /tmp/GRCh38.gtf
    echo "[COPY] to $GTF ..."
    cp /tmp/GRCh38.gtf "$GTF"
    rm -f /tmp/GRCh38.gtf.gz /tmp/GRCh38.gtf
    echo "[OK] $GTF ($(du -h "$GTF" | cut -f1))"
fi

# -- VEP cache (~27 GB compressed, ~50 GB decompressed) --
VEP_DIR="$DATA_DIR/vep_cache"
VEP_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/variation/indexed_vep_cache/homo_sapiens_vep_${ENSEMBL_RELEASE}_GRCh38.tar.gz"

if [ -d "$VEP_DIR" ] && ls -A "$VEP_DIR" 2>/dev/null | grep -q .; then
    echo "[SKIP] VEP cache exists in $VEP_DIR"
else
    echo "[DOWNLOAD] VEP cache for GRCh38 (~27 GB) ..."
    if command -v aria2c &>/dev/null; then
        aria2c --dir="$VEP_DIR" --out=vep_cache.tar.gz \
            --continue=true --max-connection-per-server=4 --split=4 \
            --min-split-size=50M --file-allocation=none \
            --connect-timeout=30 --timeout=600 --max-tries=0 \
            --retry-wait=60 --auto-file-renaming=false \
            --allow-overwrite=true --summary-interval=30 \
            "$VEP_URL"
    else
        $DL "$VEP_DIR/vep_cache.tar.gz" "$VEP_URL"
    fi
    echo "[EXTRACT] to $VEP_DIR ..."
    tar -xzf "$VEP_DIR/vep_cache.tar.gz" -C "$VEP_DIR" --strip-components=1
    rm -f "$VEP_DIR/vep_cache.tar.gz"
    echo "[OK] VEP cache extracted to $VEP_DIR"
fi

# -- Transcriptome extraction (for salmon index) --
TRANSCRIPTS="$DATA_DIR/indexes/transcripts.fa"

if [ -f "$TRANSCRIPTS" ] && [ -s "$TRANSCRIPTS" ]; then
    echo "[SKIP] $TRANSCRIPTS exists"
else
    echo "[EXTRACT] transcript sequences from GTF + FASTA (gffread) ..."
    if command -v docker &>/dev/null; then
        docker run --rm --entrypoint gffread \
            -v "$DATA_DIR:/data" \
            quay.io/biocontainers/gffread:0.12.9--hf426362_0 \
            -w /data/indexes/transcripts.fa \
            -g /data/genome/GRCh38.fa \
            /data/genome/GRCh38.gtf
    else
        echo "[WARN] docker not available, skipping transcriptome extraction"
        echo "  Install gffread and run:"
        echo "  gffread -w $TRANSCRIPTS -g $FASTA $GTF"
    fi
    echo "[OK] $TRANSCRIPTS ($(du -h "$TRANSCRIPTS" | cut -f1))"
fi

# -- Salmon index --
SALMON_INDEX="$DATA_DIR/indexes/salmon_index"

if [ -d "$SALMON_INDEX" ] && ls -A "$SALMON_INDEX" 2>/dev/null | grep -q .; then
    echo "[SKIP] salmon index exists at $SALMON_INDEX"
else
    echo "[BUILD] salmon index ..."
    if [ -f "$TRANSCRIPTS" ] && [ -s "$TRANSCRIPTS" ]; then
        if command -v docker &>/dev/null; then
            mkdir -p "$SALMON_INDEX"
            docker run --rm -v "$DATA_DIR:/data" \
                combinelab/salmon:1.10.3 \
                salmon index -t /data/indexes/transcripts.fa \
                -i /data/indexes/salmon_index -p 4
        else
            echo "[WARN] docker not available, skipping salmon index"
        fi
    else
        echo "[WARN] transcriptome not found, skipping salmon index"
    fi
    echo "[OK] salmon index at $SALMON_INDEX"
fi

echo ""
echo "=== Done ==="
echo "GRCh38 FASTA:      $FASTA"
echo "GRCh38 GTF:        $GTF"
echo "VEP cache:         $VEP_DIR"
echo "Transcriptome:     $TRANSCRIPTS"
echo "Salmon index:      $SALMON_INDEX"
echo ""
echo "To build STAR index (requires >30 GB RAM):"
echo "  docker run --rm -v $DATA_DIR:/data trinityctat/starfusion:1.12.0 STAR \\"
echo "    --runMode genomeGenerate \\"
echo "    --genomeDir /data/indexes/star_index \\"
echo "    --genomeFastaFiles /data/genome/GRCh38.fa \\"
echo "    --sjdbGTFfile /data/genome/GRCh38.gtf \\"
echo "    --runThreadN 8"
