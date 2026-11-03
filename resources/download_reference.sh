#!/usr/bin/env bash
# download_reference.sh – fetch and prepare all reference assets for NeoFlux (SITC 2026)
#
# This script is intentionally idempotent: re‑running it will overwrite any existing
# files with the freshly downloaded versions and recompute checksums.
#
# Prerequisites (all available in the Docker builder image or on a typical macOS/Linux host):
#   * wget, tar, gunzip, gzip, md5sum (or shasum -a 256 as fallback)
#   * STAR, salmon, and samtools (installed via conda in the builder stage)
#
# The script creates the following layout under the repository's "resources" folder:
#
# resources/
# ├── genome/          # GRCh38 FASTA + GTF
# │   ├── GRCh38.fa
# │   └── GRCh38.gtf
# ├── indexes/         # STAR and Salmon indexes
# │   ├── star_index/
# │   └── salmon_index/
# └── vep/             # VEP cache (Ensembl release 111, GRCh38)
#     └── 111_GRCh38/   # full cache directory
#
# Optional tarballs can be placed in the same directory before running the script:
#   * hla_hd_data.tar.gz   – contains the HLA‑HD allele database
#   * netmhcpan_license.tar.gz – contains NetMHCpan license files
# The script will extract them into the appropriate locations.
#
# After successful execution a "checksums.md5" file is written that lists the MD5
# checksum of every asset.  This file can be committed (small) to enable reviewers to
# verify exact versions.

set -euo pipefail

# ----------- Configuration ---------------------------------------------------
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENOME_DIR="${BASE_DIR}/genome"
INDEX_DIR="${BASE_DIR}/indexes"
VEP_DIR="${BASE_DIR}/vep"
CHECKSUM_FILE="${BASE_DIR}/checksums.md5"

# Number of threads for STAR/Salmon indexing – respect $THREADS env var if set
THREADS="${THREADS:-$(nproc)}"

mkdir -p "${GENOME_DIR}" "${INDEX_DIR}" "${VEP_DIR}"

# ----------- 1. Reference FASTA & GTF --------------------------------------
FA_URL="https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz"

echo "[download] GRCh38 FASTA…"
wget -q -O "${GENOME_DIR}/GRCh38.fa.gz" "${FA_URL}"
gunzip -c "${GENOME_DIR}/GRCh38.fa.gz" > "${GENOME_DIR}/GRCh38.fa" || {
    echo "[error] failed to unzip FASTA"; exit 1;
}

echo "[download] GRCh38 GTF…"
wget -q -O "${GENOME_DIR}/GRCh38.gtf.gz" "${GTF_URL}"
gunzip -c "${GENOME_DIR}/GRCh38.gtf.gz" > "${GENOME_DIR}/GRCh38.gtf" || {
    echo "[error] failed to unzip GTF"; exit 1;
}

# ----------- 2. STAR genome index ----------------------------------------
STAR_IDX="${INDEX_DIR}/star_index"
mkdir -p "${STAR_IDX}"
echo "[index] Building STAR index (threads=${THREADS})…"
STAR \
    --runThreadN ${THREADS} \
    --runMode genomeGenerate \
    --genomeDir "${STAR_IDX}" \
    --genomeFastaFiles "${GENOME_DIR}/GRCh38.fa" \
    --sjdbGTFfile "${GENOME_DIR}/GRCh38.gtf" \
    --sjdbOverhang 149

# ----------- 3. Salmon transcriptome index -------------------------------
SALMON_IDX="${INDEX_DIR}/salmon_index"
mkdir -p "${SALMON_IDX}"
echo "[index] Building Salmon index (threads=${THREADS})…"
salmon index \
    -t "${GENOME_DIR}/GRCh38.fa" \
    -g "${GENOME_DIR}/GRCh38.gtf" \
    -p ${THREADS} \
    -i "${SALMON_IDX}" \
    --type quasi

# ----------- 4. VEP cache (Ensembl 111, GRCh38) --------------------------
VEP_URL="https://ftp.ebi.ac.uk/pub/databases/vep/vep_cache/vep/homo_sapiens/111_GRCh38.tar.gz"
mkdir -p "${VEP_DIR}"
echo "[download] VEP cache…"
wget -q -O - "${VEP_URL}" | tar -xz -C "${VEP_DIR}" --strip-components=1

# ----------- 5. Optional: HLA‑HD data & NetMHCpan licence ----------------
if [[ -f "${BASE_DIR}/hla_hd_data.tar.gz" ]]; then
    echo "[extract] HLA‑HD data…"
    tar -xz -C "${BASE_DIR}" -f "${BASE_DIR}/hla_hd_data.tar.gz"
fi

if [[ -f "${BASE_DIR}/netmhcpan_license.tar.gz" ]]; then
    echo "[extract] NetMHCpan licence…"
    tar -xz -C "${BASE_DIR}" -f "${BASE_DIR}/netmhcpan_license.tar.gz"
fi

# ----------- 6. MD5 checksums --------------------------------------------
echo "[checksum] Computing MD5 sums…"
> "${CHECKSUM_FILE}"
find "${GENOME_DIR}" -type f -exec md5sum {} + >> "${CHECKSUM_FILE}"
find "${INDEX_DIR}" -type f -exec md5sum {} + >> "${CHECKSUM_FILE}"
find "${VEP_DIR}" -type f -exec md5sum {} + >> "${CHECKSUM_FILE}"
if [[ -f "${BASE_DIR}/hla_hd_data.tar.gz" ]]; then
    md5sum "${BASE_DIR}/hla_hd_data.tar.gz" >> "${CHECKSUM_FILE}"
fi
if [[ -f "${BASE_DIR}/netmhcpan_license.tar.gz" ]]; then
    md5sum "${BASE_DIR}/netmhcpan_license.tar.gz" >> "${CHECKSUM_FILE}"
fi

echo "[done] All reference assets are ready under ${BASE_DIR}"
