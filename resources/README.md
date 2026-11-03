# Reference Data for NeoFlux (SITC 2026)

This folder contains **exact, version‑locked reference assets** required for a reproducible run of the NeoFlux pipeline.  All files are either directly downloaded from stable URLs or generated from those downloads (e.g. Salmon/STAR indexes).  The accompanying `download_reference.sh` script automates the whole process and verifies integrity via MD5 checksums.

## Directory layout

```
resources/
│   README.md               # you are here
│   download_reference.sh    # download & build script
│   checksums.md5           # MD5 values for every asset (generated after download)
│
├── genome/                 # reference FASTA & GTF
│   ├── GRCh38.fa
│   └── GRCh38.gtf
│
├── indexes/                # pre‑built aligner/index assets
│   ├── star_index/         # STAR genome index
│   └── salmon_index/       # Salmon quasi‑map index
│
└── vep/                    # VEP cache (Ensembl release 111, GRCh38)
    └── 111_GRCh38/         # full VEP cache directory
```

**Important:** the assets total several tens of gigabytes.  They are *not* stored in the Git history; instead they are generated on‑demand via the script.

## How to obtain the assets

```bash
# From the repository root (where the `resources/` folder lives):
cd resources
bash download_reference.sh
```

The script will:
1. Download the reference FASTA and GTF from Ensembl release 111.
2. Build the STAR and Salmon indexes (multithreaded, using `$THREADS` if set).
3. Download the pre‑built VEP cache for GRCh38/Ensembl 111.
4. Compute MD5 checksums for *all* files and store them in `checksums.md5`.

If any download fails, the script aborts and reports which step failed.

## Using the assets in NeoFlux

NeoFlux looks for the assets under `$NEOFUX_ROOT/resources/`.  The Dockerfile we added earlier copies the entire `resources/` directory into `/opt/neoflux/data/` and sets environment variables:

```Dockerfile
ENV REFERENCE_FA=/opt/neoflux/data/genome/GRCh38.fa
ENV REFERENCE_GTF=/opt/neoflux/data/genome/GRCh38.gtf
ENV VEP_CACHE=/opt/neoflux/data/vep/111_GRCh38
ENV STAR_INDEX=/opt/neoflux/data/indexes/star_index
ENV SALMON_INDEX=/opt/neoflux/data/indexes/salmon_index
```

Make sure the paths in your `nextflow.config` or command‑line arguments match those environment variables.

---

*If you already have these assets stored elsewhere (e.g., on a shared filesystem), you can symlink the `resources/` folder to that location and skip the download step.*
