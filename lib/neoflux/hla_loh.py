"""Coverage-based HLA LOH detection.

Uses tumor/normal BAM comparison at the MHC locus (chr6:28-34 Mb)
with autosomal normalization to detect allelic loss of heterozygosity.
This implementation follows the same coverage-ratio principle used by
HLA-LOH callers such as LOHHLA and published allele-specific LOH
studies.

Source: McGranahan et al. 2017, Science; Shukla et al. 2015, Nat
Biotechnology; Schumacher and Schreiber 2015, Science.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path

logger = logging.getLogger("neoflux")

# MHC Class I genomic coordinates (GRCh38)
# HLA-A: chr6:29,942,767-29,945,885
# HLA-B: chr6:31,353,875-31,357,201
# HLA-C: chr6:31,268,749-31,271,676
_HLA_LOCUS = {
    "HLA-A": ("chr6", 29942767, 29945885),
    "HLA-B": ("chr6", 31353875, 31357201),
    "HLA-C": ("chr6", 31268749, 31271676),
    "HLA-E": ("chr6", 30489557, 30494750),
    "HLA-F": ("chr6", 29723845, 29728529),
    "HLA-G": ("chr6", 29794602, 29800833),
}

# MHC region (for genome-wide baseline normalization)
_MHC_REGION = ("chr6", 28000000, 34000000)
_AUTOSOMAL_REGIONS = [
    ("chr1", 1000000, 5000000),
    ("chr2", 1000000, 5000000),
    ("chr12", 1000000, 5000000),
    ("chr22", 1000000, 5000000),
]


def _region_coverage(
    bam_path: str, chrom: str, start: int, end: int, min_mapq: int = 20
) -> float:
    """Compute mean coverage over a genomic region."""
    import pysam

    bam = pysam.AlignmentFile(bam_path, "rb")
    total_depth = 0
    n_bases = 0
    for pileup_col in bam.pileup(chrom, start, end, min_mapping_quality=min_mapq, truncate=True):
        total_depth += pileup_col.nsegments
        n_bases += 1
    bam.close()
    return total_depth / max(n_bases, 1)


def compute_hla_loh(
    tumor_bam: str,
    normal_bam: str,
    hla_alleles: list[str],
    min_coverage: int = 20,
) -> dict[str, dict[str, object]]:
    """Compute HLA LOH status from tumor/normal BAM coverage.

    For each HLA gene present in the allele list, compares the tumor
    to normal coverage ratio against a genome-wide baseline. A
    significant reduction in the tumor HLA ratio indicates LOH.

    Args:
        tumor_bam: Path to tumor BAM.
        normal_bam: Path to normal BAM.
        hla_alleles: List of HLA alleles (e.g., HLA-A*02:01).
        min_coverage: Minimum normal coverage for reliable LOH calls.

    Returns:
        Dict mapping allele string to dict with keys:
            is_loh (bool), loh_pvalue (float), copy_number (float).
    """
    if not Path(tumor_bam).exists() or not Path(normal_bam).exists():
        logger.warning("Tumor or normal BAM missing; returning all-False LOH")
        return {a: {"is_loh": False, "loh_pvalue": 1.0, "copy_number": 2.0}
                for a in hla_alleles}

    try:
        import pysam
        # Quick validation
        pysam.AlignmentFile(tumor_bam, "rb").close()
        pysam.AlignmentFile(normal_bam, "rb").close()
    except Exception:
        logger.warning("Cannot open BAM files for LOH detection; returning all-False")
        return {a: {"is_loh": False, "loh_pvalue": 1.0, "copy_number": 2.0}
                for a in hla_alleles}

    # Compute genome-wide baseline coverage ratio
    normal_baseline = sum(
        _region_coverage(normal_bam, c, s, e) for c, s, e in _AUTOSOMAL_REGIONS
    ) / len(_AUTOSOMAL_REGIONS)
    tumor_baseline = sum(
        _region_coverage(tumor_bam, c, s, e) for c, s, e in _AUTOSOMAL_REGIONS
    ) / len(_AUTOSOMAL_REGIONS)

    baseline_ratio = tumor_baseline / max(normal_baseline, 1)

    logger.info(
        "LOH baseline: normal_depth=%.1f tumor_depth=%.1f ratio=%.2f",
        normal_baseline, tumor_baseline, baseline_ratio,
    )

    results: dict[str, dict[str, object]] = {}
    for allele in hla_alleles:
        gene = allele.split("*")[0] if "*" in allele else allele.split(":")[0]
        gene = gene.replace("HLA-", "")

        locus = _HLA_LOCUS.get(f"HLA-{gene}")
        if locus is None:
            results[allele] = {"is_loh": False, "loh_pvalue": 1.0, "copy_number": 2.0}
            continue

        chrom, start, end = locus
        normal_cov = _region_coverage(normal_bam, chrom, start, end)
        tumor_cov = _region_coverage(tumor_bam, chrom, start, end)

        # Expected tumor coverage = normal_coverage * baseline_ratio (if no LOH)
        expected_tumor = normal_cov * baseline_ratio
        observed_ratio = tumor_cov / max(normal_cov, 1)
        normalized_ratio = observed_ratio / max(baseline_ratio, 0.01)

        # Approximate uncertainty from coverage and use a one-sided LOH z-score.
        # This mirrors published allele-specific LOH callers that compare HLA
        # locus depletion against genome-wide or local coverage baselines.
        variance = max(1.0, expected_tumor)
        stderr_ratio = math.sqrt(variance) / max(normal_cov, 1)
        zscore = (normalized_ratio - 1.0) / max(stderr_ratio, 0.05)

        is_loh = normalized_ratio < 0.70 and zscore < -1.65 and normal_cov >= min_coverage
        if is_loh:
            copy_number = max(0, min(2, round(normalized_ratio * 2)))
        else:
            copy_number = 2

        # One-sided lower-tail p-value for reduced tumor coverage.
        pvalue = 0.5 * (1.0 + math.erf(zscore / math.sqrt(2)))
        pvalue = max(1e-6, min(1.0, pvalue))

        logger.info(
            "LOH %s: normal_cov=%.1f tumor_cov=%.1f norm_ratio=%.2f loh=%s",
            allele, normal_cov, tumor_cov, normalized_ratio, is_loh,
        )

        results[allele] = {
            "is_loh": is_loh,
            "loh_pvalue": pvalue,
            "copy_number": copy_number,
        }

    return results
