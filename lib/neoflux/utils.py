"""
Shared utilities for the NeoFlux pipeline.

Contains helper functions for file I/O, sequence operations, HLA parsing,
and other common operations used across modules.
"""

from __future__ import annotations

import csv
import hashlib
import logging
import re
from pathlib import Path
from typing import Optional

logger = logging.getLogger("neoflux")

# Standard amino acid alphabet
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# Kyte-Doolittle hydrophobicity scale
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# BLOSUM62 scoring matrix (subset -- full matrix loaded lazily)
# Used for self-similarity computation between WT and mutant peptides
_BLOSUM62: Optional[dict[tuple[str, str], int]] = None

# Curated cancer driver gene list (COSMIC Cancer Gene Census, tier 1)
DRIVER_GENES = frozenset({
    "TP53", "KRAS", "BRAF", "PIK3CA", "PTEN", "APC", "EGFR", "NRAS",
    "RB1", "NF1", "CDKN2A", "SMAD4", "CTNNB1", "IDH1", "IDH2",
    "BRCA1", "BRCA2", "ATM", "MYC", "ERBB2", "ALK", "ROS1", "MET",
    "FGFR1", "FGFR2", "FGFR3", "KIT", "PDGFRA", "RET", "ABL1",
    "JAK2", "NOTCH1", "FBXW7", "ARID1A", "KMT2D", "KMT2C",
    "SETD2", "BAP1", "VHL", "SMARCA4", "STK11", "KEAP1",
    "NFE2L2", "MAP2K1", "MAP3K1", "CDH1", "GATA3", "SF3B1",
    "U2AF1", "SRSF2", "DNMT3A", "TET2", "NPM1", "FLT3",
})

# Key antigen processing machinery genes
APM_GENES = frozenset({
    "TAP1", "TAP2", "TAPBP", "PSMB5", "PSMB8", "PSMB9",
    "B2M", "CALR", "CANX", "ERAP1", "ERAP2",
    "HLA-A", "HLA-B", "HLA-C",
})


def parse_hla_allele(allele_str: str) -> tuple[str, str]:
    """
    Parse an HLA allele string into (gene, allele) components.

    Handles formats like:
        HLA-A*02:01, A*02:01, HLA-A02:01

    Returns:
        Tuple of (gene, full_standardized_allele).
    """
    allele_str = allele_str.strip()

    # Standardize prefix
    if not allele_str.startswith("HLA-"):
        if allele_str[0] in "ABCDG" and ("*" in allele_str or ":" in allele_str):
            allele_str = f"HLA-{allele_str}"

    # Extract gene and specificity
    match = re.match(
        r"(HLA-[A-Z]+[0-9]?)\*?(\d{2}):?(\d{2}):?(\d{0,2}):?(\d{0,2})",
        allele_str
    )
    if not match:
        raise ValueError(f"Cannot parse HLA allele: {allele_str}")

    gene = match.group(1)
    field1 = match.group(2)
    field2 = match.group(3)
    standardized = f"{gene}*{field1}:{field2}"

    return gene, standardized


def compute_hydrophobicity(peptide: str) -> float:
    """
    Compute mean Kyte-Doolittle hydrophobicity for a peptide.

    Args:
        peptide: Amino acid sequence (uppercase single-letter).

    Returns:
        Mean hydrophobicity value.
    """
    values = [KYTE_DOOLITTLE.get(aa, 0.0) for aa in peptide.upper()]
    if not values:
        return 0.0
    return sum(values) / len(values)


def compute_self_similarity(mutant: str, wildtype: str) -> float:
    """
    Compute sequence similarity between mutant and wildtype peptides
    using BLOSUM62 alignment score, normalized by wildtype self-score.

    A score of 1.0 means identical; lower values indicate more
    dissimilar (potentially more immunogenic) mutations.

    Args:
        mutant: Mutant peptide sequence.
        wildtype: Wildtype peptide sequence.

    Returns:
        Normalized similarity score [0, 1+]. Values > 1 are possible
        for conservative substitutions that score higher than self.
    """
    if len(mutant) != len(wildtype):
        # For length-changing variants, similarity is low by definition
        return 0.0

    blosum = _get_blosum62()

    mut_score = sum(
        blosum.get((m, w), 0) for m, w in zip(mutant, wildtype)
    )
    wt_self_score = sum(
        blosum.get((w, w), 0) for w in wildtype
    )

    if wt_self_score == 0:
        return 0.0

    return mut_score / wt_self_score


def compute_agretopicity(mt_affinity: float, wt_affinity: float) -> float:
    """
    Compute agretopicity: ratio of wildtype to mutant binding affinity.

    Higher values indicate stronger differential binding (mutant binds
    better than wildtype), suggesting higher immunogenicity.

    Args:
        mt_affinity: Mutant peptide IC50 (nM).
        wt_affinity: Wildtype peptide IC50 (nM).

    Returns:
        Agretopicity ratio (wt_IC50 / mt_IC50). Values > 1 indicate
        the mutant binds more strongly.
    """
    if mt_affinity <= 0:
        return 0.0
    return wt_affinity / mt_affinity


def generate_candidate_id(
    sample_id: str,
    variant_id: str,
    peptide: str,
    hla_allele: str,
) -> str:
    """
    Generate a deterministic unique ID for a neoantigen candidate.

    Uses a hash of the key identifiers to create a short, reproducible ID.
    """
    raw = f"{sample_id}|{variant_id}|{peptide}|{hla_allele}"
    digest = hashlib.sha256(raw.encode()).hexdigest()[:12]
    return f"NF-{digest}"


def read_samplesheet(path: str | Path) -> list[dict[str, str]]:
    """
    Read a NeoFlux samplesheet CSV.

    Expected columns: sample_id, tumor_bam, normal_bam, rna_fastq_r1,
    rna_fastq_r2, hla_file, ms_mzml, ms_search_results, tumor_purity,
    tumor_ploidy

    Returns:
        List of dicts, one per sample.
    """
    samples = []
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Samplesheet not found: {path}")

    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Strip whitespace from keys and values
            cleaned = {k.strip(): v.strip() for k, v in row.items()}
            samples.append(cleaned)

    logger.info("Read %d samples from samplesheet: %s", len(samples), path)
    return samples


def validate_peptide(peptide: str) -> bool:
    """Check if a peptide contains only valid amino acid characters."""
    return all(aa in AMINO_ACIDS for aa in peptide.upper())


def _get_blosum62() -> dict[tuple[str, str], int]:
    """Lazily load a minimal BLOSUM62 matrix."""
    global _BLOSUM62
    if _BLOSUM62 is not None:
        return _BLOSUM62

    # Minimal BLOSUM62 diagonal + common substitutions
    # In production, load from a full matrix file
    _BLOSUM62 = {}
    diag = {
        "A": 4, "R": 5, "N": 6, "D": 6, "C": 9,
        "Q": 5, "E": 5, "G": 6, "H": 8, "I": 4,
        "L": 4, "K": 5, "M": 5, "F": 6, "P": 7,
        "S": 4, "T": 5, "W": 11, "Y": 7, "V": 4,
    }
    for aa, score in diag.items():
        _BLOSUM62[(aa, aa)] = score

    # Common substitutions (symmetric)
    subs = {
        ("A", "V"): 0, ("A", "G"): 0, ("D", "E"): 2, ("D", "N"): 1,
        ("F", "Y"): 3, ("I", "V"): 3, ("I", "L"): 2, ("K", "R"): 2,
        ("L", "M"): 2, ("L", "V"): 1, ("S", "T"): 1, ("N", "Q"): 0,
    }
    for (a1, a2), score in subs.items():
        _BLOSUM62[(a1, a2)] = score
        _BLOSUM62[(a2, a1)] = score

    return _BLOSUM62
