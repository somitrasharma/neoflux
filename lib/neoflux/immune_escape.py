"""
Immune escape feature extraction for NeoFlux.

Detects and quantifies tumor immune escape mechanisms that affect
neoantigen presentation:
    - HLA loss of heterozygosity (LOH)
    - Beta-2-microglobulin (B2M) mutations
    - Antigen processing machinery (APM) gene expression defects
    - Proteasomal cleavage and TAP transport prediction
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from neoflux.models import (
    AntigenProcessingFeatures,
    NeoantigenCandidate,
)
from neoflux.utils import APM_GENES

logger = logging.getLogger("neoflux")

# Amino acid weights for C-terminal proteasomal cleavage prediction.
# Based on the MAPPP/NetChop simplified model: positive weights favor
# cleavage after that residue.
PROTEASOME_CLEAVAGE_WEIGHTS = {
    "A": 0.2, "R": -0.3, "N": -0.1, "D": -0.2, "C": 0.1,
    "Q": -0.1, "E": -0.3, "G": 0.1, "H": -0.2, "I": 0.3,
    "L": 0.4, "K": -0.4, "M": 0.2, "F": 0.3, "P": -0.5,
    "S": 0.0, "T": 0.0, "W": 0.1, "Y": 0.2, "V": 0.3,
}

# TAP transport efficiency weights (simplified from IEDB TAP predictor).
# Higher values indicate better transport through TAP.
TAP_TRANSPORT_WEIGHTS = {
    "A": 0.1, "R": 0.3, "N": -0.2, "D": -0.3, "C": 0.0,
    "Q": -0.1, "E": -0.2, "G": -0.1, "H": 0.0, "I": 0.3,
    "L": 0.4, "K": 0.2, "M": 0.2, "F": 0.5, "P": -0.3,
    "S": -0.1, "T": -0.1, "W": 0.3, "Y": 0.4, "V": 0.2,
}


class ImmuneEscapeAnalyzer:
    """
    Analyzes tumor immune escape features relevant to neoantigen presentation.

    This module integrates data from multiple sources:
    1. LOHHLA output (HLA allele-specific copy number and LOH calls)
    2. VCF annotations (B2M mutations)
    3. RNA-seq expression (APM gene expression levels)
    4. Peptide sequence (proteasome cleavage and TAP transport predictions)
    """

    def __init__(self):
        self._hla_loh_calls: dict[str, bool] = {}
        self._apm_expression: dict[str, float] = {}
        self._b2m_mutated: bool = False

    def load_lohhla_results(self, path: str | Path) -> dict[str, bool]:
        """
        Parse LOHHLA output to determine HLA allele LOH status.

        LOHHLA outputs a TSV with columns including:
            HLA_A_type1, HLA_A_type2, pval_unique, ...
            copy_number_type1, copy_number_type2
            LOH (TRUE/FALSE)

        This method extracts per-allele LOH calls.

        Args:
            path: Path to LOHHLA output TSV.

        Returns:
            Dict mapping HLA allele string to LOH status (True/False).

        NOTE: This process is called by modules/hla/lohhla.nf. The Nextflow
        process invokes LOHHLA with tumor BAM, normal BAM, HLA types, and
        purity/ploidy estimates. See docs/tool_integration.md for details.
        """
        path = Path(path)
        if not path.exists():
            logger.warning("LOHHLA results not found: %s", path)
            return {}

        with open(path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                # LOHHLA reports pairs of alleles per gene
                for type_col in ["HLA_type1", "HLA_type2"]:
                    allele = row.get(type_col, "").strip()
                    if not allele:
                        continue

                    # Determine LOH status
                    pval = float(row.get("PVal_unique", row.get("pval", 1.0)))
                    loh_flag = row.get("LOH", "FALSE").strip().upper() == "TRUE"

                    # Consider LOH if p-value < 0.01 and LOH flag is set
                    is_loh = loh_flag and pval < 0.01
                    self._hla_loh_calls[allele] = is_loh

        logger.info(
            "LOHHLA results: %d alleles, %d with LOH",
            len(self._hla_loh_calls),
            sum(1 for v in self._hla_loh_calls.values() if v),
        )
        return self._hla_loh_calls

    def check_b2m_mutation(self, vcf_path: str | Path) -> bool:
        """
        Check whether B2M is somatically mutated in the tumor VCF.

        B2M (beta-2-microglobulin) is essential for MHC-I surface
        expression. Loss-of-function mutations in B2M are a major
        immune escape mechanism.

        Args:
            vcf_path: Path to annotated somatic VCF.

        Returns:
            True if B2M has a somatic mutation.

        NOTE: In the full pipeline, VEP annotations are used to identify
        B2M mutations. This is a simplified parser that searches for
        B2M in the INFO/ANN field. Integrate with pysam for production use.
        """
        path = Path(vcf_path)
        if not path.exists():
            logger.warning("VCF not found for B2M check: %s", path)
            return False

        with open(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                # Check if B2M appears in the VEP annotation field
                if "B2M" in line:
                    # Look for functional consequences
                    for term in [
                        "stop_gained", "frameshift", "splice_donor",
                        "splice_acceptor", "start_lost",
                    ]:
                        if term in line:
                            self._b2m_mutated = True
                            logger.info("B2M loss-of-function mutation detected")
                            return True

        return False

    def load_apm_expression(
        self,
        expression_path: str | Path,
        gene_column: str = "Name",
        tpm_column: str = "TPM",
    ) -> dict[str, float]:
        """
        Load expression levels of antigen processing machinery genes.

        Args:
            expression_path: Path to Salmon quant.sf or similar gene
                expression TSV.
            gene_column: Column name for gene symbol.
            tpm_column: Column name for TPM values.

        Returns:
            Dict mapping APM gene symbol to TPM.

        NOTE: Low expression of key APM genes (TAP1, TAP2, B2M, ERAP1)
        indicates impaired antigen processing and presentation.
        """
        path = Path(expression_path)
        if not path.exists():
            logger.warning("Expression file not found: %s", path)
            return {}

        with open(path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                gene = row.get(gene_column, "").strip()
                if gene.upper() in APM_GENES:
                    tpm = float(row.get(tpm_column, 0.0))
                    self._apm_expression[gene.upper()] = tpm

        logger.info(
            "Loaded expression for %d/%d APM genes",
            len(self._apm_expression), len(APM_GENES),
        )
        return self._apm_expression

    def predict_proteasome_cleavage(self, peptide: str) -> float:
        """
        Predict C-terminal proteasomal cleavage score for a peptide.

        Uses a simplified position-weighted scoring based on the
        C-terminal residue and the two preceding residues (P1, P2, P3).
        Higher scores indicate more likely cleavage.

        Args:
            peptide: Amino acid sequence.

        Returns:
            Cleavage score (arbitrary units, higher = more cleavable).

        NOTE: For production, consider integrating NetChop or netchop-3.1
        server predictions. This simplified model captures the main
        trends (preference for hydrophobic C-terminal residues).
        """
        if len(peptide) < 3:
            return 0.0

        # C-terminal residue (P1) is the cleavage site
        p1 = peptide[-1]
        p2 = peptide[-2]
        p3 = peptide[-3]

        score = (
            PROTEASOME_CLEAVAGE_WEIGHTS.get(p1, 0.0) * 1.0
            + PROTEASOME_CLEAVAGE_WEIGHTS.get(p2, 0.0) * 0.5
            + PROTEASOME_CLEAVAGE_WEIGHTS.get(p3, 0.0) * 0.3
        )
        return score

    def predict_tap_transport(self, peptide: str) -> float:
        """
        Predict TAP (Transporter associated with Antigen Processing)
        transport efficiency for a peptide.

        TAP preferentially transports peptides with hydrophobic or
        basic C-terminal residues and certain N-terminal residues.

        Args:
            peptide: Amino acid sequence.

        Returns:
            TAP transport score (arbitrary units, higher = better transport).

        NOTE: For production, consider integrating IEDB TAP prediction
        tools. This simplified model captures main trends.
        """
        if len(peptide) < 2:
            return 0.0

        # N-terminal and C-terminal residues are most important for TAP
        n_term = peptide[0]
        c_term = peptide[-1]
        p2 = peptide[1] if len(peptide) > 1 else ""
        penult = peptide[-2] if len(peptide) > 1 else ""

        score = (
            TAP_TRANSPORT_WEIGHTS.get(c_term, 0.0) * 1.0
            + TAP_TRANSPORT_WEIGHTS.get(n_term, 0.0) * 0.6
            + TAP_TRANSPORT_WEIGHTS.get(p2, 0.0) * 0.3
            + TAP_TRANSPORT_WEIGHTS.get(penult, 0.0) * 0.3
        )
        return score

    def annotate_candidates(
        self, candidates: list[NeoantigenCandidate]
    ) -> list[NeoantigenCandidate]:
        """
        Annotate candidates with immune escape and antigen processing features.

        For each candidate:
        1. Check if the best-binding HLA allele has LOH
        2. Compute proteasome cleavage and TAP transport scores
        3. Attach APM gene expression data
        4. Flag B2M mutation status

        Args:
            candidates: List of NeoantigenCandidate objects (modified in place).

        Returns:
            The same list with processing_features and presenting_allele_loh set.
        """
        loh_affected = 0
        for cand in candidates:
            # Check HLA LOH for the best-binding allele
            if cand.best_binding:
                allele = cand.best_binding.hla_allele
                if self._hla_loh_calls.get(allele, False):
                    cand.presenting_allele_loh = True
                    loh_affected += 1

            # Compute antigen processing features
            cand.processing_features = AntigenProcessingFeatures(
                proteasome_cleavage_score=self.predict_proteasome_cleavage(
                    cand.mutant_peptide
                ),
                tap_transport_score=self.predict_tap_transport(
                    cand.mutant_peptide
                ),
                b2m_mutated=self._b2m_mutated,
                apm_gene_expression=(
                    dict(self._apm_expression)
                    if self._apm_expression else None
                ),
            )

        logger.info(
            "Immune escape annotation: %d/%d candidates affected by HLA LOH, "
            "B2M mutated: %s",
            loh_affected, len(candidates), self._b2m_mutated,
        )
        return candidates
