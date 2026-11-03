"""
Feature engineering for the NeoFlux ML ranking model.

Assembles a feature matrix from NeoantigenCandidate objects, combining
binding predictions, expression, clonality, antigen processing, immune
escape, immunopeptidomics evidence, and sequence-level features into
a single numeric vector per candidate.
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

from neoflux.models import (
    BindingLevel,
    ClonalStatus,
    NeoantigenCandidate,
    VariantType,
)
from neoflux.utils import (
    DRIVER_GENES,
    compute_agretopicity,
    compute_hydrophobicity,
    compute_self_similarity,
)

logger = logging.getLogger("neoflux")

# Ordered list of feature names for the ML model. This defines the column
# order of the feature matrix and must stay synchronized with the model
# training code.
FEATURE_COLUMNS = [
    # Binding features
    "best_ic50",
    "best_percent_rank",
    "is_strong_binder",
    "is_weak_binder",
    "presentation_score",
    "pmhc_half_life",
    "stability_rank",
    "binding_promiscuity",

    # Expression features
    "gene_tpm",
    "transcript_tpm",
    "vaf",
    "rna_vaf",
    "mutant_expression",

    # Clonality features
    "is_clonal",
    "cancer_cell_fraction",

    # Variant type features (one-hot)
    "is_snv",
    "is_indel",
    "is_frameshift",
    "is_fusion",
    "is_splice",
    "is_ncorf",

    # Sequence features
    "peptide_length",
    "hydrophobicity",
    "self_similarity",
    "agretopicity",

    # Antigen processing features
    "proteasome_score",
    "tap_score",

    # Immune escape features
    "presenting_allele_loh",
    "b2m_mutated",

    # Immunopeptidomics features
    "ms_detected",
    "pepquery_score",
    "spectral_count",

    # Driver gene annotation
    "is_driver_gene",
]


class FeatureEngineer:
    """
    Transforms NeoantigenCandidate objects into a numeric feature matrix
    suitable for ML model input.

    The feature set is designed to capture the key biological determinants
    of neoantigen immunogenicity:
    1. MHC binding and presentation (necessary but not sufficient)
    2. Expression and clonality (tumor-specific availability)
    3. Antigen processing (can it reach the cell surface?)
    4. Immune escape (is the presenting allele still present?)
    5. Immunopeptidomics (empirical evidence of presentation)
    6. Sequence properties (foreignness, hydrophobicity)
    """

    def __init__(self):
        self.feature_columns = FEATURE_COLUMNS

    def extract_features(
        self, candidate: NeoantigenCandidate
    ) -> dict[str, Optional[float]]:
        """
        Extract all features from a single NeoantigenCandidate.

        Missing values are represented as None (which XGBoost/LightGBM
        handle natively as missing).

        Args:
            candidate: A fully annotated NeoantigenCandidate.

        Returns:
            Dict mapping feature name to float value (or None).
        """
        features: dict[str, Optional[float]] = {}

        # -- Binding features --
        bp = candidate.best_binding
        if bp:
            features["best_ic50"] = bp.ic50
            features["best_percent_rank"] = bp.percent_rank
            features["is_strong_binder"] = float(
                bp.binding_level == BindingLevel.STRONG_BINDER
            ) if bp.binding_level else None
            features["is_weak_binder"] = float(
                bp.binding_level == BindingLevel.WEAK_BINDER
            ) if bp.binding_level else None
            features["presentation_score"] = bp.presentation_score
            features["pmhc_half_life"] = bp.half_life
            features["stability_rank"] = bp.stability_rank
        else:
            for col in [
                "best_ic50", "best_percent_rank", "is_strong_binder",
                "is_weak_binder", "presentation_score", "pmhc_half_life",
                "stability_rank",
            ]:
                features[col] = None

        # Binding promiscuity: number of alleles with rank < 2%
        n_binding = 0
        for pred in candidate.binding_predictions:
            if pred.percent_rank is not None and pred.percent_rank < 2.0:
                n_binding += 1
        features["binding_promiscuity"] = float(n_binding)

        # -- Expression features --
        expr = candidate.expression
        if expr:
            features["gene_tpm"] = expr.gene_tpm
            features["transcript_tpm"] = expr.transcript_tpm
            features["vaf"] = expr.vaf
            features["rna_vaf"] = expr.rna_vaf
            features["mutant_expression"] = expr.mutant_expression
        else:
            for col in [
                "gene_tpm", "transcript_tpm", "vaf", "rna_vaf",
                "mutant_expression",
            ]:
                features[col] = None

        # -- Clonality features --
        features["is_clonal"] = float(
            candidate.clonal_status == ClonalStatus.CLONAL
        )
        features["cancer_cell_fraction"] = candidate.cancer_cell_fraction

        # -- Variant type (one-hot) --
        vt = candidate.peptide_source
        features["is_snv"] = float(vt == VariantType.SNV)
        features["is_indel"] = float(
            vt in (VariantType.INSERTION, VariantType.DELETION)
        )
        features["is_frameshift"] = float(vt == VariantType.FRAMESHIFT)
        features["is_fusion"] = float(vt == VariantType.FUSION)
        features["is_splice"] = float(vt == VariantType.SPLICE)
        features["is_ncorf"] = float(vt == VariantType.NON_CANONICAL_ORF)

        # -- Sequence features --
        features["peptide_length"] = float(candidate.peptide_length)
        features["hydrophobicity"] = (
            candidate.hydrophobicity
            if candidate.hydrophobicity is not None
            else compute_hydrophobicity(candidate.mutant_peptide)
        )
        features["self_similarity"] = candidate.self_similarity
        features["agretopicity"] = candidate.agretopicity

        # -- Antigen processing --
        proc = candidate.processing_features
        if proc:
            features["proteasome_score"] = proc.proteasome_cleavage_score
            features["tap_score"] = proc.tap_transport_score
        else:
            features["proteasome_score"] = None
            features["tap_score"] = None

        # -- Immune escape --
        features["presenting_allele_loh"] = float(
            candidate.presenting_allele_loh
        )
        features["b2m_mutated"] = (
            float(proc.b2m_mutated) if proc else None
        )

        # -- Immunopeptidomics --
        ms = candidate.ms_evidence
        if ms:
            features["ms_detected"] = float(ms.detected)
            features["pepquery_score"] = ms.pepquery_score
            features["spectral_count"] = float(ms.spectral_count)
        else:
            features["ms_detected"] = None
            features["pepquery_score"] = None
            features["spectral_count"] = None

        # -- Driver gene --
        features["is_driver_gene"] = float(candidate.is_driver_gene)

        return features

    def build_feature_matrix(
        self, candidates: list[NeoantigenCandidate]
    ) -> pd.DataFrame:
        """
        Build a feature matrix (DataFrame) from a list of candidates.

        Each row corresponds to one candidate. Columns follow the order
        defined in FEATURE_COLUMNS. Missing values are represented as
        NaN (compatible with XGBoost/LightGBM missing value handling).

        Args:
            candidates: List of NeoantigenCandidate objects.

        Returns:
            DataFrame with shape (n_candidates, n_features).
        """
        rows = []
        for cand in candidates:
            feat = self.extract_features(cand)
            rows.append(feat)

        df = pd.DataFrame(rows, columns=self.feature_columns)
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        logger.info(
            "Built feature matrix: %d candidates x %d features",
            len(df), len(self.feature_columns),
        )
        return df

    def enrich_candidates(
        self, candidates: list[NeoantigenCandidate]
    ) -> list[NeoantigenCandidate]:
        """
        Compute derived features and attach them to candidate objects.

        This fills in hydrophobicity, self_similarity, agretopicity, and
        is_driver_gene if they are not already set.

        Args:
            candidates: List of NeoantigenCandidate objects (modified in place).

        Returns:
            The same list, with derived fields populated.
        """
        for cand in candidates:
            # Hydrophobicity
            if cand.hydrophobicity is None:
                cand.hydrophobicity = compute_hydrophobicity(
                    cand.mutant_peptide
                )

            # Self-similarity
            if cand.self_similarity is None:
                if cand.wildtype_peptide:
                    cand.self_similarity = compute_self_similarity(
                        cand.mutant_peptide, cand.wildtype_peptide
                    )
                else:
                    # No wildtype = maximally foreign (fusion, nCORF).
                    # These peptides are entirely novel sequences with no
                    # central tolerance — a key immunogenicity advantage.
                    cand.self_similarity = 0.0

            # Agretopicity
            if (
                cand.agretopicity is None
                and cand.best_binding
                and cand.best_binding.ic50 is not None
                and cand.wildtype_peptide
            ):
                wt_ic50 = None
                for pred in cand.binding_predictions:
                    if pred.peptide == cand.wildtype_peptide and pred.ic50 is not None:
                        wt_ic50 = pred.ic50
                        break
                if wt_ic50 is not None:
                    cand.agretopicity = compute_agretopicity(
                        cand.best_binding.ic50,
                        wt_ic50,
                    )
            # Fallback: estimate agretopicity from sequence dissimilarity
            if cand.agretopicity is None and cand.self_similarity is not None:
                cand.agretopicity = max(1.0, (1.0 / max(cand.self_similarity, 0.01)) - 0.5)

            # Driver gene
            if cand.variant.gene_symbol:
                cand.is_driver_gene = (
                    cand.variant.gene_symbol.upper() in DRIVER_GENES
                )

        return candidates
