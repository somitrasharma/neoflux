"""
Data models for NeoFlux neoantigen representation.

All core data structures use Pydantic for validation, serialization, and
schema enforcement. These models flow through the entire pipeline: from
variant calling output through feature engineering to final ranked output.
"""

from __future__ import annotations

import enum
from typing import Optional

from pydantic import BaseModel, Field


class VariantType(str, enum.Enum):
    """Enumeration of supported somatic variant classes."""
    SNV = "SNV"
    INSERTION = "insertion"
    DELETION = "deletion"
    FRAMESHIFT = "frameshift"
    FUSION = "fusion"
    SPLICE = "splice"
    NON_CANONICAL_ORF = "non_canonical_orf"


class VariantConsequence(str, enum.Enum):
    """Variant effect prediction consequence types (VEP-style)."""
    MISSENSE = "missense_variant"
    FRAMESHIFT = "frameshift_variant"
    INFRAME_INSERTION = "inframe_insertion"
    INFRAME_DELETION = "inframe_deletion"
    STOP_GAINED = "stop_gained"
    SPLICE_DONOR = "splice_donor_variant"
    SPLICE_ACCEPTOR = "splice_acceptor_variant"
    FUSION = "fusion"
    NON_CANONICAL = "non_canonical"
    OTHER = "other"


class HLAClass(str, enum.Enum):
    """HLA molecule class."""
    CLASS_I = "I"
    CLASS_II = "II"


class BindingLevel(str, enum.Enum):
    """MHC binding strength classification."""
    STRONG_BINDER = "SB"
    WEAK_BINDER = "WB"
    NON_BINDER = "NB"


class ClonalStatus(str, enum.Enum):
    """Clonal vs subclonal classification of a variant."""
    CLONAL = "clonal"
    SUBCLONAL = "subclonal"
    INDETERMINATE = "indeterminate"


class GenomicVariant(BaseModel):
    """A somatic variant identified in the tumor genome or transcriptome."""

    sample_id: str = Field(
        ..., description="Patient/sample identifier"
    )
    chromosome: str = Field(
        ..., description="Chromosome (e.g. chr1, chrX)"
    )
    position: int = Field(
        ..., description="1-based genomic position"
    )
    ref_allele: str = Field(
        ..., description="Reference allele"
    )
    alt_allele: str = Field(
        ..., description="Alternate (tumor) allele"
    )
    variant_type: VariantType = Field(
        ..., description="Classification of variant"
    )
    gene_symbol: Optional[str] = Field(
        None, description="HUGO gene symbol"
    )
    transcript_id: Optional[str] = Field(
        None, description="Ensembl transcript ID"
    )
    consequence: Optional[VariantConsequence] = Field(
        None, description="VEP-predicted consequence"
    )
    protein_change: Optional[str] = Field(
        None, description="HGVS protein notation (e.g. p.V600E)"
    )
    caller: Optional[str] = Field(
        None, description="Variant caller that identified this (mutect2, strelka2, etc.)"
    )
    is_consensus: bool = Field(
        False, description="True if called by multiple callers"
    )

    # Fusion-specific fields
    fusion_partner_5p: Optional[str] = Field(
        None, description="5-prime fusion partner gene"
    )
    fusion_partner_3p: Optional[str] = Field(
        None, description="3-prime fusion partner gene"
    )
    fusion_breakpoint: Optional[str] = Field(
        None, description="Breakpoint coordinates"
    )

    # Splice-specific fields
    splice_event_type: Optional[str] = Field(
        None, description="rMATS event type (SE, A5SS, A3SS, MXE, RI)"
    )
    delta_psi: Optional[float] = Field(
        None, description="Differential splicing (delta PSI)"
    )

    @property
    def variant_id(self) -> str:
        """Unique identifier for this variant."""
        return f"{self.chromosome}:{self.position}:{self.ref_allele}>{self.alt_allele}"


class HLAAllele(BaseModel):
    """An HLA allele with optional LOH status."""

    allele: str = Field(
        ..., description="4-digit HLA allele (e.g. HLA-A*02:01)"
    )
    hla_class: HLAClass = Field(
        ..., description="Class I or Class II"
    )
    is_loh: bool = Field(
        False, description="True if this allele is lost via LOH in the tumor"
    )
    loh_pvalue: Optional[float] = Field(
        None, description="LOHHLA p-value for LOH call"
    )
    copy_number: Optional[float] = Field(
        None, description="Allele-specific copy number"
    )

    @property
    def gene(self) -> str:
        """Extract the HLA gene (e.g. HLA-A from HLA-A*02:01)."""
        return self.allele.split("*")[0]


class BindingPrediction(BaseModel):
    """MHC binding and presentation prediction for a single peptide-allele pair."""

    peptide: str = Field(
        ..., description="Peptide sequence"
    )
    hla_allele: str = Field(
        ..., description="HLA allele"
    )
    predictor: str = Field(
        ..., description="Tool used (netmhcpan, mhcflurry, netmhciipan, netmhcstabpan)"
    )

    # NetMHCpan-style outputs
    ic50: Optional[float] = Field(
        None, description="Predicted binding affinity (nM)"
    )
    percent_rank: Optional[float] = Field(
        None, description="Percent rank among random peptides"
    )
    binding_level: Optional[BindingLevel] = Field(
        None, description="SB/WB/NB classification"
    )

    # MHCflurry presentation score
    presentation_score: Optional[float] = Field(
        None, description="MHCflurry presentation score [0-1]"
    )

    # Stability
    half_life: Optional[float] = Field(
        None, description="Predicted pMHC half-life (hours)"
    )
    stability_rank: Optional[float] = Field(
        None, description="NetMHCstabpan stability percent rank"
    )


class MSEvidence(BaseModel):
    """Mass spectrometry evidence for a candidate neoantigen peptide."""

    peptide: str = Field(
        ..., description="Peptide sequence"
    )
    detected: bool = Field(
        False, description="Whether the peptide was detected by LC-MS/MS"
    )
    pepquery_score: Optional[float] = Field(
        None, description="PepQuery2 confidence score"
    )
    pepquery_pvalue: Optional[float] = Field(
        None, description="PepQuery2 p-value"
    )
    spectral_count: int = Field(
        0, description="Number of supporting PSMs"
    )
    best_search_engine_score: Optional[float] = Field(
        None, description="Best score from FragPipe/MSFragger search"
    )
    fdr: Optional[float] = Field(
        None, description="False discovery rate for this identification"
    )


class ExpressionData(BaseModel):
    """Gene and transcript-level expression for a variant."""

    gene_symbol: str = Field(
        ..., description="HUGO gene symbol"
    )
    transcript_id: Optional[str] = Field(
        None, description="Ensembl transcript ID"
    )
    gene_tpm: Optional[float] = Field(
        None, description="Gene-level TPM from Salmon"
    )
    transcript_tpm: Optional[float] = Field(
        None, description="Transcript-level TPM from Salmon"
    )
    vaf: Optional[float] = Field(
        None, description="Variant allele frequency in tumor DNA"
    )
    rna_vaf: Optional[float] = Field(
        None, description="Variant allele frequency in RNA-seq"
    )
    mutant_expression: Optional[float] = Field(
        None, description="Estimated mutant allele expression (VAF x TPM)"
    )


class AntigenProcessingFeatures(BaseModel):
    """Features related to antigen processing machinery."""

    proteasome_cleavage_score: Optional[float] = Field(
        None, description="C-terminal proteasomal cleavage prediction score"
    )
    tap_transport_score: Optional[float] = Field(
        None, description="TAP transport efficiency score"
    )
    b2m_mutated: bool = Field(
        False, description="Whether B2M is mutated in this tumor"
    )
    apm_gene_expression: Optional[dict[str, float]] = Field(
        None,
        description=(
            "Expression levels of key antigen processing genes "
            "(TAP1, TAP2, TAPBP, PSMB5, PSMB8, PSMB9, B2M, CALR, CANX, ERAP1, ERAP2)"
        )
    )


class NeoantigenCandidate(BaseModel):
    """
    A fully annotated neoantigen candidate -- the central data model.

    Aggregates genomic variant information, peptide sequences, binding
    predictions, expression, clonality, antigen processing, immune escape
    features, MS evidence, and the final ML-derived prioritization score.
    """

    # Identity
    candidate_id: str = Field(
        ..., description="Unique candidate identifier"
    )
    sample_id: str = Field(
        ..., description="Patient/sample identifier"
    )

    # Source variant
    variant: GenomicVariant = Field(
        ..., description="The underlying somatic variant"
    )

    # Peptide sequences
    mutant_peptide: str = Field(
        ..., description="Mutant peptide sequence"
    )
    wildtype_peptide: Optional[str] = Field(
        None, description="Corresponding wildtype peptide (None for fusions/ncORFs)"
    )
    peptide_length: int = Field(
        ..., description="Length of the mutant peptide"
    )
    peptide_source: VariantType = Field(
        ..., description="Source class of this peptide (SNV, fusion, splice, etc.)"
    )

    # HLA and binding
    hla_alleles: list[HLAAllele] = Field(
        default_factory=list, description="Patient HLA alleles"
    )
    binding_predictions: list[BindingPrediction] = Field(
        default_factory=list, description="All binding/presentation predictions"
    )
    best_binding: Optional[BindingPrediction] = Field(
        None, description="Best binding prediction across alleles and tools"
    )

    # Expression and clonality
    expression: Optional[ExpressionData] = Field(
        None, description="Expression and VAF data"
    )
    clonal_status: ClonalStatus = Field(
        ClonalStatus.INDETERMINATE, description="Clonal vs subclonal"
    )
    cancer_cell_fraction: Optional[float] = Field(
        None, description="Cancer cell fraction from PyClone-VI"
    )

    # Antigen processing
    processing_features: Optional[AntigenProcessingFeatures] = Field(
        None, description="Antigen processing machinery features"
    )

    # Immune escape
    presenting_allele_loh: bool = Field(
        False,
        description="True if the best-binding HLA allele has undergone LOH"
    )

    # Immunopeptidomics
    ms_evidence: Optional[MSEvidence] = Field(
        None, description="Mass spectrometry evidence"
    )

    # Sequence features
    self_similarity: Optional[float] = Field(
        None, description="BLOSUM62 similarity between mutant and wildtype peptide"
    )
    agretopicity: Optional[float] = Field(
        None, description="Ratio of mutant to wildtype binding affinity"
    )
    hydrophobicity: Optional[float] = Field(
        None, description="Kyte-Doolittle hydrophobicity of mutant peptide"
    )

    # ML scoring
    neoflux_score: Optional[float] = Field(
        None, description="Final NeoFlux ML prioritization score [0-1]"
    )
    neoflux_rank: Optional[int] = Field(
        None, description="Rank within sample (1 = best)"
    )
    shap_top_features: Optional[dict[str, float]] = Field(
        None, description="Top contributing features from SHAP analysis"
    )

    # Driver gene annotation
    is_driver_gene: bool = Field(
        False, description="Whether the gene is in a curated driver gene list"
    )

    def to_output_row(self) -> dict:
        """Flatten this candidate into a single dict for TSV/CSV output."""
        row = {
            "candidate_id": self.candidate_id,
            "sample_id": self.sample_id,
            "gene": self.variant.gene_symbol,
            "chromosome": self.variant.chromosome,
            "position": self.variant.position,
            "ref": self.variant.ref_allele,
            "alt": self.variant.alt_allele,
            "variant_type": self.variant.variant_type.value,
            "consequence": (
                self.variant.consequence.value if self.variant.consequence else None
            ),
            "protein_change": self.variant.protein_change,
            "mutant_peptide": self.mutant_peptide,
            "wildtype_peptide": self.wildtype_peptide,
            "peptide_length": self.peptide_length,
            "peptide_source": self.peptide_source.value,
        }

        # Best binding prediction
        if self.best_binding:
            row.update({
                "best_hla_allele": self.best_binding.hla_allele,
                "best_ic50": self.best_binding.ic50,
                "best_percent_rank": self.best_binding.percent_rank,
                "best_binding_level": (
                    self.best_binding.binding_level.value
                    if self.best_binding.binding_level else None
                ),
                "presentation_score": self.best_binding.presentation_score,
                "pMHC_half_life": self.best_binding.half_life,
            })

        # Expression
        if self.expression:
            row.update({
                "gene_tpm": self.expression.gene_tpm,
                "vaf": self.expression.vaf,
                "rna_vaf": self.expression.rna_vaf,
                "mutant_expression": self.expression.mutant_expression,
            })

        # Clonality
        row["clonal_status"] = self.clonal_status.value
        row["cancer_cell_fraction"] = self.cancer_cell_fraction

        # Immune escape
        row["presenting_allele_loh"] = self.presenting_allele_loh
        if self.processing_features:
            row["b2m_mutated"] = self.processing_features.b2m_mutated
            row["proteasome_score"] = (
                self.processing_features.proteasome_cleavage_score
            )
            row["tap_score"] = self.processing_features.tap_transport_score

        # MS evidence
        if self.ms_evidence:
            row.update({
                "ms_detected": self.ms_evidence.detected,
                "pepquery_score": self.ms_evidence.pepquery_score,
                "spectral_count": self.ms_evidence.spectral_count,
            })

        # Sequence features
        row["self_similarity"] = self.self_similarity
        row["agretopicity"] = self.agretopicity
        row["hydrophobicity"] = self.hydrophobicity
        row["is_driver_gene"] = self.is_driver_gene

        # ML score
        row["neoflux_score"] = self.neoflux_score
        row["neoflux_rank"] = self.neoflux_rank

        # SHAP
        if self.shap_top_features:
            for i, (feat, val) in enumerate(self.shap_top_features.items()):
                row[f"shap_feature_{i+1}"] = feat
                row[f"shap_value_{i+1}"] = val

        return row


class PatientSample(BaseModel):
    """A patient sample with all associated data paths and results."""

    sample_id: str = Field(
        ..., description="Unique patient/sample identifier"
    )
    tumor_bam: Optional[str] = Field(
        None, description="Path to tumor WES/WGS BAM/CRAM"
    )
    normal_bam: Optional[str] = Field(
        None, description="Path to matched normal BAM/CRAM"
    )
    rna_fastq_r1: Optional[str] = Field(
        None, description="Path to RNA-seq FASTQ read 1"
    )
    rna_fastq_r2: Optional[str] = Field(
        None, description="Path to RNA-seq FASTQ read 2"
    )
    rna_bam: Optional[str] = Field(
        None, description="Path to RNA-seq BAM (if pre-aligned)"
    )
    hla_file: Optional[str] = Field(
        None, description="Path to HLA allele list file (one allele per line)"
    )
    ms_mzml: Optional[str] = Field(
        None, description="Path to LC-MS/MS mzML file"
    )
    ms_search_results: Optional[str] = Field(
        None, description="Path to MS search results (FragPipe output)"
    )
    tumor_purity: Optional[float] = Field(
        None, description="Estimated tumor purity (0-1)"
    )
    tumor_ploidy: Optional[float] = Field(
        None, description="Estimated tumor ploidy"
    )

    # Results
    hla_alleles: list[HLAAllele] = Field(
        default_factory=list, description="Typed HLA alleles"
    )
    candidates: list[NeoantigenCandidate] = Field(
        default_factory=list, description="All neoantigen candidates"
    )


class CohortSummary(BaseModel):
    """Cohort-level aggregation of neoantigen results across patients."""

    cohort_name: str = Field(
        ..., description="Name of the cohort/study"
    )
    n_patients: int = Field(
        ..., description="Number of patients in the cohort"
    )
    samples: list[str] = Field(
        default_factory=list, description="Sample IDs in this cohort"
    )

    # Aggregate statistics
    total_candidates: int = Field(0, description="Total candidates across cohort")
    median_candidates_per_patient: Optional[float] = Field(None)
    variant_type_counts: Optional[dict[str, int]] = Field(
        None, description="Count of candidates by variant type"
    )
    ms_detection_rate: Optional[float] = Field(
        None, description="Fraction of top candidates confirmed by MS"
    )
    hla_loh_prevalence: Optional[float] = Field(
        None, description="Fraction of patients with at least one HLA LOH event"
    )
