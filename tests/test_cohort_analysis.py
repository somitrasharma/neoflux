from neoflux.cohort_analysis import CohortAnalyzer
from neoflux.models import (
    GenomicVariant,
    MSEvidence,
    NeoantigenCandidate,
    PatientSample,
    VariantType,
)


def _patient(sample_id: str, score: float, ms_detected: bool) -> PatientSample:
    variant = GenomicVariant(
        sample_id=sample_id,
        chromosome="chr1",
        position=100,
        ref_allele="A",
        alt_allele="T",
        variant_type=VariantType.SNV,
        gene_symbol="BRAF",
    )
    candidate = NeoantigenCandidate(
        candidate_id=f"{sample_id}-1",
        sample_id=sample_id,
        variant=variant,
        mutant_peptide="GLYVGLG",
        peptide_length=7,
        peptide_source=VariantType.SNV,
        neoflux_score=score,
        neoflux_rank=1,
        ms_evidence=MSEvidence(peptide="GLYVGLG", detected=ms_detected),
    )
    return PatientSample(sample_id=sample_id, candidates=[candidate])


def test_cohort_summary_and_correlation():
    analyzer = CohortAnalyzer(cohort_name="demo")
    analyzer.add_patient(_patient("S1", 0.9, True))
    analyzer.add_patient(_patient("S2", 0.2, False))

    summary = analyzer.compute_summary()
    assert summary.n_patients == 2
    assert summary.total_candidates == 2
    assert summary.variant_type_counts["SNV"] == 2

    correlation = analyzer.compute_score_vs_ms_correlation()
    assert not correlation.empty
    assert "detection_rate" in correlation.columns

    burden = analyzer.compute_neoantigen_burden()
    assert burden["total"].sum() == 2
