
from neoflux.features import FEATURE_COLUMNS, FeatureEngineer
from neoflux.models import (
    BindingLevel,
    BindingPrediction,
    ClonalStatus,
    ExpressionData,
    GenomicVariant,
    MSEvidence,
    NeoantigenCandidate,
    VariantType,
)


def _candidate() -> NeoantigenCandidate:
    variant = GenomicVariant(
        sample_id="S1",
        chromosome="chr7",
        position=140453136,
        ref_allele="A",
        alt_allele="T",
        variant_type=VariantType.SNV,
        gene_symbol="BRAF",
    )
    best_binding = BindingPrediction(
        peptide="GLYVGLG",
        hla_allele="HLA-A*02:01",
        predictor="mhcflurry",
        ic50=120.0,
        percent_rank=0.8,
        binding_level=BindingLevel.STRONG_BINDER,
        presentation_score=0.75,
        half_life=4.0,
        stability_rank=1.2,
    )
    wt_binding = BindingPrediction(
        peptide="GLYVGLA",
        hla_allele="HLA-A*02:01",
        predictor="mhcflurry",
        ic50=2400.0,
        percent_rank=12.0,
        binding_level=BindingLevel.NON_BINDER,
    )
    return NeoantigenCandidate(
        candidate_id="NF-feature-test",
        sample_id="S1",
        variant=variant,
        mutant_peptide="GLYVGLG",
        wildtype_peptide="GLYVGLA",
        peptide_length=7,
        peptide_source=VariantType.SNV,
        binding_predictions=[best_binding, wt_binding],
        best_binding=best_binding,
        expression=ExpressionData(
            gene_symbol="BRAF",
            gene_tpm=40.0,
            transcript_tpm=20.0,
            vaf=0.35,
            rna_vaf=0.3,
            mutant_expression=12.0,
        ),
        clonal_status=ClonalStatus.CLONAL,
        cancer_cell_fraction=0.7,
        ms_evidence=MSEvidence(peptide="GLYVGLG", detected=True, spectral_count=2),
        presenting_allele_loh=False,
    )


def test_feature_matrix_shape_and_columns():
    engineer = FeatureEngineer()
    matrix = engineer.build_feature_matrix([_candidate()])

    assert matrix.shape == (1, len(FEATURE_COLUMNS))
    assert list(matrix.columns) == FEATURE_COLUMNS


def test_enrich_candidates_sets_driver_and_agretopicity():
    engineer = FeatureEngineer()
    candidate = _candidate()
    enriched = engineer.enrich_candidates([candidate])[0]

    assert enriched.is_driver_gene is True
    assert enriched.hydrophobicity is not None
    assert enriched.self_similarity is not None
    assert enriched.agretopicity == 20.0

    row = engineer.extract_features(enriched)
    assert row["is_snv"] == 1.0
    assert row["ms_detected"] == 1.0
    assert row["is_clonal"] == 1.0
