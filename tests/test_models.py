from neoflux.ml_ranker import heuristic_score_candidates
from neoflux.models import GenomicVariant, NeoantigenCandidate, VariantType


def test_heuristic_ranker_assigns_scores_and_ranks():
    variant = GenomicVariant(
        sample_id="S1",
        chromosome="chr1",
        position=100,
        ref_allele="A",
        alt_allele="T",
        variant_type=VariantType.SNV,
    )
    candidate = NeoantigenCandidate(
        candidate_id="NF-test",
        sample_id="S1",
        variant=variant,
        mutant_peptide="AAAAAAAAA",
        peptide_length=9,
        peptide_source=VariantType.SNV,
    )

    ranked = heuristic_score_candidates([candidate])

    assert ranked[0].neoflux_rank == 1
    assert ranked[0].neoflux_score is not None
    assert 0.0 <= ranked[0].neoflux_score <= 1.0
