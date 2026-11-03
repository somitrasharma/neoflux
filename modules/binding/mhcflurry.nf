// MHC class I binding + presentation prediction with MHCflurry Python API.
// Uses neoflux predict-mhcflurry command for direct Python API access.

process MHCFLURRY {
    tag "${meta.id}"

    input:
    tuple val(meta), path(candidates_jsonl), path(peptides_fasta), path(peptides_tsv), path(hla_tsv)

    output:
    tuple val(meta), path("${meta.id}.mhcflurry.tsv"), emit: predictions

    script:
    """
    python -m neoflux.cli predict-mhcflurry \\
        --peptides ${peptides_tsv} \\
        --hla ${hla_tsv} \\
        --output ${meta.id}.mhcflurry.tsv
    """

    stub:
    """
    printf "candidate_id\\tpeptide\\thla_allele\\tpredictor\\tic50\\tpercent_rank\\tbinding_level\\tpresentation_score\\thalf_life\\tstability_rank\\n" > ${meta.id}.mhcflurry.tsv
    """
}
