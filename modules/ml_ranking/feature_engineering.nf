// Assemble candidate-level ML feature matrix.

process FEATURE_ENGINEERING {
    tag "${meta.id}"

    input:
    tuple val(meta), path(candidates_jsonl), path(peptides_fasta), path(peptides_tsv), path(binding_predictions), path(expression_quant), path(clonality), path(loh_calls), path(ms_evidence), path(annotated_vcf), path(event_manifest)

    output:
    tuple val(meta), path("${meta.id}.features.tsv"), path("${meta.id}.features.candidates.jsonl"), emit: feature_matrix

    script:
    def binding_args = binding_predictions instanceof List ? binding_predictions.join(' ') : binding_predictions
    """
    python -m neoflux.cli make-features \\
        --candidates ${candidates_jsonl} \\
        --binding ${binding_args} \\
        --expression ${expression_quant} \\
        --clonality ${clonality} \\
        --loh ${loh_calls} \\
        --ms-evidence ${ms_evidence} \\
        --annotated-vcf ${annotated_vcf} \\
        --features ${meta.id}.features.tsv \\
        --enriched-candidates ${meta.id}.features.candidates.jsonl
    """

    stub:
    """
    printf "candidate_id\\tsample_id\\tbest_ic50\\tbest_percent_rank\\tneoflux_score\\n" > ${meta.id}.features.tsv
    cp ${candidates_jsonl} ${meta.id}.features.candidates.jsonl
    """
}
