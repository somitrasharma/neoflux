// Per-patient ranked output table and lightweight HTML report.

process PATIENT_REPORT {
    tag "${meta.id}"

    input:
    tuple val(meta), path(scored_candidates)

    output:
    tuple val(meta), path("${meta.id}.ranked_neoantigens.tsv"), emit: ranked_table
    tuple val(meta), path("${meta.id}.patient_report.html"), emit: html

    script:
    """
    python -m neoflux.cli patient-report \\
        --sample-id ${meta.id} \\
        --candidates ${scored_candidates} \\
        --ranked-tsv ${meta.id}.ranked_neoantigens.tsv \\
        --html ${meta.id}.patient_report.html
    """

    stub:
    """
    printf "candidate_id\\tsample_id\\tneoflux_score\\tneoflux_rank\\n" > ${meta.id}.ranked_neoantigens.tsv
    printf "<html><body><h1>%s</h1></body></html>\\n" "${meta.id}" > ${meta.id}.patient_report.html
    """
}
