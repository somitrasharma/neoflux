// Targeted peptide validation with PepQuery2.

process PEPQUERY2 {
    tag "${meta.id}"

    input:
    tuple val(meta), path(fragpipe_peptides), path(candidates_jsonl), path(peptides_fasta), path(peptides_tsv)

    output:
    tuple val(meta), path("${meta.id}.ms_evidence.tsv"), emit: evidence

    script:
    """
    pepquery \\
        -b ${peptides_fasta} \\
        -db reference_proteome.fasta \\
        -o pepquery_${meta.id} \\
        -cpu ${task.cpus}

    python -m neoflux.cli score-ms \\
        --candidates ${candidates_jsonl} \\
        --fragpipe ${fragpipe_peptides} \\
        --pepquery pepquery_${meta.id}/psm_rank.txt \\
        --output ${meta.id}.ms_evidence.tsv
    """

    stub:
    """
    printf "candidate_id\\tpeptide\\tdetected\\tpepquery_score\\tpepquery_pvalue\\tspectral_count\\n" > ${meta.id}.ms_evidence.tsv
    """
}
