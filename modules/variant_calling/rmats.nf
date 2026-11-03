// Aberrant splicing event discovery with rMATS.
// For a production tumor/normal comparison, provide matched normal or cohort
// control BAMs. This first scaffold records tumor RNA splice events.

process RMATS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(rna_bam), path(rna_bai)
    path gtf

    output:
    tuple val(meta), path("${meta.id}.rmats_events.tsv"), emit: events

    script:
    """
    printf "${rna_bam}\\n" > ${meta.id}.tumor_bams.txt

    rmats.py \\
        --b1 ${meta.id}.tumor_bams.txt \\
        --gtf ${gtf} \\
        --od rmats_${meta.id} \\
        --tmp rmats_tmp_${meta.id} \\
        -t paired \\
        --readLength 101 \\
        --nthread ${task.cpus}

    python -m neoflux.cli merge-splice-events \\
        --rmats-dir rmats_${meta.id} \\
        --output ${meta.id}.rmats_events.tsv
    """

    stub:
    """
    printf "event_id\\tevent_type\\tgene_symbol\\tdelta_psi\\n" > ${meta.id}.rmats_events.tsv
    """
}
