// RNA expression quantification with Salmon.

process SALMON_QUANT {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads_r1), path(reads_r2)
    path salmon_index

    output:
    tuple val(meta), path("${meta.id}.quant.sf"), emit: quant

    script:
    """
    salmon quant \\
        -i ${salmon_index} \\
        -l A \\
        -1 ${reads_r1} \\
        -2 ${reads_r2} \\
        -p ${task.cpus} \\
        --validateMappings \\
        -o salmon_${meta.id}

    cp salmon_${meta.id}/quant.sf ${meta.id}.quant.sf
    """

    stub:
    """
    printf "Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\n" > ${meta.id}.quant.sf
    """
}
