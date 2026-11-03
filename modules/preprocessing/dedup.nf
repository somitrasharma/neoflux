// BAM deduplication using sambamba
// Tool: sambamba v1.0+
// Container: quay.io/biocontainers/sambamba

process DEDUP {
    tag "${meta.id}"
    label "process_medium"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.dedup.bam"), path("${meta.id}.dedup.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.dedup.metrics.txt"),                            emit: metrics

    script:
    """
    sambamba markdup \\
        --nthreads ${task.cpus} \\
        --tmpdir tmp_dedup \\
        ${bam} \\
        ${meta.id}.dedup.bam \\
        2> ${meta.id}.dedup.metrics.txt

    sambamba index ${meta.id}.dedup.bam
    """

    stub:
    """
    touch ${meta.id}.dedup.bam ${meta.id}.dedup.bam.bai ${meta.id}.dedup.metrics.txt
    """
}
