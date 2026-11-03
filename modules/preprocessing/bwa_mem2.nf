// BWA-MEM2 alignment for WES/WGS reads
// Tool: bwa-mem2 v2.2+
// Container: quay.io/biocontainers/bwa-mem2

process BWA_MEM2 {
    tag "${meta.id}"
    label "process_high"

    input:
    tuple val(meta), path(reads_r1), path(reads_r2)
    path  genome_index

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam

    script:
    // BWA-MEM2 alignment followed by samtools sort and index.
    // Read group information is derived from the sample metadata.
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:lib1" \\
        ${genome_index}/genome \\
        ${reads_r1} ${reads_r2} \\
    | samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam -

    samtools index ${meta.id}.sorted.bam
    """

    stub:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai
    """
}
