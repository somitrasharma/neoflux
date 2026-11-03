// Gene fusion discovery with STAR-Fusion.

process STAR_FUSION {
    tag "${meta.id}"

    input:
    tuple val(meta), path(rna_bam), path(rna_bai)

    output:
    tuple val(meta), path("${meta.id}.star_fusion.tsv"), emit: fusions

    script:
    """
    samtools fastq ${rna_bam} > rna_${meta.id}.fastq
    STAR-Fusion \\
        --genome_lib_dir ${params.star_fusion_lib ?: "ctat_genome_lib"} \\
        --left_fq rna_${meta.id}.fastq \\
        --CPU ${task.cpus} \\
        --output_dir star_fusion_${meta.id}

    cp star_fusion_${meta.id}/star-fusion.fusion_predictions.tsv ${meta.id}.star_fusion.tsv
    """

    stub:
    """
    printf "#FusionName\\tJunctionReadCount\\tSpanningFragCount\\n" > ${meta.id}.star_fusion.tsv
    """
}
