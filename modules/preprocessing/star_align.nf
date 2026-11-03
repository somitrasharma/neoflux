// STAR 2-pass alignment for RNA-seq
// Tool: STAR v2.7+
// Container: trinityctat/starfusion (includes STAR)
//
// Two-pass alignment is used for improved splice junction detection,
// which is critical for the splice variant neoantigen module.

process STAR_ALIGN {
    tag "${meta.id}"
    label "process_high"

    input:
    tuple val(meta), path(reads_r1), path(reads_r2)
    path  star_index
    path  gtf

    output:
    tuple val(meta), path("${meta.id}.Aligned.sortedByCoord.out.bam"),
                     path("${meta.id}.Aligned.sortedByCoord.out.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.SJ.out.tab"),                        emit: splice_junctions
    tuple val(meta), path("${meta.id}.Chimeric.out.junction"),             emit: chimeric_junctions

    script:
    // STAR 2-pass mode with chimeric junction detection for STAR-Fusion.
    // --outSAMtype BAM SortedByCoordinate produces a coordinate-sorted BAM.
    // --chimSegmentMin 12 enables chimeric read detection for fusion calling.
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --sjdbGTFfile ${gtf} \\
        --readFilesIn ${reads_r1} ${reads_r2} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${meta.id}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --twopassMode Basic \\
        --chimSegmentMin 12 \\
        --chimJunctionOverhangMin 8 \\
        --chimOutType Junctions \\
        --quantMode GeneCounts

    samtools index ${meta.id}.Aligned.sortedByCoord.out.bam
    """

    stub:
    """
    touch ${meta.id}.Aligned.sortedByCoord.out.bam
    touch ${meta.id}.Aligned.sortedByCoord.out.bam.bai
    touch ${meta.id}.SJ.out.tab
    touch ${meta.id}.Chimeric.out.junction
    """
}
