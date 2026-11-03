// FastQC quality control for FASTQ/BAM inputs
// Tool: FastQC v0.12+
// Container: biocontainers/fastqc:v0.12.1_cv2

process FASTQC {
    tag "${meta.id}"
    label "process_low"

    input:
    tuple val(meta), path(reads_r1), path(reads_r2)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${reads_r1} ${reads_r2}
    """

    stub:
    """
    touch ${meta.id}_R1_fastqc.html ${meta.id}_R1_fastqc.zip
    touch ${meta.id}_R2_fastqc.html ${meta.id}_R2_fastqc.zip
    """
}
