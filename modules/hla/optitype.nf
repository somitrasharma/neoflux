// HLA class I typing with OptiType.
// Input can be tumor DNA/RNA BAM. For production, prefer RNA-seq or WES reads
// extracted to FASTQ according to the OptiType documentation.

process OPTITYPE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(tumor_bam)

    output:
    tuple val(meta), path("${meta.id}.hla.tsv"), emit: alleles

    script:
    """
    OptiTypePipeline.py \\
        -i ${tumor_bam} \\
        --dna \\
        -p ${meta.id} \\
        -o optitype_${meta.id}

    python -m neoflux.cli parse-hla \\
        --input optitype_${meta.id}/${meta.id}_result.tsv \\
        --format optitype \\
        --output ${meta.id}.hla.tsv
    """

    stub:
    """
    printf "allele\\thla_class\\nHLA-A*02:01\\tI\\nHLA-B*07:02\\tI\\nHLA-C*07:02\\tI\\n" > ${meta.id}.hla.tsv
    """
}
