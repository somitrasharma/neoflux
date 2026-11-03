// High-resolution HLA typing with HLA-HD.
// This module is a scaffold for centers that prefer HLA-HD over OptiType.

process HLA_HD {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads_r1), path(reads_r2)

    output:
    tuple val(meta), path("${meta.id}.hla.tsv"), emit: alleles

    script:
    """
    hlahd.sh \\
        -t ${task.cpus} \\
        -m 100 \\
        -c 0.95 \\
        ${reads_r1} ${reads_r2} \\
        HLA_gene.split.txt \\
        dictionary \\
        ${meta.id} \\
        hlahd_${meta.id}

    python -m neoflux.cli parse-hla \\
        --input hlahd_${meta.id}/${meta.id}/result/${meta.id}_final.result.txt \\
        --format hlahd \\
        --output ${meta.id}.hla.tsv
    """

    stub:
    """
    printf "allele\\thla_class\\nHLA-A*02:01\\tI\\nHLA-B*07:02\\tI\\nHLA-C*07:02\\tI\\n" > ${meta.id}.hla.tsv
    """
}
