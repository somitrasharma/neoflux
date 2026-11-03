// Variant annotation with Ensembl VEP.

process VEP {
    tag "${meta.id}"

    input:
    tuple val(meta), path(consensus_vcf), path(event_manifest)
    path genome_fasta
    path vep_cache

    output:
    tuple val(meta), path("${meta.id}.vep.vcf"), path(event_manifest), emit: vcf

    script:
    """
    vep \\
        --input_file ${consensus_vcf} \\
        --output_file ${meta.id}.vep.vcf \\
        --vcf \\
        --cache \\
        --dir_cache ${vep_cache} \\
        --fasta ${genome_fasta} \\
        --assembly ${params.genome} \\
        --symbol \\
        --canonical \\
        --protein \\
        --hgvs \\
        --force_overwrite
    """

    stub:
    """
    cp ${consensus_vcf} ${meta.id}.vep.vcf
    """
}
