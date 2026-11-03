// Extract tumor DNA/RNA VAF features from annotated variants.

process VAF_EXTRACT {
    tag "${meta.id}"

    input:
    tuple val(meta), path(tumor_bam), path(annotated_vcf), path(event_manifest)

    output:
    tuple val(meta), path("${meta.id}.vaf.tsv"), emit: vaf_table

    script:
    """
    python -m neoflux.cli extract-vaf \\
        --sample-id ${meta.id} \\
        --tumor-bam ${tumor_bam} \\
        --annotated-vcf ${annotated_vcf} \\
        --output ${meta.id}.vaf.tsv
    """

    stub:
    """
    printf "variant_id\\tvaf\\tdepth\\trna_vaf\\n" > ${meta.id}.vaf.tsv
    """
}
