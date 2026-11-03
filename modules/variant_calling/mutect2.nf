// Somatic SNV/indel calling with GATK Mutect2.

process MUTECT2 {
    tag "${meta.id}"

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam)
    path genome_fasta

    output:
    tuple val(meta), path("${meta.id}.mutect2.filtered.vcf.gz"), path("${meta.id}.mutect2.filtered.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk Mutect2 \\
        -R ${genome_fasta} \\
        -I ${tumor_bam} \\
        -I ${normal_bam} \\
        -tumor ${meta.id}_tumor \\
        -normal ${meta.id}_normal \\
        -O ${meta.id}.mutect2.raw.vcf.gz

    gatk FilterMutectCalls \\
        -R ${genome_fasta} \\
        -V ${meta.id}.mutect2.raw.vcf.gz \\
        -O ${meta.id}.mutect2.filtered.vcf.gz

    tabix -p vcf ${meta.id}.mutect2.filtered.vcf.gz
    """

    stub:
    """
    printf "##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\nchr7\\t140453136\\t.\\tA\\tT\\t.\\tPASS\\tGENE=BRAF;HGVS_P=p.V600E;NEOFLUX_CALLERS=mutect2\\n" > ${meta.id}.mutect2.filtered.vcf.gz
    touch ${meta.id}.mutect2.filtered.vcf.gz.tbi
    """
}
