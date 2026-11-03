// Somatic SNV/indel calling with Strelka2.

process STRELKA2 {
    tag "${meta.id}"

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam)
    path genome_fasta

    output:
    tuple val(meta), path("${meta.id}.strelka2.somatic.vcf.gz"), path("${meta.id}.strelka2.somatic.vcf.gz.tbi"), emit: vcf

    script:
    """
    configureStrelkaSomaticWorkflow.py \\
        --tumorBam ${tumor_bam} \\
        --normalBam ${normal_bam} \\
        --referenceFasta ${genome_fasta} \\
        --runDir strelka_${meta.id}

    strelka_${meta.id}/runWorkflow.py -m local -j ${task.cpus}
    cp strelka_${meta.id}/results/variants/somatic.snvs.vcf.gz ${meta.id}.strelka2.somatic.vcf.gz
    tabix -p vcf ${meta.id}.strelka2.somatic.vcf.gz
    """

    stub:
    """
    printf "##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\nchr7\\t140453136\\t.\\tA\\tT\\t.\\tPASS\\tGENE=BRAF;HGVS_P=p.V600E;NEOFLUX_CALLERS=strelka2\\n" > ${meta.id}.strelka2.somatic.vcf.gz
    touch ${meta.id}.strelka2.somatic.vcf.gz.tbi
    """
}
