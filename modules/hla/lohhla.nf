// HLA loss of heterozygosity with LOHHLA.
// Expects tumor BAM, matched normal BAM, and normalized HLA allele TSV.
//
// Production use (requires R + LOHHLA):
//   conda install -c bioconda bioconductor-lohhla
//   git clone https://github.com/ANHIG/LOHHLA.git
//
// Without R/LOHHLA, outputs all-False placeholder.

process LOHHLA {
    tag "${meta.id}"

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(hla_tsv)

    output:
    tuple val(meta), path("${meta.id}.hla_loh.tsv"), emit: loh_calls

    script:
    """
    # Coverage-based LOH detection using tumor/normal BAM comparison.
    # Replaces all-False placeholder with biologically motivated calls.
    python -m neoflux.cli detect-loh \\
        --tumor-bam ${tumor_bam} \\
        --normal-bam ${normal_bam} \\
        --hla ${hla_tsv} \\
        --output ${meta.id}.hla_loh.tsv

    # Production: replace with LOHHLA R script for higher accuracy:
    # if command -v Rscript &> /dev/null && [ -f "${projectDir}/scripts/lohhla_wrapper.R" ]; then
    #     Rscript ${projectDir}/scripts/lohhla_wrapper.R \\
    #         --tumor_bam ${tumor_bam} \\
    #         --normal_bam ${normal_bam} \\
    #         --hla_tsv ${hla_tsv} \\
    #         --sample_id ${meta.id} \\
    #         --output ${meta.id}.hla_loh.tsv
    # fi
    """

    stub:
    """
    printf "allele\\tis_loh\\tloh_pvalue\\tcopy_number\\n" > ${meta.id}.hla_loh.tsv
    """
}
