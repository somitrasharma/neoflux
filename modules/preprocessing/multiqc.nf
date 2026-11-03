process MULTIQC {
    tag "multiqc"

    // RESTORE THE INPUTS: collect all FastQC HTML and ZIP files (and any other QC logs)
    input:
    path qc_files

    output:
    path "results/multiqc/multiqc_report.html", emit: report
    path "results/multiqc/multiqc_data", emit: data

    script:
    """
    mkdir -p results/multiqc
    # MultiQC will scan the staged qc_files (which are symlinked into the work dir)
    multiqc . --force --title "NeoFlux: ${workflow.runName}" --outdir results/multiqc
    """

    stub:
    """
    mkdir -p results/multiqc
    printf "<html><body><h1>MultiQC Report (stub)</h1></body></html>" > results/multiqc/multiqc_report.html
    mkdir -p results/multiqc/multiqc_data
    touch results/multiqc/multiqc_data/multiqc.log
    """
}
