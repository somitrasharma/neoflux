// Cohort-level summaries, metrics, and SITC-style figure tables.

process COHORT_REPORT {
    tag "cohort"

    input:
    path ranked_tables

    output:
    path "cohort_summary.tsv", emit: summary
    path "cohort_report.html", emit: html
    path "figures", emit: figures

    script:
    def ranked_args = ranked_tables instanceof List ? ranked_tables.join(' ') : ranked_tables
    """
    mkdir -p figures
    python -m neoflux.cli cohort-report \\
        --ranked-tables ${ranked_args} \\
        --summary cohort_summary.tsv \\
        --html cohort_report.html \\
        --figures-dir figures \\
        --cohort-name NeoFlux_Cohort
    """

    stub:
    """
    mkdir -p figures
    printf "metric\\tvalue\\n" > cohort_summary.tsv
    printf "<html><body><h1>NeoFlux Cohort Report</h1></body></html>\\n" > cohort_report.html
    """
}
