// Patient and cohort reporting subworkflow.

include { PATIENT_REPORT } from '../modules/reporting/patient_report'
include { COHORT_REPORT  } from '../modules/reporting/cohort_report'

workflow RANKING_AND_REPORTING {
    take:
    ch_scored_candidates

    main:
    PATIENT_REPORT(ch_scored_candidates)
    COHORT_REPORT(PATIENT_REPORT.out.ranked_table.map { meta, table -> table }.collect())

    emit:
    ranked_tables = PATIENT_REPORT.out.ranked_table
    patient_html = PATIENT_REPORT.out.html
    cohort_summary = COHORT_REPORT.out.summary
    cohort_html = COHORT_REPORT.out.html
    cohort_figures = COHORT_REPORT.out.figures
}
