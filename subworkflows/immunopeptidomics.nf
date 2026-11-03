// Optional immunopeptidomics integration.

include { FRAGPIPE  } from '../modules/immunopeptidomics/fragpipe'
include { PEPQUERY2 } from '../modules/immunopeptidomics/pepquery2'

workflow IMMUNOPEPTIDOMICS {
    take:
    ch_ms_data
    ch_peptides

    main:
    FRAGPIPE(ch_ms_data)

    ch_pepquery_input = FRAGPIPE.out.search_results.join(ch_peptides, by: 0)
    PEPQUERY2(ch_pepquery_input)

    emit:
    evidence = PEPQUERY2.out.evidence
}
