// HLA typing and immune escape subworkflow.
// Binding itself is run in main.nf after peptides exist; this subworkflow emits
// normalized HLA alleles and LOHHLA-compatible immune escape calls.

include { OPTITYPE } from '../modules/hla/optitype'
include { LOHHLA   } from '../modules/hla/lohhla'
// HLA_HD is available at modules/hla/hla_hd.nf as an alternative to OPTITYPE

workflow HLA_AND_BINDING {
    take:
    ch_tumor_bam
    ch_normal_bam
    ch_hla_file
    run_hla_typing

    main:
    if (run_hla_typing) {
        OPTITYPE(ch_tumor_bam)
        ch_hla_alleles = OPTITYPE.out.alleles
    } else {
        ch_hla_alleles = ch_hla_file
    }

    ch_loh_input = ch_tumor_bam
        .join(ch_normal_bam, by: 0)
        .join(ch_hla_alleles, by: 0)

    LOHHLA(ch_loh_input)

    emit:
    hla_alleles = ch_hla_alleles
    loh_calls = LOHHLA.out.loh_calls
}
