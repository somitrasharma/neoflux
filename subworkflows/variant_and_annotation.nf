// Variant calling, unconventional event discovery, VEP annotation, and peptide
// generation.

include { MUTECT2            } from '../modules/variant_calling/mutect2'
include { STRELKA2           } from '../modules/variant_calling/strelka2'
include { STAR_FUSION        } from '../modules/variant_calling/star_fusion'
include { RMATS              } from '../modules/variant_calling/rmats'
include { CONSENSUS_MERGE    } from '../modules/variant_calling/consensus_merge'
include { VEP                } from '../modules/annotation/vep'
include { PEPTIDE_GENERATION } from '../modules/annotation/peptide_generation'

workflow VARIANT_AND_ANNOTATION {
    take:
    ch_tumor_bam
    ch_normal_bam
    ch_rna_bam
    genome_fasta
    vep_cache
    gtf

    main:
    ch_tn_pairs = ch_tumor_bam.join(ch_normal_bam, by: 0)

    MUTECT2(ch_tn_pairs, genome_fasta)
    STRELKA2(ch_tn_pairs, genome_fasta)
    STAR_FUSION(ch_rna_bam)
    RMATS(ch_rna_bam, gtf)

    ch_mutect_vcf = MUTECT2.out.vcf.map { meta, vcf, tbi -> [meta, vcf] }
    ch_strelka_vcf = STRELKA2.out.vcf.map { meta, vcf, tbi -> [meta, vcf] }
    ch_small_variants = ch_mutect_vcf.join(ch_strelka_vcf, by: 0)
    ch_rna_events = STAR_FUSION.out.fusions.join(RMATS.out.events, by: 0)
    ch_all_events = ch_small_variants.join(ch_rna_events, by: 0)

    CONSENSUS_MERGE(ch_all_events)
    VEP(CONSENSUS_MERGE.out.merged, genome_fasta, vep_cache)
    PEPTIDE_GENERATION(VEP.out.vcf, gtf)

    emit:
    vcf_annotated = VEP.out.vcf
    peptides = PEPTIDE_GENERATION.out.peptides
}
