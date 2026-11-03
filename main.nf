#!/usr/bin/env nextflow

// NeoFlux: Unified Neoantigen Discovery and Prioritization Pipeline
// SITC 2026 -- Neoantigen Identification and Characterization
//
// A modular, end-to-end pipeline integrating genomics, transcriptomics,
// immunopeptidomics, and ML-based ranking for solid tumor neoantigen
// discovery.

nextflow.enable.dsl = 2

// -- Include modules --

include { FASTQC             } from './modules/preprocessing/fastqc'
include { MULTIQC            } from './modules/preprocessing/multiqc'
include { STAR_ALIGN         } from './modules/preprocessing/star_align'
include { DEDUP              } from './modules/preprocessing/dedup'
include { NETMHCPAN          } from './modules/binding/netmhcpan'
include { MHCFLURRY          } from './modules/binding/mhcflurry'
include { NETMHCIIPAN        } from './modules/binding/netmhciipan'
include { NETMHCSTABPAN      } from './modules/binding/netmhcstabpan'
include { SALMON_QUANT       } from './modules/expression/salmon_quant'
include { VAF_EXTRACT        } from './modules/expression/vaf_extract'
include { PYCLONE            } from './modules/expression/pyclone'
include { FEATURE_ENGINEERING} from './modules/ml_ranking/feature_engineering'
include { ML_SCORE           } from './modules/ml_ranking/ml_score'

// -- Include subworkflows --

include { VARIANT_AND_ANNOTATION } from './subworkflows/variant_and_annotation'
include { HLA_AND_BINDING        } from './subworkflows/hla_and_binding'
include { IMMUNOPEPTIDOMICS      } from './subworkflows/immunopeptidomics'
include { RANKING_AND_REPORTING  } from './subworkflows/ranking_and_reporting'


// -- Helper: parse samplesheet --

def parse_samplesheet(samplesheet_path) {
    Channel
        .fromPath(samplesheet_path)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id:           row.sample_id,
                tumor_purity: row.tumor_purity ? row.tumor_purity as float : null,
                tumor_ploidy: row.tumor_ploidy ? row.tumor_ploidy as float : null,
            ]
            def tumor_bam    = row.tumor_bam    ? file(row.tumor_bam)    : null
            def normal_bam   = row.normal_bam   ? file(row.normal_bam)   : null
            def rna_r1       = row.rna_fastq_r1 ? file(row.rna_fastq_r1) : null
            def rna_r2       = row.rna_fastq_r2 ? file(row.rna_fastq_r2) : null
            def hla_file     = row.hla_file     ? file(row.hla_file)     : null
            def ms_mzml      = row.ms_mzml      ? file(row.ms_mzml)      : null
            def ms_results   = row.ms_search_results ? file(row.ms_search_results) : file("${projectDir}/assets/empty_ms_search.tsv")
            return [meta, tumor_bam, normal_bam, rna_r1, rna_r2, hla_file, ms_mzml, ms_results]
        }
}


// ============================================================
//  MAIN WORKFLOW
// ============================================================

workflow {

    // Validate required parameters
    if (!params.samplesheet) {
        error "Please provide a samplesheet via --samplesheet"
    }

    // Parse samplesheet. DSL2 auto-forks channels for multiple consumers.
    ch_samplesheet = parse_samplesheet(params.samplesheet)

    ch_tumor_bam  = ch_samplesheet.map { meta, tb, nb, r1, r2, hla, ms, msr -> [meta, tb] }.filter { it[1] != null }
    ch_normal_bam = ch_samplesheet.map { meta, tb, nb, r1, r2, hla, ms, msr -> [meta, nb] }.filter { it[1] != null }
    ch_rna_reads  = ch_samplesheet.map { meta, tb, nb, r1, r2, hla, ms, msr -> [meta, r1, r2] }.filter { it[1] != null }
    ch_hla_file   = ch_samplesheet.map { meta, tb, nb, r1, r2, hla, ms, msr -> [meta, hla] }.filter { it[1] != null }
    ch_ms_data    = ch_samplesheet.map { meta, tb, nb, r1, r2, hla, ms, msr -> [meta, ms, msr] }.filter { it[1] != null }

    // -----------------------------------------------------------
    //  MODULE 1: PREPROCESSING
    // -----------------------------------------------------------

    if (!params.skip_qc) {
        FASTQC(ch_rna_reads)
        ch_qc_html = FASTQC.out.html.map { meta, report -> report }.collect()
        ch_qc_zip  = FASTQC.out.zip.map { meta, archive -> archive }.collect()
        MULTIQC(ch_qc_html.mix(ch_qc_zip))
    }

    STAR_ALIGN(ch_rna_reads, params.star_index, params.genome_gtf)
    ch_rna_bam = STAR_ALIGN.out.bam

    DEDUP(ch_tumor_bam.mix(ch_rna_bam.map { meta, bam, bai -> [meta, bam] }))

    // -----------------------------------------------------------
    //  MODULE 2: HLA TYPING + IMMUNE ESCAPE
    // -----------------------------------------------------------

    HLA_AND_BINDING(ch_tumor_bam, ch_normal_bam, ch_hla_file, params.run_hla_typing)

    // -----------------------------------------------------------
    //  MODULE 3-4: VARIANT CALLING + ANNOTATION + PEPTIDE GEN
    // -----------------------------------------------------------

    VARIANT_AND_ANNOTATION(
        ch_tumor_bam, ch_normal_bam, ch_rna_bam,
        params.genome_fasta, params.vep_cache_dir, params.genome_gtf
    )

    // DSL2 auto-forks process outputs, so VARIANT_AND_ANNOTATION.out.peptides
    // and .out.vcf_annotated can be consumed by multiple downstream operators.
    ch_peptides = VARIANT_AND_ANNOTATION.out.peptides
    ch_vcf_ann  = VARIANT_AND_ANNOTATION.out.vcf_annotated

    // -----------------------------------------------------------
    //  MODULE 5: BINDING PREDICTION
    // -----------------------------------------------------------

    ch_hla = HLA_AND_BINDING.out.hla_alleles

    NETMHCPAN(ch_peptides.combine(ch_hla, by: 0))
    MHCFLURRY(ch_peptides.combine(ch_hla, by: 0))
    NETMHCSTABPAN(ch_peptides.combine(ch_hla, by: 0))

    if (params.run_class_ii) {
        NETMHCIIPAN(ch_peptides.combine(ch_hla, by: 0))
    }

    // -----------------------------------------------------------
    //  MODULE 6: EXPRESSION + CLONALITY
    // -----------------------------------------------------------

    SALMON_QUANT(ch_rna_reads, params.salmon_index)
    VAF_EXTRACT(ch_tumor_bam.join(ch_vcf_ann, by: 0))
    PYCLONE(VAF_EXTRACT.out.vaf_table)

    // -----------------------------------------------------------
    //  MODULE 7: IMMUNOPEPTIDOMICS (optional)
    // -----------------------------------------------------------

    if (params.run_ms) {
        IMMUNOPEPTIDOMICS(ch_ms_data, ch_peptides)
    }

    // -----------------------------------------------------------
    //  MODULE 8: FEATURE ENGINEERING + ML RANKING
    // -----------------------------------------------------------

    ch_binding_all = NETMHCPAN.out.predictions
        .mix(MHCFLURRY.out.predictions)
        .mix(NETMHCSTABPAN.out.predictions)

    if (params.run_class_ii) {
        ch_binding_all = ch_binding_all.mix(NETMHCIIPAN.out.predictions)
    }

    ch_binding = ch_binding_all.groupTuple(by: 0)

    ch_ms_evidence = params.run_ms
        ? IMMUNOPEPTIDOMICS.out.evidence
        : ch_peptides.map { meta, candidates_jsonl, peptides_fasta, peptides_tsv ->
            [meta, file("${projectDir}/assets/empty_ms_evidence.tsv")]
        }

    ch_feature_input = ch_peptides
        .join(ch_binding, by: 0)
        .join(SALMON_QUANT.out.quant, by: 0)
        .join(PYCLONE.out.clusters, by: 0)
        .join(HLA_AND_BINDING.out.loh_calls, by: 0)
        .join(ch_ms_evidence, by: 0)
        .join(ch_vcf_ann, by: 0)

    FEATURE_ENGINEERING(ch_feature_input)

    ML_SCORE(
        FEATURE_ENGINEERING.out.feature_matrix,
        params.model_path
    )

    // -----------------------------------------------------------
    //  MODULE 9: REPORTING
    // -----------------------------------------------------------

    RANKING_AND_REPORTING(
        ML_SCORE.out.scored_candidates
    )
}



