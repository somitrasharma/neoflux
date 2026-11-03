#!/usr/bin/env Rscript
# LOHHLA wrapper for NeoFlux pipeline.
#
# Usage:
#   Rscript lohhla_wrapper.R \
#     --tumor_bam <path> \
#     --normal_bam <path> \
#     --hla_tsv <path> \
#     --sample_id <id> \
#     --output <path>
#
# Requirements:
#   - R >= 4.0 with data.table and optparse
#   - LOHHLA installed: https://github.com/ANHIG/LOHHLA
#   - samtools in PATH
#   - bedtools in PATH
#
# Installation (conda):
#   conda install -c bioconda bioconductor-lohhla
#
# Alternatively, clone and source manually:
#   git clone https://github.com/ANHIG/LOHHLA.git
#   # then set LOHHLA_DIR below

library(optparse)
library(data.table)

option_list <- list(
  make_option("--tumor_bam", type="character", help="Tumor BAM file"),
  make_option("--normal_bam", type="character", help="Matched normal BAM file"),
  make_option("--hla_tsv", type="character", help="HLA allele TSV (allele column)"),
  make_option("--sample_id", type="character", default="sample", help="Sample ID"),
  make_option("--output", type="character", help="Output TSV path"),
  make_option("--lohhla_dir", type="character", default="LOHHLA",
              help="Path to LOHHLA repository root"),
  make_option("--purity", type="numeric", default=NA, help="Tumor purity (0-1)"),
  make_option("--ploidy", type="numeric", default=2, help="Tumor ploidy")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Read HLA alleles
hla_table <- fread(opt$hla_tsv)
hla_alleles <- hla_table$allele

# Source LOHHLA functions (adjust path to your LOHHLA install)
lohhla_script <- file.path(opt$lohhla_dir, "LOHHLA.R")
if (!file.exists(lohhla_script)) {
  stop("LOHHLA.R not found at: ", lohhla_script,
       "\n  Install: git clone https://github.com/ANHIG/LOHHLA.git")
}
source(lohhla_script)

# Build HLA FASTA from alleles using LOHHLA utility
hla_fasta <- paste0(opt$sample_id, "_hla.fa")
# LOHHLA uses HLA_seq_ref for each typed allele

run_LOHHLA(
  sample_id = opt$sample_id,
  patient = opt$sample_id,
  normal_bam = opt$normal_bam,
  tumor_bam = opt$tumor_bam,
  hla_alleles = hla_alleles,
  hla_ref_fasta = hla_fasta,
  hla_ref_dir = file.path(opt$lohhla_dir, "data", "ref"),
  output_dir = dirname(opt$output),
  purity = opt$purity,
  ploidy = opt$ploidy,
  copy_num = 2,
  max_cut = 0.5,
  min_cov = 20
)

# LOHHLA produces a result table with columns:
#   Allele, Copy_number, LOH_status, Pval, Telomeric, Centromeric
# Read and write in NeoFlux contract format
result_file <- file.path(dirname(opt$output), paste0(opt$sample_id, "_LOH_result.tsv"))
if (file.exists(result_file)) {
  result <- fread(result_file)
  output <- data.table(
    allele = result$Allele,
    is_loh = ifelse(result$LOH_status == "LOH", "True", "False"),
    loh_pvalue = as.character(result$Pval),
    copy_number = as.character(result$Copy_number)
  )
} else {
  # Fallback: all False
  output <- data.table(
    allele = hla_alleles,
    is_loh = "False",
    loh_pvalue = "",
    copy_number = ""
  )
  warning("LOHHLA result file not found; writing all-False placeholder")
}

fwrite(output, opt$output, sep="\t")
