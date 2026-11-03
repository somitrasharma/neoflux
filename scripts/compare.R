# compare.R
# Compare top-10 neo-epitopes from NeoFlux and pVACseq per sample

library(dplyr)
library(readr)
library(stringr)

# Robust helper to read top-10 per sample by variant binding affinity
read_top10 <- function(file, id_col, peptide_col, ic50_col) {
  if (!file.exists(file)) {
    warning(paste("File missing:", file))
    return(NULL)
  }
  
  df <- read_tsv(file, col_types = cols(), lazy = FALSE)
  # Normalize column names to lowercase to prevent casing mismatches
  names(df) <- tolower(names(df))
  
  id_col_lower <- tolower(id_col)
  peptide_col_lower <- tolower(peptide_col)
  ic50_col_lower <- tolower(ic50_col)
  
  if (!ic50_col_lower %in% names(df) || !peptide_col_lower %in% names(df)) {
    warning(paste("Required columns missing in:", file))
    return(NULL)
  }
  
  df %>%
    mutate(parsed_ic50 = as.numeric(!!sym(ic50_col_lower))) %>%
    filter(!is.na(parsed_ic50)) %>%
    arrange(parsed_ic50) %>%
    slice_head(n = 10) %>%
    select(
      sample_id = any_of(id_col_lower),
      peptide = !!sym(peptide_col_lower),
      ic50 = parsed_ic50
    )
}

# Strict pathing relative to execution context
neo_dir <- "demo_results"
pvac_dir <- "pvacseq_results"

# 1. Collect NeoFlux top-10 mutations
neo_files <- list.files(neo_dir, pattern = "_ranked\\.tsv$", full.names = TRUE, recursive = TRUE)
neo_list <- lapply(neo_files, function(f) {
  # Robustly extract whatever sits before '_ranked.tsv' as the sample name
  sample_name <- str_remove(basename(f), "_ranked\\.tsv$")
  df <- read_top10(f, id_col = "sample_id", peptide_col = "mutant_peptide", ic50_col = "best_ic50")
  if (!is.null(df)) df$sample_id <- sample_name
  return(df)
})
neo_top10 <- bind_rows(neo_list)

# 2. Collect pVACseq top-10 mutations (looking for subfolder filtered files)
pvac_files <- list.files(pvac_dir, pattern = "filtered\\.tsv$", full.names = TRUE, recursive = TRUE)
pvac_list <- lapply(pvac_files, function(f) {
  # Deduce sample name from the parent directory name if the file is generic
  sample_name <- basename(dirname(f))
  df <- read_top10(f, id_col = "sample", peptide_col = "peptide", ic50_col = "ic50")
  if (!is.null(df)) df$sample_id <- sample_name
  return(df)
})
pvac_top10 <- bind_rows(pvac_list)

# 3. Validation checkpoint
if (nrow(neo_top10) == 0 || nrow(pvac_top10) == 0) {
  stop("CRITICAL ERROR: Zero rows loaded. Verify file names and locations in demo_results/ and pvacseq_results/.")
}

# 4. Process Intersection Analysis per Sample Cohort
samples <- unique(c(neo_top10$sample_id, pvac_top10$sample_id))
metrics_list <- list()

for (s in samples) {
  neo_peps <- neo_top10 %>% filter(sample_id == s) %>% pull(peptide)
  pvac_peps <- pvac_top10 %>% filter(sample_id == s) %>% pull(peptide)
  
  intersect_count <- length(intersect(neo_peps, pvac_peps))
  
  metrics_list[[s]] <- tibble(
    sample_id = s,
    neoflux_top10_count = length(neo_peps),
    pvacseq_top10_count = length(pvac_peps),
    shared_peptides = intersect_count,
    precision_at_10 = ifelse(length(pvac_peps) > 0, intersect_count / length(pvac_peps), 0),
    recall_at_10 = ifelse(length(neo_peps) > 0, intersect_count / length(neo_peps), 0)
  )
}

overlap_metrics <- bind_rows(metrics_list)
write_csv(overlap_metrics, "results/comparison_metrics.tsv")
cat("\n================ VALIDATION STATISTICS SUCCESS ================\n")
print(overlap_metrics)
cat("===============================================================\n")
