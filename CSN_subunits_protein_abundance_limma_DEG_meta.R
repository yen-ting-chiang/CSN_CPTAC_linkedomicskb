## =========================================================================
##  CSN Subunits Protein Abundance Limma DEG Meta-Analysis
##  Data source: Limma DEG results from CPTAC datasets
##
##  Purpose:
##    Perform cross-dataset meta-analysis using Stouffer's z-score method 
##    with Benjamini-Hochberg FDR correction.
##    This script reads the output from CSN_subunits_protein_abundance_limma_DEG.R
##    and combines p-values and effect directions (logFC) across datasets
##    for each stratum and each predictor using sample-size weighted Z-scores.
##
##  Methodology Reference:
##    - Willer, C. J., Li, Y., & Abecasis, G. R. (2010). METAL: fast and efficient 
##      meta-analysis of genomewide association scans. Bioinformatics, 26(17), 2190-2191.
##      https://doi.org/10.1093/bioinformatics/btq340
##    - Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: 
##      a practical and powerful approach to multiple testing. J. R. Stat. Soc. Ser. B, 289-300.
## =========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ---- 1. Global configuration --------------------------------------------

BASE_DIR <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb"
INPUT_DIR <- file.path(BASE_DIR, "CSN_subunits_protein_abundance_limma_DEG")
OUT_ROOT <- file.path(BASE_DIR, "CSN_subunits_protein_abundance_limma_DEG_meta")

# Minimum number of datasets a gene must be present in to be included in meta-analysis
MIN_DATASETS <- 3

dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Find and process input files ------------------------------------

message("Scanning for Limma DEG result files in: ", INPUT_DIR)
all_files <- list.files(INPUT_DIR, pattern = "_limma_DEG_predictor_.*\\.csv$", recursive = TRUE, full.names = TRUE)

if (length(all_files) == 0) {
  stop("No Limma DEG result files found in the input directory.")
}

message(sprintf("Found %d result files.", length(all_files)))

# Read all files into a single data.table
message("Loading all result files...")
dt_list <- lapply(all_files, data.table::fread)
full_dt <- data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# Ensure required columns are present
req_cols <- c("dataset", "stratum", "csn_subunit", "ensembl_id", "ensembl_id_base", 
              "gene_symbol", "logFC", "p_value", "n_samples", "is_predictor")
missing_cols <- setdiff(req_cols, names(full_dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns in input files: ", paste(missing_cols, collapse = ", "))
}

# Remove rows with NA p-values or logFCs
full_dt <- full_dt[!is.na(p_value) & !is.na(logFC)]

# Ensure p-values are not exactly 0 to avoid Inf Z-scores
min_p <- .Machine$double.xmin
full_dt[p_value < min_p, p_value := min_p]

# ---- 3. Meta-analysis function ------------------------------------------

run_stouffer_meta <- function(df) {
  # Stouffer's z-score method (Sample-size weighted Liptak-Stouffer method)
  
  # 1. Convert two-sided p-value to one-sided Z-score, keeping the sign of logFC
  # Use lower.tail = FALSE to avoid precision loss for extremely small p-values
  z_scores <- sign(df$logFC) * qnorm(df$p_value / 2, lower.tail = FALSE)
  
  # 2. Weights based on square root of sample size (standard approach in meta-analysis)
  if (all(!is.na(df$n_samples))) {
    w <- sqrt(df$n_samples)
  } else {
    w <- rep(1, nrow(df))
  }
  
  # 3. Calculate meta Z-score
  z_meta <- sum(w * z_scores) / sqrt(sum(w^2))
  
  # 4. Convert meta Z-score back to two-sided p-value
  p_meta <- 2 * pnorm(abs(z_meta), lower.tail = FALSE)
  
  # 5. Summarize other metrics
  mean_logfc <- mean(df$logFC)
  n_datasets <- nrow(df)
  total_samples <- sum(df$n_samples, na.rm = TRUE)
  datasets_included <- paste(df$dataset, collapse = "|")
  is_pred <- any(df$is_predictor, na.rm = TRUE)
  
  list(
    meta_Z_score = z_meta,
    meta_p_value = p_meta,
    mean_logFC = mean_logfc,
    n_datasets = n_datasets,
    total_samples = total_samples,
    datasets = datasets_included,
    is_predictor = is_pred
  )
}

# ---- 4. Execute meta-analysis per stratum, per predictor, per gene -------

strata <- unique(full_dt$stratum)

for (strat in strata) {
  strat_dt <- full_dt[stratum == strat]
  predictors <- unique(strat_dt$csn_subunit)
  
  strat_out_dir <- file.path(OUT_ROOT, strat)
  dir.create(strat_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pred in predictors) {
    message(sprintf("Processing Stratum: %s | Predictor: %s", strat, pred))
    
    pred_dt <- strat_dt[csn_subunit == pred]
    
    # Group by gene and apply meta-analysis
    meta_res <- pred_dt[, run_stouffer_meta(.SD), by = .(ensembl_id, ensembl_id_base, gene_symbol)]
    
    # Filter by minimum datasets
    meta_res <- meta_res[n_datasets >= MIN_DATASETS]
    
    if (nrow(meta_res) == 0) {
      message(sprintf("  No genes passed the MIN_DATASETS >= %d filter. Skipping.", MIN_DATASETS))
      next
    }
    
    # Apply Benjamini-Hochberg FDR correction
    meta_res[, meta_BH_FDR := p.adjust(meta_p_value, method = "BH")]
    
    # Reorder columns for readability
    setcolorder(meta_res, c("ensembl_id", "ensembl_id_base", "gene_symbol", 
                            "meta_p_value", "meta_BH_FDR", "meta_Z_score", 
                            "mean_logFC", "n_datasets", "total_samples", 
                            "is_predictor", "datasets"))
    
    # Sort by meta p-value
    meta_res <- meta_res[order(meta_p_value)]
    
    # Format output file name
    out_file <- file.path(strat_out_dir, paste0(strat, "_meta_limma_DEG_predictor_", pred, ".csv"))
    
    # Save to CSV
    data.table::fwrite(meta_res, out_file)
    
    n_sig_005 <- sum(meta_res$meta_BH_FDR < 0.05, na.rm = TRUE)
    message(sprintf("  Saved %s (Genes: %d | FDR < 0.05: %d)", basename(out_file), nrow(meta_res), n_sig_005))
  }
}

message("\nAll meta-analyses completed successfully! Results saved in: ", OUT_ROOT)
