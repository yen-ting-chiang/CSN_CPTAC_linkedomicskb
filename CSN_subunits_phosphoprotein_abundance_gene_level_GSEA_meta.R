## =========================================================================
##  CSN Subunits phosphoprotein abundance gene level GSEA Meta-Analysis
##  Data source: GSEA results from CSN_subunits_phosphoprotein_abundance_gene_level_GSEA
##
##  Purpose:
##    Perform cross-dataset meta-analysis using Stouffer's z-score method 
##    with Benjamini-Hochberg FDR correction.
##    This script reads the output from CSN_subunits_phosphoprotein_abundance_gene_level_GSEA.R
##    and combines p-values and effect directions (NES) across datasets
##    for each stratum, each predictor, and each gene set collection using Z-scores.
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
INPUT_DIR <- file.path(BASE_DIR, "CSN_subunits_phosphoprotein_abundance_gene_level_GSEA")
OUT_ROOT <- file.path(BASE_DIR, "CSN_subunits_phosphoprotein_abundance_gene_level_GSEA_meta")

# Minimum number of datasets a pathway must be present in to be included in meta-analysis
MIN_DATASETS <- 3

dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Find and process input files ------------------------------------

message("Searching for GSEA grand combined result file in: ", INPUT_DIR)
grand_file <- file.path(INPUT_DIR, "ALL_datasets_ALL_strata_GSEA_all_collections_all_CSN_subunits.csv")

if (file.exists(grand_file)) {
  message("Loading grand combined GSEA file...")
  full_dt <- data.table::fread(grand_file)
} else {
  message("Grand combined file not found. Scanning for per-dataset result files...")
  # Fallback to reading individual files if the grand file doesn't exist
  all_files <- list.files(INPUT_DIR, pattern = "_GSEA_.*_predictor_.*\\.csv$", recursive = TRUE, full.names = TRUE)
  
  # Exclude previously combined files to avoid duplication
  all_files <- all_files[!grepl("ALL_datasets", basename(all_files))]
  
  if (length(all_files) == 0) {
    stop("No GSEA result files found in the input directory.")
  }
  
  message(sprintf("Found %d result files. Loading all...", length(all_files)))
  dt_list <- lapply(all_files, data.table::fread)
  full_dt <- data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE)
}

# Ensure required columns are present
req_cols <- c("dataset", "stratum", "csn_subunit", "collection", "pathway", "NES", "pval", "size")
missing_cols <- setdiff(req_cols, names(full_dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns in input files: ", paste(missing_cols, collapse = ", "))
}

# Remove rows with NA p-values or NES
full_dt <- full_dt[!is.na(pval) & !is.na(NES)]

# Ensure p-values are not exactly 0 to avoid Inf Z-scores
min_p <- .Machine$double.xmin
full_dt[pval < min_p, pval := min_p]

# ---- 3. Meta-analysis function ------------------------------------------

run_stouffer_meta_gsea <- function(df) {
  # Stouffer's z-score method (Sample-size weighted Liptak-Stouffer method if n_samples exists, otherwise unweighted)
  
  # 1. Convert two-sided p-value to one-sided Z-score, keeping the sign of NES
  # Use lower.tail = FALSE to avoid precision loss for extremely small p-values
  z_scores <- sign(df$NES) * qnorm(df$pval / 2, lower.tail = FALSE)
  
  # 2. Weights
  if ("n_samples" %in% names(df) && all(!is.na(df$n_samples))) {
    w <- sqrt(df$n_samples)
  } else {
    w <- rep(1, nrow(df))
  }
  
  # 3. Calculate meta Z-score
  z_meta <- sum(w * z_scores) / sqrt(sum(w^2))
  
  # 4. Convert meta Z-score back to two-sided p-value
  p_meta <- 2 * pnorm(abs(z_meta), lower.tail = FALSE)
  
  # 5. Summarize other metrics
  mean_nes <- mean(df$NES)
  mean_size <- round(mean(df$size))
  n_datasets <- nrow(df)
  datasets_included <- paste(df$dataset, collapse = "|")
  
  # 6. Significant datasets calculation
  up_idx <- which(df$NES > 0 & df$pval < 0.05)
  down_idx <- which(df$NES < 0 & df$pval < 0.05)
  
  n_up_sig <- length(up_idx)
  up_sig_ds <- if (n_up_sig > 0) paste(df$dataset[up_idx], collapse = "|") else ""
  
  n_down_sig <- length(down_idx)
  down_sig_ds <- if (n_down_sig > 0) paste(df$dataset[down_idx], collapse = "|") else ""
  
  list(
    meta_Z_score = z_meta,
    meta_p_value = p_meta,
    mean_NES = mean_nes,
    mean_size = mean_size,
    n_datasets = n_datasets,
    datasets = datasets_included,
    n_up_sig_datasets = n_up_sig,
    up_sig_datasets = up_sig_ds,
    n_down_sig_datasets = n_down_sig,
    down_sig_datasets = down_sig_ds
  )
}

# ---- 4. Execute meta-analysis per stratum, per predictor, per collection --

strata <- unique(full_dt$stratum)

for (strat in strata) {
  strat_dt <- full_dt[stratum == strat]
  predictors <- unique(strat_dt$csn_subunit)
  collections <- unique(strat_dt$collection)
  
  strat_out_dir <- file.path(OUT_ROOT, strat)
  dir.create(strat_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pred in predictors) {
    for (coll in collections) {
      message(sprintf("Processing Stratum: %s | Predictor: %s | Collection: %s", strat, pred, coll))
      
      pred_coll_dt <- strat_dt[csn_subunit == pred & collection == coll]
      if (nrow(pred_coll_dt) == 0) next
      
      # Group by pathway and apply meta-analysis
      meta_res <- pred_coll_dt[, run_stouffer_meta_gsea(.SD), by = .(pathway)]
      
      # Filter by minimum datasets
      meta_res <- meta_res[n_datasets >= MIN_DATASETS]
      
      if (nrow(meta_res) == 0) {
        message(sprintf("  No pathways passed the MIN_DATASETS >= %d filter. Skipping.", MIN_DATASETS))
        next
      }
      
      # Get wide format of NES and pval per dataset
      wide_dt <- data.table::dcast(pred_coll_dt[pathway %in% meta_res$pathway], 
                                   pathway ~ dataset, 
                                   value.var = c("NES", "pval"))
      
      # Merge wide data into meta_res
      meta_res <- merge(meta_res, wide_dt, by = "pathway", all.x = TRUE)
      
      # Apply Benjamini-Hochberg FDR correction
      meta_res[, meta_BH_FDR := p.adjust(meta_p_value, method = "BH")]
      
      # Add metadata columns back for clarity
      meta_res[, csn_subunit := pred]
      meta_res[, collection := coll]
      
      # Reorder columns for readability
      fixed_cols <- c("pathway", 
                      "n_up_sig_datasets", "up_sig_datasets",
                      "n_down_sig_datasets", "down_sig_datasets",
                      "collection", "csn_subunit",
                      "meta_p_value", "meta_BH_FDR", "meta_Z_score", 
                      "mean_NES", "mean_size", "n_datasets", "datasets")
      
      other_cols <- setdiff(names(meta_res), fixed_cols)
      # Sort other_cols so NES and pval columns are ordered alphabetically by dataset
      other_cols <- sort(other_cols)
      
      setcolorder(meta_res, c(fixed_cols, other_cols))
      
      # Sort by meta Z-score descending
      meta_res <- meta_res[order(-meta_Z_score)]
      
      # Format output file name
      out_file <- file.path(strat_out_dir, paste0(strat, "_meta_GSEA_", coll, "_predictor_", pred, ".csv"))
      
      # Save to CSV
      data.table::fwrite(meta_res, out_file)
      
      n_sig_005 <- sum(meta_res$meta_BH_FDR < 0.05, na.rm = TRUE)
      n_sig_025 <- sum(meta_res$meta_BH_FDR < 0.25, na.rm = TRUE)
      message(sprintf("  Saved %s (Pathways: %d | FDR < 0.05: %d | FDR < 0.25: %d)", 
                      basename(out_file), nrow(meta_res), n_sig_005, n_sig_025))
    }
  }
}

message("\nAll GSEA meta-analyses completed successfully! Results saved in: ", OUT_ROOT)
