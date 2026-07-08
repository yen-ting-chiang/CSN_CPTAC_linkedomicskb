## =========================================================================
##  CSN Subunits Protein Abundance vs Cell Type Correlation Meta-Analysis
##  Data source: Correlation results from
##               CSN_subunits_protein_abundance_cell_type_correlation.R
##
##  Purpose:
##    Perform cross-dataset meta-analysis of Pearson partial correlation
##    results (CSN subunit protein abundance vs cell type deconvolution
##    scores) using Stouffer's z-score method with sample-size weighting,
##    followed by Benjamini-Hochberg FDR correction.
##
##  Methodology:
##    Stouffer's z-score method (Liptak-Stouffer):
##      1. Convert each per-dataset two-sided p-value to a one-sided
##         z-score, preserving the sign of the Pearson r.
##      2. Weight z-scores by sqrt(n_samples) per dataset.
##      3. Combine: Z_meta = sum(w_i * z_i) / sqrt(sum(w_i^2))
##      4. Convert Z_meta back to a two-sided p-value.
##      5. Apply Benjamini-Hochberg FDR correction.
##
##    Additionally, the script reports the sample-size-weighted mean
##    Pearson r across datasets as an effect-size summary.
##
##  Methodology Reference:
##    - Stouffer, S. A. et al. (1949). The American Soldier: Adjustment
##      During Army Life. Princeton University Press.
##    - Willer, C. J., Li, Y., & Abecasis, G. R. (2010). METAL: fast
##      and efficient meta-analysis of genomewide association scans.
##      Bioinformatics, 26(17), 2190-2191.
##      https://doi.org/10.1093/bioinformatics/btq340
##    - Benjamini, Y., & Hochberg, Y. (1995). Controlling the false
##      discovery rate: a practical and powerful approach to multiple
##      testing. J. R. Stat. Soc. Ser. B, 57(1), 289-300.
##
##  Output structure:
##    CSN_subunits_protein_abundance_cell_type_correlation_meta/
##      {TP53_all, TP53_MUT, TP53_WT}/
##        {CIBERSORT, ESTIMATE, xCell}/
##          {stratum}_meta_{prefix}_correlation_predictor_{pred}.csv
##        {stratum}_meta_ALL_prefixes_correlation_all_predictors.csv
##      ALL_strata_ALL_prefixes_correlation_meta_all_predictors.csv
##      ALL_strata_ALL_prefixes_correlation_meta_summary.csv
## =========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ---- 1. Global configuration --------------------------------------------

BASE_DIR  <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb"
INPUT_DIR <- file.path(BASE_DIR, "CSN_subunits_protein_abundance_cell_type_correlation")
OUT_ROOT  <- file.path(BASE_DIR, "CSN_subunits_protein_abundance_cell_type_correlation_meta")

# Minimum number of datasets a cell type must be present in for meta-analysis
MIN_DATASETS <- 3

dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Find and load input files ----------------------------------------

message("Scanning for correlation result files in: ", INPUT_DIR)
all_files <- list.files(
  INPUT_DIR,
  pattern    = "_correlation_predictor_.*\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)

if (length(all_files) == 0) {
  stop("No correlation result files found in the input directory.")
}

message(sprintf("Found %d result files.", length(all_files)))

# Read all files into a single data.table
message("Loading all result files...")
dt_list <- lapply(all_files, data.table::fread)
full_dt <- data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# Ensure required columns are present
req_cols <- c("dataset", "stratum", "prefix", "predictor",
              "cell_type", "cell_type_col",
              "pearson_r", "t_statistic", "p_value", "n_samples")
missing_cols <- setdiff(req_cols, names(full_dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns in input files: ",
       paste(missing_cols, collapse = ", "))
}

# Remove rows with NA p-values or Pearson r
full_dt <- full_dt[!is.na(p_value) & !is.na(pearson_r)]

# Ensure p-values are not exactly 0 to avoid Inf Z-scores
min_p <- .Machine$double.xmin
full_dt[p_value < min_p, p_value := min_p]

message(sprintf("Total records loaded: %d", nrow(full_dt)))
message(sprintf("Strata found: %s", paste(unique(full_dt$stratum), collapse = ", ")))
message(sprintf("Prefixes found: %s", paste(unique(full_dt$prefix), collapse = ", ")))
message(sprintf("Predictors found: %s", paste(unique(full_dt$predictor), collapse = ", ")))

# ---- 3. Meta-analysis function: Stouffer's z-score method ----------------

run_stouffer_meta_corr <- function(df) {
  # Stouffer's z-score method (sample-size weighted)

  # 1. Convert two-sided p-value to one-sided Z-score,
  #    preserving the sign of the Pearson r
  z_scores <- sign(df$pearson_r) * qnorm(df$p_value / 2, lower.tail = FALSE)

  # 2. Weights based on square root of sample size
  if (all(!is.na(df$n_samples))) {
    w <- sqrt(df$n_samples)
  } else {
    w <- rep(1, nrow(df))
  }

  # 3. Calculate meta Z-score
  z_meta <- sum(w * z_scores) / sqrt(sum(w^2))

  # 4. Convert meta Z-score back to two-sided p-value
  p_meta <- 2 * pnorm(abs(z_meta), lower.tail = FALSE)

  # 5. Sample-size-weighted mean Pearson r as effect-size summary
  w_n <- df$n_samples
  w_n[is.na(w_n)] <- 1
  weighted_mean_r <- sum(w_n * df$pearson_r) / sum(w_n)

  # 6. Summary metrics
  n_datasets     <- nrow(df)
  total_samples  <- sum(df$n_samples, na.rm = TRUE)
  datasets_incl  <- paste(df$dataset, collapse = "|")
  n_positive     <- sum(df$pearson_r > 0, na.rm = TRUE)
  n_negative     <- sum(df$pearson_r < 0, na.rm = TRUE)
  direction_consistency <- max(n_positive, n_negative) / n_datasets

  list(
    meta_Z_score             = z_meta,
    meta_p_value             = p_meta,
    weighted_mean_pearson_r  = weighted_mean_r,
    mean_pearson_r           = mean(df$pearson_r),
    n_datasets               = n_datasets,
    total_samples            = total_samples,
    n_positive               = n_positive,
    n_negative               = n_negative,
    direction_consistency    = direction_consistency,
    datasets                 = datasets_incl
  )
}

# ---- 4. Execute meta-analysis per stratum, per prefix, per predictor -----

strata   <- unique(full_dt$stratum)
prefixes <- unique(full_dt$prefix)

all_meta_results <- list()

for (strat in strata) {
  for (pfx in prefixes) {

    sub_dt <- full_dt[stratum == strat & prefix == pfx]
    if (nrow(sub_dt) == 0) next

    predictors <- unique(sub_dt$predictor)

    pfx_out_dir <- file.path(OUT_ROOT, strat, pfx)
    dir.create(pfx_out_dir, recursive = TRUE, showWarnings = FALSE)

    for (pred in predictors) {
      message(sprintf("Processing Stratum: %s | Prefix: %s | Predictor: %s",
                      strat, pfx, pred))

      pred_dt <- sub_dt[predictor == pred]

      # Group by cell_type_col (the full column name, which is unique per
      # cell type even across prefixes) and apply meta-analysis
      meta_res <- pred_dt[,
        run_stouffer_meta_corr(.SD),
        by = .(cell_type, cell_type_col)
      ]

      # Filter by minimum number of datasets
      meta_res <- meta_res[n_datasets >= MIN_DATASETS]

      if (nrow(meta_res) == 0) {
        message(sprintf("  No cell types passed MIN_DATASETS >= %d filter. Skipping.",
                        MIN_DATASETS))
        next
      }

      # Apply Benjamini-Hochberg FDR correction
      meta_res[, meta_BH_FDR := p.adjust(meta_p_value, method = "BH")]

      # Add metadata columns
      meta_res[, `:=`(stratum = strat, prefix = pfx, predictor = pred)]

      # Reorder columns for readability
      setcolorder(meta_res, c(
        "stratum", "prefix", "predictor",
        "cell_type", "cell_type_col",
        "meta_p_value", "meta_BH_FDR", "meta_Z_score",
        "weighted_mean_pearson_r", "mean_pearson_r",
        "n_datasets", "total_samples",
        "n_positive", "n_negative", "direction_consistency",
        "datasets"
      ))

      # Sort by meta p-value
      meta_res <- meta_res[order(meta_p_value)]

      # Save per-predictor CSV
      out_file <- file.path(pfx_out_dir,
        paste0(strat, "_meta_", pfx, "_correlation_predictor_", pred, ".csv"))
      data.table::fwrite(meta_res, out_file)

      n_sig_005 <- sum(meta_res$meta_BH_FDR < 0.05, na.rm = TRUE)
      n_sig_010 <- sum(meta_res$meta_BH_FDR < 0.10, na.rm = TRUE)
      message(sprintf("  Saved %s (Cell types: %d | FDR < 0.05: %d | FDR < 0.10: %d)",
                      basename(out_file), nrow(meta_res), n_sig_005, n_sig_010))

      # Accumulate results
      all_meta_results[[paste(strat, pfx, pred, sep = "|")]] <- meta_res
    }
  }
}

# ---- 5. Grand combined output -------------------------------------------

if (length(all_meta_results) > 0) {

  grand_combined <- data.table::rbindlist(all_meta_results, use.names = TRUE, fill = TRUE)

  # --- Per-stratum x per-prefix combined files ---
  for (strat in strata) {
    for (pfx in prefixes) {
      sub_data <- grand_combined[stratum == strat & prefix == pfx]
      if (nrow(sub_data) == 0) next

      sub_dir <- file.path(OUT_ROOT, strat, pfx)
      dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)

      combined_csv <- file.path(sub_dir,
        paste0(strat, "_meta_", pfx, "_correlation_all_predictors.csv"))
      data.table::fwrite(sub_data, combined_csv)
      message(sprintf("\n[%s|%s] Combined CSV saved: %s", strat, pfx, basename(combined_csv)))
    }
  }

  # --- Per-stratum combined files (all prefixes) ---
  for (strat in strata) {
    st_data <- grand_combined[stratum == strat]
    if (nrow(st_data) == 0) next

    st_dir <- file.path(OUT_ROOT, strat)
    dir.create(st_dir, showWarnings = FALSE, recursive = TRUE)

    combined_csv <- file.path(st_dir,
      paste0(strat, "_meta_ALL_prefixes_correlation_all_predictors.csv"))
    data.table::fwrite(st_data, combined_csv)
    message(sprintf("\n[%s] Combined CSV (all prefixes) saved: %s", strat, basename(combined_csv)))
  }

  # --- Grand combined file ---
  grand_file <- file.path(OUT_ROOT,
    "ALL_strata_ALL_prefixes_correlation_meta_all_predictors.csv")
  data.table::fwrite(grand_combined, grand_file)
  message(sprintf("\nGrand combined CSV saved: %s", grand_file))

  # --- Summary table ---
  summary_df <- grand_combined[, .(
    n_cell_types_tested = .N,
    n_sig_FDR_0.05      = sum(meta_BH_FDR < 0.05, na.rm = TRUE),
    n_sig_FDR_0.10      = sum(meta_BH_FDR < 0.10, na.rm = TRUE),
    n_sig_FDR_0.25      = sum(meta_BH_FDR < 0.25, na.rm = TRUE),
    n_pos_sig_FDR_0.05  = sum(meta_BH_FDR < 0.05 & weighted_mean_pearson_r > 0, na.rm = TRUE),
    n_neg_sig_FDR_0.05  = sum(meta_BH_FDR < 0.05 & weighted_mean_pearson_r < 0, na.rm = TRUE),
    mean_n_datasets     = round(mean(n_datasets), 1),
    mean_total_samples  = round(mean(total_samples), 0)
  ), by = .(stratum, prefix, predictor)]

  summary_file <- file.path(OUT_ROOT,
    "ALL_strata_ALL_prefixes_correlation_meta_summary.csv")
  data.table::fwrite(summary_df, summary_file)
  message(sprintf("Summary CSV saved: %s", summary_file))

  # Per-stratum summary files
  for (strat in strata) {
    st_summary <- summary_df[stratum == strat]
    if (nrow(st_summary) == 0) next
    st_dir <- file.path(OUT_ROOT, strat)
    st_summary_file <- file.path(st_dir,
      paste0(strat, "_meta_correlation_summary.csv"))
    data.table::fwrite(st_summary, st_summary_file)
    message(sprintf("[%s] Summary CSV saved: %s", strat, basename(st_summary_file)))
  }

  # Print summary to console
  message("\n========== Correlation Meta-Analysis Summary ==========")
  print(as.data.frame(summary_df), row.names = FALSE)

} else {
  message("\nNo meta-analysis results were generated.")
}


message("\n============================================================")
message("  All meta-analyses completed!")
message(sprintf("  Output directory: %s", OUT_ROOT))
message("============================================================\n")
