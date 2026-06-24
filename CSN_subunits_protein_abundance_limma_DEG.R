## =========================================================================
##  CSN Subunits Protein Abundance Limma DEG Analysis
##  Data source: LinkedOmicsKB CPTAC datasets
##
##  Purpose:
##    For each CPTAC dataset and each CSN subunit, use the subunit's protein
##    abundance as a continuous predictor in limma to perform genome-wide
##    Differentially Expressed Protein (DEP) analysis across the proteome.
##
##  Statistical rationale:
##    Using a continuous predictor in limma is a well-established approach
##    (limma User's Guide Section 9.4). It identifies genes whose protein
##    abundance changes linearly with the CSN subunit level, after adjusting
##    for covariates. limma's empirical Bayes moderation stabilizes variance
##    estimates and yields well-calibrated p-values.
##    Reference: Ritchie ME et al. (2015) Nucleic Acids Res 43(7):e47.
##               https://doi.org/10.1093/nar/gkv007
##
##  Covariates adjusted: Sex, Age, Tumor Purity (WES_purity)
##  Multiple testing correction: Benjamini-Hochberg FDR
##  Proteomics data: NOT imputed (genes with excessive missing values are
##                   filtered out; limma handles remaining NAs internally)
##  Covariate data: Imputed per covariate using the most appropriate method:
##                  - Sex (categorical): mode imputation
##                  - Age (continuous):  median imputation
##                  - Tumor purity (continuous): median imputation
##  Stratification: TP53 mutation status (TP53_all, TP53_MUT, TP53_WT)
## =========================================================================

# ---- 0. Load / install required packages ---------------------------------

required_cran <- c("data.table", "dplyr", "readr", "stringr")
required_bioc <- c("limma", "AnnotationDbi", "org.Hs.eg.db")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ---- 1. Global configuration --------------------------------------------

BASE_DIR <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb"
OUT_ROOT <- file.path(BASE_DIR, "CSN_subunits_protein_abundance_limma_DEG")
TP53_DIR <- file.path(BASE_DIR, "TP53_mutation_classification")
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

DATASETS <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC",
               "LSCC", "LUAD", "OV", "PDAC", "UCEC")

CSN_SUBUNITS <- c("GPS1", "COPS2", "COPS3", "COPS4", "COPS5",
                   "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9")

# TP53 stratification strata
TP53_STRATA <- c("TP53_all", "TP53_MUT", "TP53_WT")

# Datasets excluded from TP53 stratification:
#   CCRCC: too few TP53 mutant samples (5/103)
#   LSCC:  too few TP53 wild-type samples (5/108)
#   OV:    too few TP53 wild-type samples (5/82)
DS_SKIP_TP53_STRATIFICATION <- c("CCRCC", "LSCC", "OV")

# Minimum number of non-NA samples per gene to include in limma
MIN_SAMPLES_PER_GENE <- 10L

# Minimum proportion of non-NA values per gene (relative to valid samples)
MIN_NONMISSING_PROP <- 0.3

# ---- 2. Build Ensembl-to-gene-symbol mapping via org.Hs.eg.db -----------

message("[Mapping] Building comprehensive Ensembl gene ID to symbol mapping ...")

## Map ALL Ensembl IDs to gene symbols for genome-wide annotation
all_ensembl_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = keys(org.Hs.eg.db, keytype = "ENSEMBL"),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
)
all_ensembl_map <- all_ensembl_map[!is.na(all_ensembl_map$SYMBOL), ]
# Keep first mapping for duplicates
all_ensembl_map <- all_ensembl_map[!duplicated(all_ensembl_map$ENSEMBL), ]

## Named vector: ENSEMBL (no version) -> SYMBOL
ensembl_to_symbol <- setNames(all_ensembl_map$SYMBOL, all_ensembl_map$ENSEMBL)

message(sprintf("[Mapping] Built mapping for %d Ensembl IDs to gene symbols",
                length(ensembl_to_symbol)))

## Build CSN-specific mapping (SYMBOL -> ENSEMBL)
csn_ensembl_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = CSN_SUBUNITS,
  columns = c("SYMBOL", "ENSEMBL"),
  keytype = "SYMBOL"
)
csn_ensembl_map <- csn_ensembl_map[!is.na(csn_ensembl_map$ENSEMBL), ]
csn_ensembl_map <- csn_ensembl_map[!duplicated(csn_ensembl_map$SYMBOL), ]
symbol_to_ensembl <- setNames(csn_ensembl_map$ENSEMBL, csn_ensembl_map$SYMBOL)

message(sprintf("[Mapping] Mapped %d / %d CSN subunits to Ensembl IDs: %s",
                nrow(csn_ensembl_map), length(CSN_SUBUNITS),
                paste(csn_ensembl_map$SYMBOL, collapse = ", ")))

# ---- 3. Helper: impute covariates (per covariate best-practice) ----------

impute_covariates <- function(cov_df, ds_id, stratum = "TP53_all") {
  ## cov_df: data.frame with columns sex, age, purity
  ## Returns: list(imputed_df, imputation_log)
  log_entries <- list()

  # --- Sex: mode imputation (categorical) ---
  if ("sex" %in% names(cov_df)) {
    n_miss <- sum(is.na(cov_df$sex))
    n_total <- nrow(cov_df)
    if (n_miss > 0 && n_miss < n_total) {
      mode_val <- names(sort(table(cov_df$sex), decreasing = TRUE))[1]
      cov_df$sex[is.na(cov_df$sex)] <- mode_val
      log_entries[["sex"]] <- sprintf(
        "Imputed %d / %d missing Sex values with mode = '%s'",
        n_miss, n_total, mode_val
      )
    } else if (n_miss == n_total) {
      log_entries[["sex"]] <- "All Sex values missing; column will be dropped"
    } else {
      log_entries[["sex"]] <- sprintf(
        "No missing Sex values (%d samples)", n_total
      )
    }
  }

  # --- Age: median imputation (continuous) ---
  if ("age" %in% names(cov_df)) {
    n_miss <- sum(is.na(cov_df$age))
    n_total <- nrow(cov_df)
    if (n_miss > 0 && n_miss < n_total) {
      med_val <- median(cov_df$age, na.rm = TRUE)
      cov_df$age[is.na(cov_df$age)] <- med_val
      log_entries[["age"]] <- sprintf(
        "Imputed %d / %d missing Age values with median = %.1f",
        n_miss, n_total, med_val
      )
    } else if (n_miss == n_total) {
      log_entries[["age"]] <- "All Age values missing; column will be dropped"
    } else {
      log_entries[["age"]] <- sprintf(
        "No missing Age values (%d samples)", n_total
      )
    }
  }

  # --- Tumor purity: median imputation (continuous) ---
  if ("purity" %in% names(cov_df)) {
    n_miss <- sum(is.na(cov_df$purity))
    n_total <- nrow(cov_df)
    if (n_miss > 0 && n_miss < n_total) {
      med_val <- median(cov_df$purity, na.rm = TRUE)
      cov_df$purity[is.na(cov_df$purity)] <- med_val
      log_entries[["purity"]] <- sprintf(
        "Imputed %d / %d missing WES_purity values with median = %.4f",
        n_miss, n_total, med_val
      )
    } else if (n_miss == n_total) {
      log_entries[["purity"]] <- "All WES_purity values missing; column will be dropped"
    } else {
      log_entries[["purity"]] <- sprintf(
        "No missing WES_purity values (%d samples)", n_total
      )
    }
  }

  list(imputed_df = cov_df, imputation_log = log_entries)
}

# ---- 4. Helper: read TP53 classification and split sample IDs ------------

get_tp53_sample_lists <- function(ds_id, tp53_dir = TP53_DIR) {
  ## Returns a named list:
  ##   TP53_all = all sample IDs
  ##   TP53_MUT = TP53 mutant sample IDs (mt == 1)
  ##   TP53_WT  = TP53 wild-type sample IDs (wt == 1)
  tp53_file <- file.path(tp53_dir, paste0(ds_id, "_TP53_classification.csv"))
  if (!file.exists(tp53_file)) {
    message(sprintf("[TP53] %s: classification file not found: %s",
                    ds_id, basename(tp53_file)))
    return(NULL)
  }

  tp53_df <- data.table::fread(tp53_file, header = TRUE, check.names = FALSE)

  all_ids <- as.character(tp53_df$sample_id)
  mut_ids <- as.character(tp53_df$sample_id[tp53_df$mt == 1])
  wt_ids  <- as.character(tp53_df$sample_id[tp53_df$wt == 1])

  message(sprintf("[TP53] %s: Total=%d, TP53mut=%d, TP53wt=%d",
                  ds_id, length(all_ids), length(mut_ids), length(wt_ids)))

  list(
    TP53_all = all_ids,
    TP53_MUT = mut_ids,
    TP53_WT  = wt_ids
  )
}

# ---- 5. Main function: limma DEG for one dataset + stratum ---------------

run_limma_deg_one_dataset <- function(
    ds_id,
    stratum       = "TP53_all",
    sample_subset = NULL,
    base_dir      = BASE_DIR,
    out_root      = OUT_ROOT,
    subunits      = CSN_SUBUNITS,
    ens2sym       = ensembl_to_symbol,
    sym2ens       = symbol_to_ensembl,
    min_n_gene    = MIN_SAMPLES_PER_GENE,
    min_prop      = MIN_NONMISSING_PROP
) {

  ds_dir <- file.path(base_dir, ds_id)
  ds_out <- file.path(out_root, stratum, ds_id)
  dir.create(ds_out, recursive = TRUE, showWarnings = FALSE)

  message(sprintf("\n========== Processing dataset: %s | Stratum: %s ==========",
                  ds_id, stratum))

  # --- 5a. Read proteomics matrix ------------------------------------------

  prot_file <- file.path(ds_dir,
    paste0(ds_id,
           "_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"))

  if (!file.exists(prot_file)) {
    message(sprintf("[%s|%s] Proteomics file not found: %s",
                    ds_id, stratum, basename(prot_file)))
    return(invisible(NULL))
  }

  prot_dt <- data.table::fread(prot_file, header = TRUE, check.names = FALSE)
  idx_col <- names(prot_dt)[1]  # "idx"
  ensembl_ids_raw <- prot_dt[[idx_col]]

  # Strip version suffix from Ensembl IDs (e.g., ENSG00000008083.17)
  ensembl_ids_base <- sub("\\.[0-9]+$", "", ensembl_ids_raw)

  # Map ALL Ensembl IDs to gene symbols
  gene_symbols <- ens2sym[ensembl_ids_base]
  # For unmapped Ensembl IDs, keep the raw Ensembl ID as identifier
  gene_symbols[is.na(gene_symbols)] <- ensembl_ids_raw[is.na(gene_symbols)]
  names(gene_symbols) <- ensembl_ids_raw

  # Sample IDs (all columns except the index column)
  all_sample_ids <- setdiff(names(prot_dt), idx_col)

  # Build numeric matrix (genes x samples)
  expr_mat <- as.matrix(prot_dt[, ..all_sample_ids])
  rownames(expr_mat) <- ensembl_ids_raw
  storage.mode(expr_mat) <- "double"

  # --- 5a2. Subset samples by TP53 stratum ---------------------------------

  if (!is.null(sample_subset)) {
    sam_stratum <- intersect(sample_subset, all_sample_ids)
    if (length(sam_stratum) < 10) {
      message(sprintf("[%s|%s] Only %d samples after TP53 subsetting (need >= 10), skip",
                      ds_id, stratum, length(sam_stratum)))
      return(invisible(NULL))
    }
    expr_mat <- expr_mat[, sam_stratum, drop = FALSE]
    all_sample_ids <- sam_stratum
  }

  message(sprintf("[%s|%s] Proteomics matrix: %d genes x %d samples",
                  ds_id, stratum, nrow(expr_mat), ncol(expr_mat)))

  # --- 5b. Identify CSN subunit rows in the proteomics matrix ---------------

  csn_row_map <- list()  # subunit_symbol -> row index in expr_mat
  for (sub in subunits) {
    ens_id <- sym2ens[sub]
    if (is.na(ens_id)) {
      message(sprintf("[%s|%s] CSN subunit %s: no Ensembl mapping available, skip",
                      ds_id, stratum, sub))
      next
    }
    # Match by base Ensembl ID (ignoring version suffix)
    hit <- which(ensembl_ids_base == ens_id)
    if (length(hit) == 0) {
      message(sprintf("[%s|%s] CSN subunit %s (%s): not found in proteomics data",
                      ds_id, stratum, sub, ens_id))
      next
    }
    csn_row_map[[sub]] <- hit[1]
  }

  if (length(csn_row_map) == 0) {
    message(sprintf("[%s|%s] No CSN subunits found in proteomics data, skip",
                    ds_id, stratum))
    return(invisible(NULL))
  }

  message(sprintf("[%s|%s] Found %d / %d CSN subunits: %s",
                  ds_id, stratum, length(csn_row_map), length(subunits),
                  paste(names(csn_row_map), collapse = ", ")))

  # --- 5c. Read covariates -------------------------------------------------

  ## Sex & Age from meta.txt
  meta_file <- file.path(ds_dir, paste0(ds_id, "_meta.txt"))
  meta_dt   <- data.table::fread(meta_file, header = TRUE, check.names = FALSE)
  # Remove the "data_type" descriptor row if present
  if (nrow(meta_dt) > 0 && meta_dt[[1]][1] == "data_type") {
    meta_dt <- meta_dt[-1, ]
  }
  meta_id_col <- names(meta_dt)[1]  # "case_id"
  meta_df <- as.data.frame(meta_dt, check.names = FALSE)
  rownames(meta_df) <- meta_df[[meta_id_col]]

  ## Tumor purity from phenotype.txt
  pheno_file <- file.path(ds_dir, paste0(ds_id, "_phenotype.txt"))
  pheno_dt   <- data.table::fread(pheno_file, header = TRUE, check.names = FALSE)
  pheno_id_col <- names(pheno_dt)[1]  # "idx"
  pheno_df <- as.data.frame(pheno_dt, check.names = FALSE)
  rownames(pheno_df) <- pheno_df[[pheno_id_col]]

  ## Build covariate data.frame aligned to sample IDs
  cov_df <- data.frame(row.names = all_sample_ids, check.names = FALSE)

  # Sex: encode as numeric (Female = 0, Male = 1)
  if ("Sex" %in% names(meta_df)) {
    sex_raw <- meta_df[all_sample_ids, "Sex"]
    sex_num <- rep(NA_real_, length(all_sample_ids))
    sex_num[tolower(sex_raw) == "female"] <- 0
    sex_num[tolower(sex_raw) == "male"]   <- 1
    cov_df$sex <- sex_num
  } else {
    message(sprintf("[%s|%s] Warning: 'Sex' column not found in meta.txt",
                    ds_id, stratum))
    cov_df$sex <- NA_real_
  }

  # Age: numeric
  if ("Age" %in% names(meta_df)) {
    age_raw <- meta_df[all_sample_ids, "Age"]
    cov_df$age <- suppressWarnings(as.numeric(age_raw))
  } else {
    message(sprintf("[%s|%s] Warning: 'Age' column not found in meta.txt",
                    ds_id, stratum))
    cov_df$age <- NA_real_
  }

  # Tumor purity (WES_purity): numeric
  if ("WES_purity" %in% names(pheno_df)) {
    pur_raw <- pheno_df[all_sample_ids, "WES_purity"]
    cov_df$purity <- suppressWarnings(as.numeric(pur_raw))
  } else {
    message(sprintf("[%s|%s] Warning: 'WES_purity' column not found in phenotype.txt",
                    ds_id, stratum))
    cov_df$purity <- NA_real_
  }

  ## Impute covariates
  imp_result  <- impute_covariates(cov_df, ds_id, stratum)
  cov_imputed <- imp_result$imputed_df
  imp_log     <- imp_result$imputation_log

  # Log imputation results
  message(sprintf("[%s|%s] Covariate imputation results:", ds_id, stratum))
  for (cv_name in names(imp_log)) {
    message(sprintf("  %s: %s", cv_name, imp_log[[cv_name]]))
  }

  # Save imputation log for this dataset and stratum
  imp_log_df <- data.frame(
    dataset   = ds_id,
    stratum   = stratum,
    covariate = names(imp_log),
    result    = unlist(imp_log),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  imp_log_file <- file.path(ds_out,
    paste0(ds_id, "_", stratum, "_covariate_imputation_log.csv"))
  data.table::fwrite(imp_log_df, imp_log_file)
  message(sprintf("[%s|%s] Imputation log saved: %s",
                  ds_id, stratum, basename(imp_log_file)))

  # --- 5d. Prepare usable covariates (drop all-NA or constant columns) -----

  cov_use <- cov_imputed
  drop_cols <- c()
  for (cn in names(cov_use)) {
    vals <- cov_use[[cn]]
    non_na_vals <- vals[!is.na(vals)]
    if (length(non_na_vals) == 0) {
      # All NA: drop
      drop_cols <- c(drop_cols, cn)
    } else if (length(unique(non_na_vals)) <= 1) {
      # Constant (e.g., single-sex cancer like BRCA, OV, UCEC): drop
      drop_cols <- c(drop_cols, cn)
    }
  }
  if (length(drop_cols) > 0) {
    cov_use <- cov_use[, !(names(cov_use) %in% drop_cols), drop = FALSE]
    message(sprintf("[%s|%s] Dropped constant/all-NA covariate columns: %s",
                    ds_id, stratum, paste(drop_cols, collapse = ", ")))
  }
  if (ncol(cov_use) > 0) {
    message(sprintf("[%s|%s] Covariates retained in model: %s",
                    ds_id, stratum, paste(names(cov_use), collapse = ", ")))
  } else {
    message(sprintf("[%s|%s] No covariates retained (all constant or all-NA)",
                    ds_id, stratum))
  }

  # --- 5e. Run limma for each CSN subunit ----------------------------------

  all_deg_results <- list()

  for (sub in names(csn_row_map)) {

    message(sprintf("\n  [%s|%s | %s] Running limma DEG analysis ...",
                    ds_id, stratum, sub))

    # Extract CSN subunit protein abundance as the continuous predictor
    csn_abund <- as.numeric(expr_mat[csn_row_map[[sub]], ])
    names(csn_abund) <- all_sample_ids

    # Identify samples with non-NA CSN abundance AND complete covariates
    ok_csn <- !is.na(csn_abund)
    if (ncol(cov_use) > 0) {
      ok_cov <- complete.cases(cov_use)
    } else {
      ok_cov <- rep(TRUE, length(all_sample_ids))
    }
    ok_samples <- ok_csn & ok_cov
    sam_use <- all_sample_ids[ok_samples]

    if (length(sam_use) < 10) {
      message(sprintf("  [%s|%s | %s] Only %d valid samples (need >= 10), skip",
                      ds_id, stratum, sub, length(sam_use)))
      next
    }

    message(sprintf("  [%s|%s | %s] Samples used: %d / %d (valid predictor + covariates)",
                    ds_id, stratum, sub, length(sam_use), length(all_sample_ids)))

    # Subset expression matrix and predictor to valid samples
    expr_sub <- expr_mat[, sam_use, drop = FALSE]
    csn_sub  <- csn_abund[sam_use]

    # Filter genes: require sufficient non-NA values across valid samples
    n_nonmissing <- rowSums(!is.na(expr_sub))
    threshold    <- max(min_n_gene, ceiling(length(sam_use) * min_prop))
    gene_keep    <- n_nonmissing >= threshold
    expr_sub     <- expr_sub[gene_keep, , drop = FALSE]

    message(sprintf("  [%s|%s | %s] Genes retained (>= %d non-NA values): %d / %d",
                    ds_id, stratum, sub, threshold,
                    sum(gene_keep), length(gene_keep)))

    if (nrow(expr_sub) == 0) {
      message(sprintf("  [%s|%s | %s] No genes passed the filter, skip",
                      ds_id, stratum, sub))
      next
    }

    # Build design matrix: ~ csn_abundance + covariates
    design_df <- data.frame(
      csn_abundance = csn_sub,
      row.names     = sam_use,
      check.names   = FALSE
    )
    if (ncol(cov_use) > 0) {
      for (cn in names(cov_use)) {
        design_df[[cn]] <- cov_use[sam_use, cn]
      }
    }

    design_mat <- model.matrix(~ ., data = design_df)

    # Safety check: rank deficiency
    qr_design <- qr(design_mat)
    if (qr_design$rank < ncol(design_mat)) {
      message(sprintf(
        "  [%s|%s | %s] Warning: design matrix rank-deficient (%d / %d cols). Removing dependent columns.",
        ds_id, stratum, sub, qr_design$rank, ncol(design_mat)))
      keep_cols  <- qr_design$pivot[seq_len(qr_design$rank)]
      design_mat <- design_mat[, keep_cols, drop = FALSE]
    }

    # Print design matrix column names for verification
    message(sprintf("  [%s|%s | %s] Design matrix columns: %s",
                    ds_id, stratum, sub,
                    paste(colnames(design_mat), collapse = ", ")))

    # Run limma: lmFit + eBayes
    fit <- limma::lmFit(expr_sub, design_mat)
    fit <- limma::eBayes(fit)

    # The coefficient of interest is "csn_abundance"
    coef_name <- "csn_abundance"
    if (!(coef_name %in% colnames(fit$coefficients))) {
      message(sprintf("  [%s|%s | %s] Coefficient '%s' not found in fitted model, skip",
                      ds_id, stratum, sub, coef_name))
      next
    }

    # Extract results: topTable with BH FDR adjustment
    tt <- limma::topTable(
      fit,
      coef         = coef_name,
      number       = Inf,
      sort.by      = "none",
      adjust.method = "BH"
    )

    # Build result data.frame
    matched_ensembl_raw  <- rownames(tt)
    matched_symbols      <- gene_symbols[matched_ensembl_raw]
    matched_ensembl_base <- sub("\\.[0-9]+$", "", matched_ensembl_raw)

    # Flag if the tested gene is the predictor CSN subunit itself
    predictor_ensembl_base <- sub("\\.[0-9]+$", "",
                                  ensembl_ids_raw[csn_row_map[[sub]]])
    is_self_predictor <- (matched_ensembl_base == predictor_ensembl_base)

    result_df <- data.frame(
      dataset          = ds_id,
      stratum          = stratum,
      csn_subunit      = sub,
      ensembl_id       = matched_ensembl_raw,
      ensembl_id_base  = matched_ensembl_base,
      gene_symbol      = matched_symbols,
      logFC            = tt$logFC,
      AveExpr          = tt$AveExpr,
      t_statistic      = tt$t,
      p_value          = tt$P.Value,
      BH_FDR           = tt$adj.P.Val,
      B_statistic      = tt$B,
      n_samples        = length(sam_use),
      is_predictor     = is_self_predictor,
      stringsAsFactors = FALSE,
      check.names      = FALSE
    )

    # Sort by p-value
    result_df <- result_df[order(result_df$p_value), ]
    rownames(result_df) <- NULL

    # Summary statistics
    n_tested   <- nrow(result_df)
    n_sig_005  <- sum(result_df$BH_FDR < 0.05, na.rm = TRUE)
    n_sig_01   <- sum(result_df$BH_FDR < 0.1, na.rm = TRUE)
    n_sig_025  <- sum(result_df$BH_FDR < 0.25, na.rm = TRUE)
    message(sprintf(
      "  [%s|%s | %s] Results: %d genes tested | BH FDR < 0.05: %d | < 0.1: %d | < 0.25: %d",
      ds_id, stratum, sub, n_tested, n_sig_005, n_sig_01, n_sig_025
    ))

    all_deg_results[[sub]] <- result_df

    # Save per-subunit CSV
    per_sub_file <- file.path(ds_out,
      paste0(ds_id, "_", stratum, "_limma_DEG_predictor_", sub, ".csv"))
    data.table::fwrite(result_df, per_sub_file)
    message(sprintf("  [%s|%s | %s] Saved: %s",
                    ds_id, stratum, sub, basename(per_sub_file)))
  }

  # --- 5f. Save combined results for this dataset + stratum ----------------

  if (length(all_deg_results) > 0) {
    combined <- do.call(rbind, all_deg_results)
    rownames(combined) <- NULL

    combined_file <- file.path(ds_out,
      paste0(ds_id, "_", stratum, "_limma_DEG_all_CSN_subunits.csv"))
    data.table::fwrite(combined, combined_file)
    message(sprintf("\n[%s|%s] Combined CSV saved: %s",
                    ds_id, stratum, basename(combined_file)))

    return(invisible(combined))
  } else {
    message(sprintf("\n[%s|%s] No DEG results generated for any CSN subunit",
                    ds_id, stratum))
    return(invisible(NULL))
  }
}


# ---- 6. Run all datasets with TP53 stratification ------------------------

message("\n============================================================")
message("  CSN Subunits Protein Abundance Limma DEG Analysis")
message("  Predictor: each CSN subunit's protein abundance (continuous)")
message("  Covariates: Sex, Age, Tumor Purity (WES_purity)")
message("  Stratification: TP53 mutation status")
message("  FDR method: Benjamini-Hochberg")
message("============================================================\n")

all_results <- list()

for (ds in DATASETS) {

  # Read TP53 classification for this dataset
  tp53_lists <- get_tp53_sample_lists(ds)

  # Determine which strata to run for this dataset
  if (ds %in% DS_SKIP_TP53_STRATIFICATION) {
    strata_to_run <- "TP53_all"
    message(sprintf("[%s] Skipping TP53 stratification (insufficient MUT or WT samples)",
                    ds))
  } else {
    strata_to_run <- TP53_STRATA
  }

  for (st in strata_to_run) {
    # Determine sample subset for each stratum
    if (st == "TP53_all") {
      sample_sub <- NULL  # use all samples in proteomics data
    } else if (!is.null(tp53_lists)) {
      sample_sub <- tp53_lists[[st]]
      if (is.null(sample_sub) || length(sample_sub) == 0) {
        message(sprintf("[%s|%s] No samples found for this stratum, skip", ds, st))
        next
      }
    } else {
      message(sprintf("[%s|%s] TP53 classification not available, skip", ds, st))
      next
    }

    res <- tryCatch(
      run_limma_deg_one_dataset(
        ds_id         = ds,
        stratum       = st,
        sample_subset = sample_sub
      ),
      error = function(e) {
        message(sprintf("[%s|%s] ERROR: %s", ds, st, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(res)) {
      key <- paste(ds, st, sep = "|")
      all_results[[key]] <- res
    }
  }
}


# ---- 7. Grand combined summary across all datasets and strata ------------

if (length(all_results) > 0) {

  grand_combined <- do.call(rbind, all_results)
  rownames(grand_combined) <- NULL

  # --- Per-stratum combined files ---
  for (st in TP53_STRATA) {
    st_data <- grand_combined[grand_combined$stratum == st, ]
    if (nrow(st_data) == 0) next

    st_dir <- file.path(OUT_ROOT, st)
    dir.create(st_dir, showWarnings = FALSE, recursive = TRUE)

    combined_csv <- file.path(st_dir,
      paste0("ALL_datasets_", st, "_limma_DEG_all_CSN_subunits.csv"))
    data.table::fwrite(st_data, combined_csv)
    message(sprintf("\n[%s] Combined CSV saved: %s", st, combined_csv))
  }

  # --- Grand combined file (all strata together) ---
  grand_file <- file.path(OUT_ROOT,
    "ALL_datasets_ALL_strata_limma_DEG_all_CSN_subunits.csv")
  data.table::fwrite(grand_combined, grand_file)
  message(sprintf("\nGrand combined CSV saved: %s", grand_file))

  # --- Summary table: number of significant DEGs per dataset x stratum x subunit ---
  summary_df <- grand_combined %>%
    dplyr::filter(!is_predictor) %>%
    dplyr::group_by(dataset, stratum, csn_subunit) %>%
    dplyr::summarise(
      n_genes_tested   = dplyr::n(),
      n_samples        = dplyr::first(n_samples),
      n_sig_FDR_0.05   = sum(BH_FDR < 0.05, na.rm = TRUE),
      n_sig_FDR_0.10   = sum(BH_FDR < 0.10, na.rm = TRUE),
      n_sig_FDR_0.25   = sum(BH_FDR < 0.25, na.rm = TRUE),
      n_up_FDR_0.05    = sum(BH_FDR < 0.05 & logFC > 0, na.rm = TRUE),
      n_down_FDR_0.05  = sum(BH_FDR < 0.05 & logFC < 0, na.rm = TRUE),
      .groups = "drop"
    )

  summary_file <- file.path(OUT_ROOT,
    "ALL_datasets_ALL_strata_limma_DEG_summary.csv")
  data.table::fwrite(summary_df, summary_file)
  message(sprintf("Summary CSV saved: %s", summary_file))

  # Per-stratum summary files
  for (st in TP53_STRATA) {
    st_summary <- summary_df[summary_df$stratum == st, ]
    if (nrow(st_summary) == 0) next
    st_dir <- file.path(OUT_ROOT, st)
    st_summary_file <- file.path(st_dir,
      paste0("ALL_datasets_", st, "_limma_DEG_summary.csv"))
    data.table::fwrite(st_summary, st_summary_file)
    message(sprintf("[%s] Summary CSV saved: %s", st, st_summary_file))
  }

  # Print summary table to console
  message("\n========== DEG Analysis Summary (excluding self-predictor) ==========")
  print(as.data.frame(summary_df), row.names = FALSE)
}


# ---- 8. Combine all imputation logs across datasets and strata -----------

imp_log_files <- list.files(
  OUT_ROOT,
  pattern     = "_covariate_imputation_log\\.csv$",
  recursive   = TRUE,
  full.names  = TRUE
)

if (length(imp_log_files) > 0) {
  all_imp_logs <- lapply(imp_log_files, data.table::fread)
  combined_imp <- do.call(rbind, all_imp_logs)

  imp_combined_file <- file.path(OUT_ROOT,
    "ALL_datasets_ALL_strata_covariate_imputation_log.csv")
  data.table::fwrite(combined_imp, imp_combined_file)
  message(sprintf("Combined imputation log saved: %s", imp_combined_file))
}


message("\n============================================================")
message("  All analyses completed!")
message(sprintf("  Output directory: %s", OUT_ROOT))
message("============================================================\n")
