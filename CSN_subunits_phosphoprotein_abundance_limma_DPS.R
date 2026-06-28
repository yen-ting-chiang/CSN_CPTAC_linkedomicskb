## =========================================================================
##  CSN Subunits Protein Abundance Limma DPS Analysis
##  Data source: LinkedOmicsKB CPTAC datasets
##
##  Purpose:
##    For each CPTAC dataset and each CSN subunit, use the subunit's protein
##    abundance as a continuous predictor in limma to perform genome-wide
##    Differentially Phosphorylated Sites (DPS) analysis across the
##    phosphoproteome. Additionally includes CSN_SCORE (PCA PC1 of z-scored
##    CSN subunits, with COPS7A/7B combined) as a composite predictor, and
##    resid_* predictors (individual subunit with CSN_SCORE partialled out
##    as an additional covariate) to isolate each subunit's unique
##    contribution.
##
##  Phosphosite normalization:
##    Phosphosite abundances (log2) are normalized by subtracting the
##    corresponding protein's abundance (log2) on a per-sample basis:
##      normalized_phospho = log2(phospho) - log2(protein)
##    This is the standard approach used in CPTAC proteogenomics studies
##    to remove the confounding effect of protein abundance changes on
##    phosphorylation level estimates.
##    References:
##      Zhang B et al. (2014) Nature 513:382-387.
##        https://doi.org/10.1038/nature13438
##      Mertins P et al. (2016) Nature 534:55-62.
##        https://doi.org/10.1038/nature18003
##      Gillette MA et al. (2020) Cell 182:200-225.
##        https://doi.org/10.1016/j.cell.2020.06.013
##
##  Statistical rationale:
##    Using a continuous predictor in limma is a well-established approach
##    (limma User's Guide Section 9.4). It identifies phospho sites whose
##    normalized abundance changes linearly with the CSN subunit level,
##    after adjusting for covariates. limma's empirical Bayes moderation
##    stabilizes variance estimates and yields well-calibrated p-values.
##    Reference: Ritchie ME et al. (2015) Nucleic Acids Res 43(7):e47.
##               https://doi.org/10.1093/nar/gkv007
##
##  Covariates adjusted: Sex, Age, Tumor Purity (WES_purity)
##  Multiple testing correction: Benjamini-Hochberg FDR
##  Phosphosite data: NOT imputed (sites with excessive missing values are
##                    filtered out; limma handles remaining NAs internally)
##  Covariate data: Imputed per covariate using the most appropriate method:
##                  - Sex (categorical): mode imputation
##                  - Age (continuous):  median imputation
##                  - Tumor purity (continuous): median imputation
##  Stratification: TP53 mutation status (TP53_all, TP53_MUT, TP53_WT,
##                  and TP53_interaction).
##                  TP53_interaction tests the predictor * TP53_status
##                  interaction term in the union of MUT and WT samples.
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

# Toggle whether to re-run the time-consuming main limma analysis (Step 6)
# If FALSE, the script skips Step 6 and loads previously saved results for Step 7.
RE_RUN_MAIN_ANALYSIS <- FALSE

# Export configurations for combined files (default to FALSE as requested)
EXPORT_PER_STRATUM_COMBINED <- FALSE
EXPORT_GRAND_COMBINED       <- FALSE

BASE_DIR <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb"
OUT_ROOT <- file.path(BASE_DIR, "CSN_subunits_phosphoprotein_abundance_limma_DPS")
TP53_DIR <- file.path(BASE_DIR, "TP53_mutation_classification")
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

DATASETS <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC",
               "LSCC", "LUAD", "OV", "PDAC", "UCEC")

CSN_SUBUNITS <- c("GPS1", "COPS2", "COPS3", "COPS4", "COPS5",
                   "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9")

# TP53 stratification strata
TP53_STRATA <- c("TP53_all", "TP53_MUT", "TP53_WT", "TP53_interaction")

# Datasets excluded from TP53 stratification:
#   CCRCC: too few TP53 mutant samples (5/103)
#   LSCC:  too few TP53 wild-type samples (5/108)
#   OV:    too few TP53 wild-type samples (5/82)
DS_SKIP_TP53_STRATIFICATION <- c("CCRCC", "LSCC", "OV")

# Minimum number of non-NA samples per phospho site to include in limma
MIN_SAMPLES_PER_GENE <- 10L

# Minimum proportion of non-NA values per phospho site (relative to valid samples)
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

# ---- 4b. Helpers: CSN_SCORE computation (PCA PC1 of z-scored subunits) ---

## Clean matrix for PCA: remove rows with too few finite values, impute
## remaining NAs with row medians
.clean_for_pca <- function(X, min_samples = 10L, min_genes = 5L) {
  X <- as.matrix(X)
  X[!is.finite(X)] <- NA
  keep_rows <- rowSums(is.finite(X)) >= min_samples
  if (!any(keep_rows)) {
    return(NULL)
  }
  X <- X[keep_rows, , drop = FALSE]
  if (nrow(X) < min_genes) {
    return(NULL)
  }
  if (anyNA(X)) {
    med <- apply(X, 1, function(r) median(r[is.finite(r)], na.rm = TRUE))
    for (i in seq_len(nrow(X))) {
      xi <- X[i, ]
      xi[!is.finite(xi)] <- med[i]
      X[i, ] <- xi
    }
  }
  X
}

## Perform PCA using the cross-sample z-values of subunits, and take PC1
## as the CSN score; adjust direction to have the same sign as the mean z.
build_csn_score <- function(mat0,
                            subunits = CSN_SUBUNITS,
                            combine_7AB = TRUE,
                            min_members = 5L) {
  present <- intersect(subunits, rownames(mat0))
  # Pre-create the return skeleton
  s <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
  if (!length(present)) {
    return(s)
  }

  # z-score (preserve sample names)
  get_z <- function(v) {
    nm <- names(v)
    v <- as.numeric(v)
    mu <- mean(v[is.finite(v)], na.rm = TRUE)
    sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    v[!is.finite(v)] <- mu
    out <- (v - mu) / sdv
    names(out) <- nm
    out
  }

  X <- do.call(rbind, lapply(present, function(g) get_z(mat0[g, ])))
  rownames(X) <- present
  colnames(X) <- colnames(mat0)

  # combine COPS7A/7B
  if (combine_7AB && all(c("COPS7A", "COPS7B") %in% rownames(X))) {
    Z7 <- colMeans(X[c("COPS7A", "COPS7B"), , drop = FALSE], na.rm = TRUE)
    X <- rbind(X[setdiff(rownames(X), c("COPS7A", "COPS7B")), , drop = FALSE],
      "COPS7*" = Z7
    )
  }

  enough <- colSums(is.finite(mat0[present, , drop = FALSE])) >= min_members
  keep_sam <- names(s)[enough]

  if (length(keep_sam) >= 10) {
    pc <- stats::prcomp(t(X[, keep_sam, drop = FALSE]), center = TRUE, scale. = FALSE)
    sc <- pc$x[, 1]
    # Direction correction: Same sign as subunit average z
    mu <- colMeans(X[, keep_sam, drop = FALSE], na.rm = TRUE)
    if (suppressWarnings(cor(sc, mu, use = "pairwise.complete.obs")) < 0) sc <- -sc
    s[keep_sam] <- sc
  }

  s
}

## Safe version: wraps build_csn_score with error handling and fallback
build_csn_score_safe <- function(mat0, subunits, combine_7AB = TRUE,
                                 min_members = 5L, pca_min_samples = 10L) {
  sub <- intersect(subunits, rownames(mat0))
  out_na <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
  if (length(sub) < min_members) {
    return(out_na)
  }

  cs_try <- try(
    {
      build_csn_score(mat0, subunits = sub, combine_7AB = combine_7AB, min_members = min_members)
    },
    silent = TRUE
  )

  if (!inherits(cs_try, "try-error") && sum(is.finite(cs_try)) >= pca_min_samples) {
    return(cs_try)
  }

  X <- .clean_for_pca(mat0[sub, , drop = FALSE], min_samples = pca_min_samples, min_genes = min_members)
  if (is.null(X)) {
    message("[CSN_SCORE-safe] Insufficient available subunits or samples; returning all NA")
    return(out_na)
  }
  pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
  if (inherits(pc, "try-error")) {
    message("[CSN_SCORE-safe] prcomp failed; returning all NA")
    return(out_na)
  }
  sc <- pc$x[, 1]
  names(sc) <- rownames(pc$x)
  ref <- colMeans(X, na.rm = TRUE)
  rr <- suppressWarnings(stats::cor(sc, ref, use = "pairwise.complete.obs"))
  if (is.finite(rr) && rr < 0) sc <- -sc
  out <- out_na
  out[names(sc)] <- as.numeric(sc)

  varpc1 <- if (!is.null(pc$sdev)) (pc$sdev[1]^2) / sum(pc$sdev^2) else NA_real_
  message(sprintf("[CSN_SCORE-safe] fallback: genes=%d; PC1%%=%.1f; nonNA=%d/%d",
                  nrow(X), 100 * varpc1, sum(is.finite(out)), length(out)))
  out
}

# ---- 5. Main function: limma DPS for one dataset + stratum ---------------

run_limma_dps_one_dataset <- function(
    ds_id,
    stratum       = "TP53_all",
    sample_subset = NULL,
    tp53_lists    = NULL,
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

  # --- 5a. Read proteomics matrix (for predictor extraction + phospho normalization) ---

  prot_file <- file.path(ds_dir,
    paste0(ds_id,
           "_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"))

  if (!file.exists(prot_file)) {
    message(sprintf("[%s|%s] Proteomics file not found: %s",
                    ds_id, stratum, basename(prot_file)))
    return(invisible(NULL))
  }

  prot_dt <- data.table::fread(prot_file, header = TRUE, check.names = FALSE)
  prot_idx_col <- names(prot_dt)[1]  # "idx"
  prot_ensembl_ids_raw <- prot_dt[[prot_idx_col]]

  # Strip version suffix from Ensembl IDs (e.g., ENSG00000008083.17)
  prot_ensembl_ids_base <- sub("\\.[0-9]+$", "", prot_ensembl_ids_raw)

  # Sample IDs from proteomics data
  prot_sample_ids <- setdiff(names(prot_dt), prot_idx_col)

  # Build numeric protein matrix (genes x samples)
  prot_mat <- as.matrix(prot_dt[, ..prot_sample_ids])
  rownames(prot_mat) <- prot_ensembl_ids_raw
  storage.mode(prot_mat) <- "double"

  # Create a lookup: Ensembl base ID -> row index in prot_mat (for normalization)
  prot_base_to_row <- setNames(seq_along(prot_ensembl_ids_base), prot_ensembl_ids_base)
  # Handle duplicates: keep first occurrence
  prot_base_to_row <- prot_base_to_row[!duplicated(names(prot_base_to_row))]

  message(sprintf("[%s|%s] Proteomics matrix: %d genes x %d samples",
                  ds_id, stratum, nrow(prot_mat), ncol(prot_mat)))

  # --- 5a2. Read phospho site abundance matrix ---

  phospho_file <- file.path(ds_dir,
    paste0(ds_id,
           "_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt"))

  if (!file.exists(phospho_file)) {
    message(sprintf("[%s|%s] Phospho site file not found: %s",
                    ds_id, stratum, basename(phospho_file)))
    return(invisible(NULL))
  }

  phospho_dt <- data.table::fread(phospho_file, header = TRUE, check.names = FALSE)
  phospho_idx_col <- names(phospho_dt)[1]  # "idx"
  phospho_ids_raw <- phospho_dt[[phospho_idx_col]]

  # Parse phospho site IDs:
  #   Format: ENSG00000048028.11|ENSP00000003302.4|S1053|PPTIRPNSPYDLCSR|1
  #   Fields: ensembl_gene | ensembl_protein | phospho_site | peptide | number
  phospho_gene_ids <- sub("\\|.*$", "", phospho_ids_raw)  # extract ENSG with version
  phospho_gene_base <- sub("\\.[0-9]+$", "", phospho_gene_ids)  # strip version

  # Extract phospho site annotation (e.g., "S1053")
  phospho_site_info <- sub("^[^|]+\\|[^|]+\\|", "", phospho_ids_raw)  # remove first two fields
  phospho_site_label <- sub("\\|.*$", "", phospho_site_info)  # extract site label

  # Map phospho site Ensembl gene IDs to gene symbols
  phospho_gene_symbols <- ens2sym[phospho_gene_base]
  phospho_gene_symbols[is.na(phospho_gene_symbols)] <- phospho_gene_ids[is.na(phospho_gene_symbols)]

  # Use common samples between phospho and proteomics data
  phospho_sample_ids <- setdiff(names(phospho_dt), phospho_idx_col)
  common_samples <- intersect(phospho_sample_ids, prot_sample_ids)

  if (length(common_samples) < 10) {
    message(sprintf("[%s|%s] Only %d common samples between phospho and proteomics (need >= 10), skip",
                    ds_id, stratum, length(common_samples)))
    return(invisible(NULL))
  }

  # Build numeric phospho matrix (phospho sites x common samples)
  phospho_mat <- as.matrix(phospho_dt[, ..common_samples])
  rownames(phospho_mat) <- phospho_ids_raw
  storage.mode(phospho_mat) <- "double"

  message(sprintf("[%s|%s] Phospho site matrix: %d sites x %d samples (common with proteomics)",
                  ds_id, stratum, nrow(phospho_mat), length(common_samples)))

  # --- 5a3. Normalize phospho by protein abundance -------------------------
  #   For log2 data: normalized = log2(phospho) - log2(protein)
  #   This removes the confounding effect of total protein level changes.

  message(sprintf("[%s|%s] Normalizing phospho site abundance by protein abundance ...",
                  ds_id, stratum))

  # Subset protein matrix to common samples
  prot_mat_common <- prot_mat[, common_samples, drop = FALSE]

  n_normalized <- 0L
  n_no_protein <- 0L

  for (i in seq_len(nrow(phospho_mat))) {
    gene_base <- phospho_gene_base[i]
    if (gene_base %in% names(prot_base_to_row)) {
      prot_row <- prot_base_to_row[[gene_base]]
      prot_vals <- prot_mat_common[prot_row, ]
      # Subtract protein abundance from phospho abundance (both log2)
      # Result is NA if either phospho or protein is NA
      phospho_mat[i, ] <- phospho_mat[i, ] - prot_vals
      n_normalized <- n_normalized + 1L
    } else {
      # No matching protein found: set to NA (cannot normalize)
      phospho_mat[i, ] <- NA_real_
      n_no_protein <- n_no_protein + 1L
    }
  }

  message(sprintf("[%s|%s] Phospho normalization: %d / %d sites normalized; %d sites set to NA (no matching protein)",
                  ds_id, stratum, n_normalized, nrow(phospho_mat), n_no_protein))

  # Remove sites that are entirely NA after normalization
  all_na_rows <- rowSums(!is.na(phospho_mat)) == 0
  if (any(all_na_rows)) {
    phospho_mat <- phospho_mat[!all_na_rows, , drop = FALSE]
    phospho_ids_raw <- phospho_ids_raw[!all_na_rows]
    phospho_gene_ids <- phospho_gene_ids[!all_na_rows]
    phospho_gene_base <- phospho_gene_base[!all_na_rows]
    phospho_site_label <- phospho_site_label[!all_na_rows]
    phospho_gene_symbols <- phospho_gene_symbols[!all_na_rows]
    message(sprintf("[%s|%s] Removed %d all-NA sites after normalization; %d sites remaining",
                    ds_id, stratum, sum(all_na_rows), nrow(phospho_mat)))
  }

  # Use common_samples as the working sample set
  all_sample_ids <- common_samples

  # --- 5a4. Subset samples by TP53 stratum ---------------------------------

  if (!is.null(sample_subset)) {
    sam_stratum <- intersect(sample_subset, all_sample_ids)
    if (length(sam_stratum) < 10) {
      message(sprintf("[%s|%s] Only %d samples after TP53 subsetting (need >= 10), skip",
                      ds_id, stratum, length(sam_stratum)))
      return(invisible(NULL))
    }
    phospho_mat <- phospho_mat[, sam_stratum, drop = FALSE]
    prot_mat_common <- prot_mat_common[, sam_stratum, drop = FALSE]
    all_sample_ids <- sam_stratum
  }

  message(sprintf("[%s|%s] Normalized phospho matrix after subsetting: %d sites x %d samples",
                  ds_id, stratum, nrow(phospho_mat), ncol(phospho_mat)))

  # --- 5b. Identify CSN subunit rows in the proteomics matrix ---------------

  csn_row_map <- list()  # subunit_symbol -> row index in prot_mat
  for (sub in subunits) {
    ens_id <- sym2ens[sub]
    if (is.na(ens_id)) {
      message(sprintf("[%s|%s] CSN subunit %s: no Ensembl mapping available, skip",
                      ds_id, stratum, sub))
      next
    }
    # Match by base Ensembl ID (ignoring version suffix)
    hit <- which(prot_ensembl_ids_base == ens_id)
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

  # --- 5b2. Build gene-symbol-keyed matrix for CSN_SCORE computation -------
  # mat0_csn: rows = gene symbols of present CSN subunits, cols = samples
  # Use protein abundance (NOT phospho) for CSN_SCORE computation
  present_sub <- names(csn_row_map)
  mat0_csn <- prot_mat_common[unlist(csn_row_map), , drop = FALSE]
  rownames(mat0_csn) <- present_sub
  storage.mode(mat0_csn) <- "double"

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

  # --- 5e. Build predictor list: individual CSN subunits + CSN_SCORE -------
  # Predictors are protein abundances (NOT phospho-normalized values)

  # Compute CSN_SCORE via PCA PC1 of z-scored CSN subunit abundances
  csn_score_vec <- build_csn_score_safe(
    mat0         = mat0_csn,
    subunits     = present_sub,
    combine_7AB  = TRUE,
    min_members  = 5L,
    pca_min_samples = 10L
  )

  csn_score_valid <- sum(is.finite(csn_score_vec))
  message(sprintf("[%s|%s] CSN_SCORE computed: %d / %d samples with valid score",
                  ds_id, stratum, csn_score_valid, length(csn_score_vec)))

  # Build predictor list: each individual CSN subunit + CSN_SCORE
  # Use protein abundance from prot_mat_common (aligned to all_sample_ids)
  predictor_list <- list()
  for (sub in present_sub) {
    v <- as.numeric(prot_mat_common[csn_row_map[[sub]], ])
    names(v) <- all_sample_ids
    predictor_list[[sub]] <- v
  }
  # Add CSN_SCORE only if it has enough valid values
  csn_score_available <- (csn_score_valid >= 10)
  if (csn_score_available) {
    predictor_list[["CSN_SCORE"]] <- csn_score_vec
  } else {
    message(sprintf("[%s|%s] CSN_SCORE has < 10 valid samples, skipping CSN_SCORE predictor",
                    ds_id, stratum))
  }

  # Add resid_* predictors: same subunit abundance as predictor, but with
  # CSN_SCORE added as an additional covariate to partial out the shared
  # CSN complex activity. Only available when CSN_SCORE is valid.
  if (csn_score_available) {
    for (sub in present_sub) {
      resid_name <- paste0("resid_", sub)
      predictor_list[[resid_name]] <- predictor_list[[sub]]
    }
    message(sprintf("[%s|%s] Added %d resid_* predictors (with CSN_SCORE as covariate)",
                    ds_id, stratum, length(present_sub)))
  } else {
    message(sprintf("[%s|%s] Skipping resid_* predictors (CSN_SCORE not available)",
                    ds_id, stratum))
  }

  # Collect Ensembl base IDs of all CSN subunit genes for is_predictor flagging
  all_csn_ensembl_base <- sub("\\.[0-9]+$", "",
                              prot_ensembl_ids_raw[unlist(csn_row_map)])

  # --- 5f. Run limma for each predictor ------------------------------------

  all_dps_results <- list()

  for (pred_name in names(predictor_list)) {

    message(sprintf("\n  [%s|%s | %s] Running limma DPS analysis ...",
                    ds_id, stratum, pred_name))

    # Detect if this is a resid_* predictor (CSN_SCORE as additional covariate)
    is_resid_predictor <- startsWith(pred_name, "resid_")

    # Extract predictor vector
    pred_vec <- predictor_list[[pred_name]]

    # Identify samples with non-NA predictor AND complete covariates
    # For resid_* predictors, also require non-NA CSN_SCORE
    ok_pred <- !is.na(pred_vec)
    if (is_resid_predictor) {
      ok_pred <- ok_pred & is.finite(csn_score_vec)
    }
    if (ncol(cov_use) > 0) {
      ok_cov <- complete.cases(cov_use)
    } else {
      ok_cov <- rep(TRUE, length(all_sample_ids))
    }
    ok_samples <- ok_pred & ok_cov
    sam_use <- all_sample_ids[ok_samples]

    if (length(sam_use) < 10) {
      message(sprintf("  [%s|%s | %s] Only %d valid samples (need >= 10), skip",
                      ds_id, stratum, pred_name, length(sam_use)))
      next
    }

    message(sprintf("  [%s|%s | %s] Samples used: %d / %d (valid predictor + covariates)",
                    ds_id, stratum, pred_name, length(sam_use), length(all_sample_ids)))

    # Subset normalized phospho matrix and predictor to valid samples
    phospho_sub <- phospho_mat[, sam_use, drop = FALSE]
    pred_sub <- pred_vec[sam_use]

    # Filter phospho sites: require sufficient non-NA values across valid samples
    n_nonmissing <- rowSums(!is.na(phospho_sub))
    threshold    <- max(min_n_gene, ceiling(length(sam_use) * min_prop))
    site_keep    <- n_nonmissing >= threshold
    phospho_sub  <- phospho_sub[site_keep, , drop = FALSE]

    # Also subset the annotation vectors for kept sites
    kept_ids_raw     <- phospho_ids_raw[site_keep]
    kept_gene_ids    <- phospho_gene_ids[site_keep]
    kept_gene_base   <- phospho_gene_base[site_keep]
    kept_site_label  <- phospho_site_label[site_keep]
    kept_gene_symbols <- phospho_gene_symbols[site_keep]

    message(sprintf("  [%s|%s | %s] Phospho sites retained (>= %d non-NA values): %d / %d",
                    ds_id, stratum, pred_name, threshold,
                    sum(site_keep), length(site_keep)))

    if (nrow(phospho_sub) == 0) {
      message(sprintf("  [%s|%s | %s] No phospho sites passed the filter, skip",
                      ds_id, stratum, pred_name))
      next
    }

    # Build design matrix base
    design_df <- data.frame(
      csn_abundance = pred_sub,
      row.names     = sam_use,
      check.names   = FALSE
    )
    # For resid_* predictors, add CSN_SCORE as an additional covariate
    if (is_resid_predictor) {
      design_df$csn_score_cov <- csn_score_vec[sam_use]
    }
    # For TP53_interaction, add TP53_status
    if (stratum == "TP53_interaction") {
      tp53_stat <- rep(NA_character_, length(sam_use))
      tp53_stat[sam_use %in% tp53_lists$TP53_MUT] <- "MUT"
      tp53_stat[sam_use %in% tp53_lists$TP53_WT]  <- "WT"
      design_df$TP53_status <- factor(tp53_stat, levels = c("WT", "MUT")) # WT is reference
    }
    if (ncol(cov_use) > 0) {
      for (cn in names(cov_use)) {
        design_df[[cn]] <- cov_use[sam_use, cn]
      }
    }

    # Construct formula dynamically
    form_str <- "~ csn_abundance"
    if (stratum == "TP53_interaction") {
      form_str <- "~ csn_abundance * TP53_status"
    }
    if (is_resid_predictor) {
      form_str <- paste(form_str, "+ csn_score_cov")
    }
    if (ncol(cov_use) > 0) {
      form_str <- paste(form_str, "+", paste(names(cov_use), collapse = " + "))
    }

    design_mat <- model.matrix(as.formula(form_str), data = design_df)

    # Safety check: rank deficiency
    qr_design <- qr(design_mat)
    if (qr_design$rank < ncol(design_mat)) {
      message(sprintf(
        "  [%s|%s | %s] Warning: design matrix rank-deficient (%d / %d cols). Removing dependent columns.",
        ds_id, stratum, pred_name, qr_design$rank, ncol(design_mat)))
      keep_cols  <- qr_design$pivot[seq_len(qr_design$rank)]
      design_mat <- design_mat[, keep_cols, drop = FALSE]
    }

    # Print design matrix column names for verification
    message(sprintf("  [%s|%s | %s] Design matrix columns: %s",
                    ds_id, stratum, pred_name,
                    paste(colnames(design_mat), collapse = ", ")))

    # Run limma: lmFit + eBayes
    fit <- limma::lmFit(phospho_sub, design_mat)
    fit <- limma::eBayes(fit)

    # Determine coefficient of interest
    if (stratum == "TP53_interaction") {
      # We want the interaction term (csn_abundance:TP53_statusMUT)
      coef_match <- grep("csn_abundance:TP53_status|TP53_status.*:csn_abundance", colnames(fit$coefficients), value = TRUE)
      if (length(coef_match) > 0) {
        coef_name <- coef_match[1]
      } else {
        coef_name <- "MISSING_INTERACTION"
      }
    } else {
      coef_name <- "csn_abundance"
    }

    if (!(coef_name %in% colnames(fit$coefficients))) {
      message(sprintf("  [%s|%s | %s] Coefficient '%s' not found in fitted model, skip",
                      ds_id, stratum, pred_name, coef_name))
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
    matched_phospho_ids  <- rownames(tt)
    # Find the indices of matched phospho sites in the kept annotation vectors
    matched_idx <- match(matched_phospho_ids, kept_ids_raw)

    matched_gene_ids     <- kept_gene_ids[matched_idx]
    matched_gene_base    <- kept_gene_base[matched_idx]
    matched_site_label   <- kept_site_label[matched_idx]
    matched_gene_symbols <- kept_gene_symbols[matched_idx]

    # Flag is_predictor:
    #   For individual subunits: flag phospho sites from the predictor gene
    #   For CSN_SCORE: flag phospho sites from ALL CSN subunit genes
    #   For resid_*: flag phospho sites from the corresponding subunit gene
    if (pred_name == "CSN_SCORE") {
      is_self_predictor <- (matched_gene_base %in% all_csn_ensembl_base)
    } else if (is_resid_predictor) {
      base_sub <- sub("^resid_", "", pred_name)
      predictor_ensembl_base <- sub("\\.[0-9]+$", "",
                                    prot_ensembl_ids_raw[csn_row_map[[base_sub]]])
      is_self_predictor <- (matched_gene_base == predictor_ensembl_base)
    } else {
      predictor_ensembl_base <- sub("\\.[0-9]+$", "",
                                    prot_ensembl_ids_raw[csn_row_map[[pred_name]]])
      is_self_predictor <- (matched_gene_base == predictor_ensembl_base)
    }

    result_df <- data.frame(
      dataset          = ds_id,
      stratum          = stratum,
      csn_subunit      = pred_name,
      phospho_site_id  = matched_phospho_ids,
      ensembl_gene_id  = matched_gene_ids,
      ensembl_gene_base = matched_gene_base,
      gene_symbol      = matched_gene_symbols,
      phospho_site     = matched_site_label,
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
      "  [%s|%s | %s] Results: %d sites tested | BH FDR < 0.05: %d | < 0.1: %d | < 0.25: %d",
      ds_id, stratum, pred_name, n_tested, n_sig_005, n_sig_01, n_sig_025
    ))

    all_dps_results[[pred_name]] <- result_df

    # Save per-predictor CSV
    per_sub_file <- file.path(ds_out,
      paste0(ds_id, "_", stratum, "_limma_DPS_predictor_", pred_name, ".csv"))
    data.table::fwrite(result_df, per_sub_file)
    message(sprintf("  [%s|%s | %s] Saved: %s",
                    ds_id, stratum, pred_name, basename(per_sub_file)))
  }

  # --- 5g. Save combined results for this dataset + stratum ----------------

  if (length(all_dps_results) > 0) {
    combined <- do.call(rbind, all_dps_results)
    rownames(combined) <- NULL

    combined_file <- file.path(ds_out,
      paste0(ds_id, "_", stratum, "_limma_DPS_all_predictors.csv"))
    data.table::fwrite(combined, combined_file)
    message(sprintf("\n[%s|%s] Combined CSV saved: %s",
                    ds_id, stratum, basename(combined_file)))

    return(invisible(combined))
  } else {
    message(sprintf("\n[%s|%s] No DPS results generated for any predictor",
                    ds_id, stratum))
    return(invisible(NULL))
  }
}


# ---- 6. Run all datasets with TP53 stratification ------------------------

all_results <- list()

if (RE_RUN_MAIN_ANALYSIS) {
  message("\n============================================================")
  message("  CSN Subunits Protein Abundance Limma DPS Analysis")
  message("  Outcome: Phospho site abundance (normalized by protein abundance)")
  message("  Predictors: each CSN subunit's protein abundance + CSN_SCORE (continuous)")
  message("  Covariates: Sex, Age, Tumor Purity (WES_purity)")
  message("  Stratification: TP53 mutation status")
  message("  FDR method: Benjamini-Hochberg")
  message("============================================================\n")

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
    } else if (st == "TP53_interaction") {
      if (!is.null(tp53_lists) && length(tp53_lists$TP53_MUT) > 0 && length(tp53_lists$TP53_WT) > 0) {
        sample_sub <- union(tp53_lists$TP53_MUT, tp53_lists$TP53_WT)
      } else {
        message(sprintf("[%s|%s] Missing MUT or WT samples, skip", ds, st))
        next
      }
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
      run_limma_dps_one_dataset(
        ds_id         = ds,
        stratum       = st,
        sample_subset = sample_sub,
        tp53_lists    = tp53_lists
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
} else {
  message("\n[Skip] RE_RUN_MAIN_ANALYSIS is FALSE. Skipping main limma execution.")
}


# ---- 7. Grand combined summary across all datasets and strata ------------

if (!RE_RUN_MAIN_ANALYSIS) {
  message("\nLoading previously saved predictor results from disk for summary...")
  saved_files <- list.files(
    OUT_ROOT,
    pattern     = "_limma_DPS_predictor_.*\\.csv$",
    recursive   = TRUE,
    full.names  = TRUE
  )
  if (length(saved_files) > 0) {
    message(sprintf("Found %d per-predictor CSV files. Reading and combining...", length(saved_files)))
    all_results <- lapply(saved_files, data.table::fread)
  } else {
    message("Warning: No previously saved per-predictor CSV files found in OUT_ROOT!")
  }
}

if (length(all_results) > 0) {

  grand_combined <- do.call(rbind, all_results)
  rownames(grand_combined) <- NULL

  # --- Per-stratum combined files ---
  if (EXPORT_PER_STRATUM_COMBINED) {
    for (st in TP53_STRATA) {
      st_data <- grand_combined[grand_combined$stratum == st, ]
      if (nrow(st_data) == 0) next
  
      st_dir <- file.path(OUT_ROOT, st)
      dir.create(st_dir, showWarnings = FALSE, recursive = TRUE)
  
      combined_csv <- file.path(st_dir,
        paste0("ALL_datasets_", st, "_limma_DPS_all_predictors.csv"))
      data.table::fwrite(st_data, combined_csv)
      message(sprintf("\n[%s] Combined CSV saved: %s", st, combined_csv))
    }
  }

  # --- Grand combined file (all strata together) ---
  if (EXPORT_GRAND_COMBINED) {
    grand_file <- file.path(OUT_ROOT,
      "ALL_datasets_ALL_strata_limma_DPS_all_predictors.csv")
    data.table::fwrite(grand_combined, grand_file)
    message(sprintf("\nGrand combined CSV saved: %s", grand_file))
  }

  # --- Summary table: number of significant DPS per dataset x stratum x subunit ---
  summary_df <- grand_combined %>%
    dplyr::filter(!is_predictor) %>%
    dplyr::group_by(dataset, stratum, csn_subunit) %>%
    dplyr::summarise(
      n_sites_tested   = dplyr::n(),
      n_samples        = dplyr::first(n_samples),
      n_sig_FDR_0.05   = sum(BH_FDR < 0.05, na.rm = TRUE),
      n_sig_FDR_0.10   = sum(BH_FDR < 0.10, na.rm = TRUE),
      n_sig_FDR_0.25   = sum(BH_FDR < 0.25, na.rm = TRUE),
      n_up_FDR_0.05    = sum(BH_FDR < 0.05 & logFC > 0, na.rm = TRUE),
      n_down_FDR_0.05  = sum(BH_FDR < 0.05 & logFC < 0, na.rm = TRUE),
      .groups = "drop"
    )

  summary_file <- file.path(OUT_ROOT,
    "ALL_datasets_ALL_strata_limma_DPS_summary.csv")
  data.table::fwrite(summary_df, summary_file)
  message(sprintf("Summary CSV saved: %s", summary_file))

  # Per-stratum summary files
  for (st in TP53_STRATA) {
    st_summary <- summary_df[summary_df$stratum == st, ]
    if (nrow(st_summary) == 0) next
    st_dir <- file.path(OUT_ROOT, st)
    st_summary_file <- file.path(st_dir,
      paste0("ALL_datasets_", st, "_limma_DPS_summary.csv"))
    data.table::fwrite(st_summary, st_summary_file)
    message(sprintf("[%s] Summary CSV saved: %s", st, st_summary_file))
  }

  # Print summary table to console
  message("\n========== DPS Analysis Summary (excluding self-predictor) ==========")
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
message("  All DPS analyses completed!")
message(sprintf("  Output directory: %s", OUT_ROOT))
message("============================================================\n")
