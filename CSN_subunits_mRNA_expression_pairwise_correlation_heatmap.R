## =========================================================================
##  CSN Subunits mRNA Expression Pairwise Correlation & Heatmap
##  Data source: LinkedOmicsKB CPTAC datasets
##  mRNA data: RNAseq gene RSEM coding UQ 1500 log2 (Tumor)
##  Covariates removed: Sex, Age, Tumor Purity (WES_purity)
##  Correlation method: Pearson correlation coefficient
##  Multiple testing correction: Benjamini-Hochberg FDR
##  Stratification: TP53 mutation status (TP53_all, TP53_MUT, TP53_WT)
## =========================================================================

# ---- 0. Load / install required packages ---------------------------------

required_pkgs <- c(
  "data.table", "dplyr", "tidyr", "readr", "stringr",
  "ggplot2", "scales", "openxlsx",
  "AnnotationDbi", "org.Hs.eg.db"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("AnnotationDbi", "org.Hs.eg.db")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(openxlsx)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ---- 1. Global configuration --------------------------------------------

BASE_DIR  <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb"
OUT_ROOT  <- file.path(BASE_DIR, "CSN_subunits_mRNA_expression_pairwise_correlation")
TP53_DIR  <- file.path(BASE_DIR, "TP53_mutation_classification")
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

DATASETS <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC")

CSN_SUBUNITS <- c("GPS1", "COPS2", "COPS3", "COPS4", "COPS5",
                   "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9")

# TP53 stratification strata
TP53_STRATA <- c("TP53_all", "TP53_MUT", "TP53_WT")

# Datasets excluded from TP53 stratification:
#   CCRCC: too few TP53 mutant samples (5/103)
#   LSCC:  too few TP53 wild-type samples (5/108)
#   OV:    too few TP53 wild-type samples (5/82)
DS_SKIP_TP53_STRATIFICATION <- c("CCRCC", "LSCC", "OV")

# Heatmap color scheme (blue-white-red diverging)
CELL_BLUE  <- "#3B4CC0"
CELL_WHITE <- "#F7F7F7"
CELL_RED   <- "#B40426"

MIN_PAIRS <- 10L   # minimum sample pairs for correlation

# ---- 2. Build Ensembl-to-gene-symbol mapping via org.Hs.eg.db -----------

message("[Mapping] Building Ensembl gene ID to symbol mapping via org.Hs.eg.db ...")

## Map gene symbols to Ensembl IDs (without version suffix)
csn_ensembl_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = CSN_SUBUNITS,
  columns = c("SYMBOL", "ENSEMBL"),
  keytype = "SYMBOL"
)
csn_ensembl_map <- csn_ensembl_map[!is.na(csn_ensembl_map$ENSEMBL), ]
csn_ensembl_map <- csn_ensembl_map[!duplicated(csn_ensembl_map$SYMBOL), ]

message(sprintf("[Mapping] Mapped %d / %d CSN subunits to Ensembl IDs",
                nrow(csn_ensembl_map), length(CSN_SUBUNITS)))

## Named vector: ENSEMBL (no version) -> SYMBOL
ensembl_to_symbol <- setNames(csn_ensembl_map$SYMBOL, csn_ensembl_map$ENSEMBL)
symbol_to_ensembl <- setNames(csn_ensembl_map$ENSEMBL, csn_ensembl_map$SYMBOL)

# ---- 3. Helper: safe filesystem name ------------------------------------

.safe_fs <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

# ---- 4. Helper: residualize expression values to covariates --------------

residualize_to_covars <- function(y, covars = NULL, min_n = 8L) {
  ## y: named numeric vector (sample -> expression)
  ## covars: data.frame with rownames = sample IDs
  ## Returns residuals after regressing out covariates
  stopifnot(!is.null(names(y)))
  sam <- names(y)
  DF  <- data.frame(row.names = sam, check.names = FALSE)

  if (!is.null(covars)) {
    C <- as.data.frame(covars, check.names = FALSE)
    rn <- rownames(C)
    if (is.null(rn)) stop("[residualize_to_covars] covars need rownames = sample IDs")
    C <- C[sam, , drop = FALSE]
    for (cn in colnames(C)) DF[[cn]] <- C[[cn]]
  }

  yv <- suppressWarnings(as.numeric(y[sam]))

  # Remove all-NA or constant columns
  if (ncol(DF) > 0) {
    all_na <- vapply(DF, function(v) all(is.na(v)), logical(1))
    if (any(all_na)) DF <- DF[, !all_na, drop = FALSE]

    is_const <- vapply(DF, function(v) {
      vv <- v[!is.na(v)]
      if (!length(vv)) {
        TRUE
      } else if (is.factor(v)) {
        nlevels(droplevels(factor(vv))) <= 1
      } else {
        stats::var(as.numeric(vv), na.rm = TRUE) == 0
      }
    }, logical(1))
    if (any(is_const)) DF <- DF[, !is_const, drop = FALSE]
  }

  ok <- is.finite(yv) & (if (ncol(DF) == 0) TRUE else stats::complete.cases(DF))
  if (sum(ok) < min_n) {
    return(setNames(rep(NA_real_, length(y)), names(y)))
  }

  des <- if (ncol(DF) == 0) {
    model.matrix(~1, data = data.frame(row.names = which(ok)))
  } else {
    model.matrix(~ 1 + ., data = DF[ok, , drop = FALSE])
  }

  fit <- lm.fit(x = des, y = yv[ok])
  out <- setNames(rep(NA_real_, length(y)), names(y))
  out[ok] <- fit$residuals
  out
}

# ---- 5. Helper: impute covariates (per covariate best-practice) ----------

impute_covariates <- function(cov_df, ds_id, stratum = "TP53_all") {
  ## cov_df: data.frame with columns sex, age, purity
  ## Returns: list(imputed_df, imputation_log)
  log_entries <- list()

  # --- Sex: mode imputation (categorical) ---
  if ("sex" %in% names(cov_df)) {
    n_miss <- sum(is.na(cov_df$sex))
    if (n_miss > 0 && n_miss < nrow(cov_df)) {
      mode_val <- names(sort(table(cov_df$sex), decreasing = TRUE))[1]
      cov_df$sex[is.na(cov_df$sex)] <- mode_val
      log_entries[["sex"]] <- sprintf(
        "Imputed %d / %d missing Sex values with mode = '%s'",
        n_miss, nrow(cov_df), mode_val
      )
    } else if (n_miss == nrow(cov_df)) {
      log_entries[["sex"]] <- "All Sex values missing; column dropped"
    } else {
      log_entries[["sex"]] <- "No missing Sex values"
    }
  }

  # --- Age: median imputation (continuous) ---
  if ("age" %in% names(cov_df)) {
    n_miss <- sum(is.na(cov_df$age))
    if (n_miss > 0 && n_miss < nrow(cov_df)) {
      med_val <- median(cov_df$age, na.rm = TRUE)
      cov_df$age[is.na(cov_df$age)] <- med_val
      log_entries[["age"]] <- sprintf(
        "Imputed %d / %d missing Age values with median = %.1f",
        n_miss, nrow(cov_df), med_val
      )
    } else if (n_miss == nrow(cov_df)) {
      log_entries[["age"]] <- "All Age values missing; column dropped"
    } else {
      log_entries[["age"]] <- "No missing Age values"
    }
  }

  # --- Tumor purity: median imputation (continuous) ---
  if ("purity" %in% names(cov_df)) {
    n_miss <- sum(is.na(cov_df$purity))
    if (n_miss > 0 && n_miss < nrow(cov_df)) {
      med_val <- median(cov_df$purity, na.rm = TRUE)
      cov_df$purity[is.na(cov_df$purity)] <- med_val
      log_entries[["purity"]] <- sprintf(
        "Imputed %d / %d missing WES_purity values with median = %.4f",
        n_miss, nrow(cov_df), med_val
      )
    } else if (n_miss == nrow(cov_df)) {
      log_entries[["purity"]] <- "All WES_purity values missing; column dropped"
    } else {
      log_entries[["purity"]] <- "No missing WES_purity values"
    }
  }

  list(imputed_df = cov_df, imputation_log = log_entries)
}

# ---- 5b. Helper: read TP53 classification and split sample IDs -----------

get_tp53_sample_lists <- function(ds_id, tp53_dir = TP53_DIR) {
  ## Returns a named list:
  ##   TP53_all = all sample IDs
  ##   TP53_MUT = TP53 mutant sample IDs (mt == 1)
  ##   TP53_WT  = TP53 wild-type sample IDs (wt == 1)
  tp53_file <- file.path(tp53_dir, paste0(ds_id, "_TP53_classification.csv"))
  if (!file.exists(tp53_file)) {
    message(sprintf("[TP53] %s: classification file not found: %s", ds_id, basename(tp53_file)))
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

# ---- 6. Main function: pairwise correlation for one dataset + stratum ----

csn_pairwise_correlation_one_ds <- function(
    ds_id,
    stratum       = "TP53_all",
    sample_subset = NULL,
    base_dir      = BASE_DIR,
    out_root      = OUT_ROOT,
    subunits      = CSN_SUBUNITS,
    min_pairs     = MIN_PAIRS,
    ensembl2sym   = ensembl_to_symbol,
    sym2ensembl   = symbol_to_ensembl
) {

  ds_dir <- file.path(base_dir, ds_id)
  ds_out <- file.path(out_root, stratum, ds_id)
  dir.create(ds_out, recursive = TRUE, showWarnings = FALSE)

  message(sprintf("\n========== Processing dataset: %s | Stratum: %s ==========",
                  ds_id, stratum))

  # --- 6a. Read mRNA expression matrix ------------------------------------
  rnaseq_file <- file.path(ds_dir,
    paste0(ds_id, "_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"))

  if (!file.exists(rnaseq_file)) {
    message(sprintf("[%s|%s] RNAseq file not found: %s",
                    ds_id, stratum, basename(rnaseq_file)))
    return(invisible(NULL))
  }

  rna_dt <- data.table::fread(rnaseq_file, header = TRUE, check.names = FALSE)
  idx_col <- names(rna_dt)[1]  # "idx"
  ensembl_ids_raw <- rna_dt[[idx_col]]

  # Strip version suffix from Ensembl IDs (e.g., ENSG00000008083.17 -> ENSG00000008083)
  ensembl_ids_base <- sub("\\.[0-9]+$", "", ensembl_ids_raw)

  # Map Ensembl IDs to gene symbols for CSN subunits
  hit_idx <- which(ensembl_ids_base %in% names(ensembl2sym))
  if (length(hit_idx) == 0) {
    message(sprintf("[%s|%s] No CSN subunit Ensembl IDs found in RNAseq data",
                    ds_id, stratum))
    return(invisible(NULL))
  }

  # Build gene-by-sample matrix for CSN subunits
  all_sample_ids <- setdiff(names(rna_dt), idx_col)
  mat0 <- as.matrix(rna_dt[hit_idx, ..all_sample_ids])
  rownames(mat0) <- ensembl2sym[ensembl_ids_base[hit_idx]]
  # Ensure numeric
  storage.mode(mat0) <- "double"

  present <- intersect(subunits, rownames(mat0))
  if (length(present) < 2) {
    message(sprintf("[%s|%s] Available CSN subunits < 2, skip", ds_id, stratum))
    return(invisible(NULL))
  }
  mat0 <- mat0[present, , drop = FALSE]

  # --- 6a2. Subset samples by stratum -------------------------------------
  if (!is.null(sample_subset)) {
    # Intersect with samples actually present in RNAseq data
    sam_use <- intersect(sample_subset, colnames(mat0))
    if (length(sam_use) < min_pairs) {
      message(sprintf("[%s|%s] Only %d samples after subsetting (need >= %d), skip",
                      ds_id, stratum, length(sam_use), min_pairs))
      return(invisible(NULL))
    }
    mat0 <- mat0[, sam_use, drop = FALSE]
  }

  sam_all <- colnames(mat0)
  message(sprintf("[%s|%s] Found %d / %d CSN subunits, %d samples: %s",
                  ds_id, stratum, length(present), length(subunits),
                  length(sam_all), paste(present, collapse = ", ")))

  # --- 6b. Read covariates ------------------------------------------------

  ## Sex & Age from meta.txt
  meta_file <- file.path(ds_dir, paste0(ds_id, "_meta.txt"))
  meta_dt   <- data.table::fread(meta_file, header = TRUE, check.names = FALSE)
  # The second row is "data_type" row; remove it
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

  ## Build covariate data.frame aligned to sam_all
  cov_df <- data.frame(row.names = sam_all, check.names = FALSE)

  # Sex: encode as numeric factor (Female=0, Male=1)
  if ("Sex" %in% names(meta_df)) {
    sex_raw <- meta_df[sam_all, "Sex"]
    sex_num <- rep(NA_real_, length(sam_all))
    sex_num[tolower(sex_raw) == "female"] <- 0
    sex_num[tolower(sex_raw) == "male"]   <- 1
    cov_df$sex <- sex_num
  } else {
    message(sprintf("[%s|%s] Warning: 'Sex' column not found in meta.txt", ds_id, stratum))
    cov_df$sex <- NA_real_
  }

  # Age: numeric
  if ("Age" %in% names(meta_df)) {
    age_raw <- meta_df[sam_all, "Age"]
    cov_df$age <- suppressWarnings(as.numeric(age_raw))
  } else {
    message(sprintf("[%s|%s] Warning: 'Age' column not found in meta.txt", ds_id, stratum))
    cov_df$age <- NA_real_
  }

  # Tumor purity (WES_purity)
  if ("WES_purity" %in% names(pheno_df)) {
    pur_raw <- pheno_df[sam_all, "WES_purity"]
    cov_df$purity <- suppressWarnings(as.numeric(pur_raw))
  } else {
    message(sprintf("[%s|%s] Warning: 'WES_purity' column not found in phenotype.txt",
                    ds_id, stratum))
    cov_df$purity <- NA_real_
  }

  ## Impute covariates
  imp_result    <- impute_covariates(cov_df, ds_id, stratum)
  cov_imputed   <- imp_result$imputed_df
  imp_log       <- imp_result$imputation_log

  # Log imputation results
  message(sprintf("[%s|%s] Covariate imputation results:", ds_id, stratum))
  for (cv_name in names(imp_log)) {
    message(sprintf("  %s: %s", cv_name, imp_log[[cv_name]]))
  }

  # --- 6c. Pairwise correlation -------------------------------------------
  pairs    <- utils::combn(present, 2, simplify = FALSE)
  all_rows <- list()

  for (p in pairs) {
    gx <- p[1]
    gy <- p[2]
    x  <- as.numeric(mat0[gx, sam_all])
    names(x) <- sam_all
    y  <- as.numeric(mat0[gy, sam_all])
    names(y) <- sam_all

    ## ---- Version 1: NoCovariate (raw Pearson) ----
    ok <- stats::complete.cases(x, y)
    if (sum(ok) >= min_pairs) {
      rr  <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "pearson"))
      fit <- stats::lm(y ~ x, data = data.frame(x = x[ok], y = y[ok]))
      all_rows[[length(all_rows) + 1L]] <- data.frame(
        dataset   = ds_id,
        stratum   = stratum,
        version   = "NoCovariate",
        gene_x    = gx,
        gene_y    = gy,
        n         = sum(ok),
        pearson_r = as.numeric(rr$estimate),
        pearson_p = rr$p.value,
        R2        = summary(fit)$r.squared,
        slope     = unname(stats::coef(fit)[["x"]]),
        intercept = unname(stats::coef(fit)[["(Intercept)"]]),
        stringsAsFactors = FALSE, check.names = FALSE
      )
    }

    ## ---- Version 2: CovariateAdj (residualized by sex, age, purity) ----
    xr <- residualize_to_covars(x, covars = cov_imputed)
    yr <- residualize_to_covars(y, covars = cov_imputed)
    ok2 <- stats::complete.cases(xr, yr)

    if (sum(ok2) >= min_pairs) {
      rr2  <- suppressWarnings(stats::cor.test(xr[ok2], yr[ok2], method = "pearson"))
      fit2 <- stats::lm(yr ~ xr, data = data.frame(xr = xr[ok2], yr = yr[ok2]))
      all_rows[[length(all_rows) + 1L]] <- data.frame(
        dataset   = ds_id,
        stratum   = stratum,
        version   = "CovariateAdj",
        gene_x    = gx,
        gene_y    = gy,
        n         = sum(ok2),
        pearson_r = as.numeric(rr2$estimate),
        pearson_p = rr2$p.value,
        R2        = summary(fit2)$r.squared,
        slope     = unname(stats::coef(fit2)[["xr"]]),
        intercept = unname(stats::coef(fit2)[["(Intercept)"]]),
        stringsAsFactors = FALSE, check.names = FALSE
      )
    }
  }  # end for pairs

  # --- 6d. Compile results and FDR correction ------------------------------
  if (!length(all_rows)) {
    message(sprintf("[%s|%s] No available pairs (insufficient samples or too many NAs)",
                    ds_id, stratum))
    return(invisible(NULL))
  }

  RES <- do.call(rbind, all_rows)

  # Benjamini-Hochberg FDR correction, per version within this stratum
  RES$pearson_padj <- ave(
    RES$pearson_p,
    interaction(RES$dataset, RES$stratum, RES$version, drop = TRUE),
    FUN = function(p) stats::p.adjust(p, method = "BH")
  )

  # --- 6e. Save CSV and XLSX ----------------------------------------------
  out_csv <- file.path(ds_out, paste0(ds_id, "_", stratum, "_pairwise_correlations.csv"))
  data.table::fwrite(RES, out_csv)
  message(sprintf("[%s|%s] CSV saved: %s", ds_id, stratum, out_csv))

  out_xlsx <- file.path(ds_out, paste0(ds_id, "_", stratum, "_pairwise_correlations.xlsx"))
  wb <- openxlsx::createWorkbook()
  for (ver in unique(RES$version)) {
    openxlsx::addWorksheet(wb, ver)
    openxlsx::writeData(wb, ver, RES[RES$version == ver, ], withFilter = TRUE)
  }
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  message(sprintf("[%s|%s] XLSX saved: %s", ds_id, stratum, out_xlsx))

  # --- 6f. Save imputation log --------------------------------------------
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
  message(sprintf("[%s|%s] Imputation log saved: %s", ds_id, stratum, imp_log_file))

  invisible(RES)
}


# ---- 7. Heatmap functions (adapted from reference script) ----------------

## 7a. Build symmetric correlation grid for heatmap
.build_corr_grid <- function(csv_path, version, sub_order = CSN_SUBUNITS) {
  stopifnot(file.exists(csv_path))
  df <- suppressMessages(readr::read_csv(csv_path, show_col_types = FALSE))
  req <- c("version", "gene_x", "gene_y", "pearson_r", "pearson_padj")
  if (!all(req %in% names(df))) {
    stop("CSV missing required columns: ", paste(setdiff(req, names(df)), collapse = ", "))
  }

  dfv <- df %>%
    dplyr::filter(.data$version == !!version) %>%
    dplyr::mutate(
      gene_x = stringr::str_trim(gene_x),
      gene_y = stringr::str_trim(gene_y)
    ) %>%
    dplyr::transmute(gene_x, gene_y, r = pearson_r, padj = pearson_padj)

  present <- intersect(sub_order, unique(c(dfv$gene_x, dfv$gene_y)))
  if (length(present) == 0L) {
    stop("No CSN subunits found for version '", version, "' in this dataset.")
  }

  # Create full symmetric grid
  grid <- tidyr::expand_grid(gene_x = present, gene_y = present) %>%
    dplyr::mutate(
      key     = paste(gene_x, gene_y, sep = "|"),
      key_rev = paste(gene_y, gene_x, sep = "|")
    )

  dfv_key <- dfv %>%
    dplyr::mutate(key = paste(gene_x, gene_y, sep = "|")) %>%
    dplyr::select(key, r, padj)

  grid2 <- grid %>%
    dplyr::left_join(dfv_key, by = "key") %>%
    dplyr::left_join(
      dplyr::rename(dfv_key, key_rev = key, r2 = r, padj2 = padj),
      by = "key_rev"
    ) %>%
    dplyr::mutate(
      r    = dplyr::coalesce(r, r2, ifelse(gene_x == gene_y, 1, NA_real_)),
      padj = dplyr::coalesce(padj, padj2, ifelse(gene_x == gene_y, NA_real_, NA_real_))
    ) %>%
    dplyr::select(gene_x, gene_y, r, padj) %>%
    dplyr::mutate(
      gene_x = factor(gene_x, levels = present),
      gene_y = factor(gene_y, levels = present),
      signif  = !is.na(padj) & padj < 0.05 & abs(r) < 0.999999
    )

  grid2
}

## 7b. Plot and save heatmap
.plot_corr_heatmap_save <- function(df_grid, title, out_base,
                                     width = 6.5, height = 6.5, dpi = 600) {
  p <- ggplot2::ggplot(df_grid, ggplot2::aes(x = gene_x, y = gene_y, fill = r)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3, na.rm = FALSE) +
    # Black dots for BH-adjusted p < 0.05
    ggplot2::geom_point(
      data = subset(df_grid, signif),
      ggplot2::aes(x = gene_x, y = gene_y),
      inherit.aes = FALSE, shape = 16, size = 2.0, color = "black", alpha = 0.9
    ) +
    ggplot2::scale_fill_gradient2(
      low      = CELL_BLUE,
      mid      = CELL_WHITE,
      high     = CELL_RED,
      midpoint = 0,
      na.value = "grey90",
      limits   = c(-1, 1),
      oob      = scales::squish,
      name     = "Pearson r"
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 0),
      axis.ticks       = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.title     = ggplot2::element_text(),
      legend.position  = "right"
    )

  # Save as TIFF, PNG, and PDF
  ggplot2::ggsave(paste0(out_base, ".tiff"), p,
                  width = width, height = height, dpi = dpi, compression = "lzw")
  ggplot2::ggsave(paste0(out_base, ".png"),  p,
                  width = width, height = height, dpi = dpi)
  ggplot2::ggsave(paste0(out_base, ".pdf"),  p,
                  width = width, height = height)

  message(sprintf("  Heatmap saved: %s (.tiff/.png/.pdf)", basename(out_base)))
  invisible(p)
}


# ---- 8. Run all datasets with TP53 stratification ------------------------

message("\n============================================================")
message("  Starting CSN subunits mRNA expression pairwise correlation")
message("  analysis with TP53 mutation status stratification")
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
      sample_sub <- NULL  # use all samples in RNAseq data
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
      csn_pairwise_correlation_one_ds(
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

# ---- 9. Combined summary across all datasets and strata ------------------

if (length(all_results) > 0) {
  combined <- do.call(rbind, all_results)
  rownames(combined) <- NULL

  # --- Per-stratum combined files ---
  for (st in TP53_STRATA) {
    st_data <- combined[combined$stratum == st, ]
    if (nrow(st_data) == 0) next

    st_dir <- file.path(OUT_ROOT, st)
    dir.create(st_dir, showWarnings = FALSE, recursive = TRUE)

    combined_csv  <- file.path(st_dir,
      paste0("ALL_datasets_", st, "_pairwise_correlations.csv"))
    combined_xlsx <- file.path(st_dir,
      paste0("ALL_datasets_", st, "_pairwise_correlations.xlsx"))

    data.table::fwrite(st_data, combined_csv)
    message(sprintf("\n[%s] Combined CSV saved: %s", st, combined_csv))

    wb_st <- openxlsx::createWorkbook()
    for (ver in unique(st_data$version)) {
      openxlsx::addWorksheet(wb_st, ver)
      openxlsx::writeData(wb_st, ver, st_data[st_data$version == ver, ], withFilter = TRUE)
    }
    openxlsx::saveWorkbook(wb_st, combined_xlsx, overwrite = TRUE)
    message(sprintf("[%s] Combined XLSX saved: %s", st, combined_xlsx))
  }

  # --- Grand combined file (all strata together) ---
  grand_csv  <- file.path(OUT_ROOT, "ALL_datasets_ALL_strata_pairwise_correlations.csv")
  grand_xlsx <- file.path(OUT_ROOT, "ALL_datasets_ALL_strata_pairwise_correlations.xlsx")

  data.table::fwrite(combined, grand_csv)
  message(sprintf("\nGrand combined CSV saved: %s", grand_csv))

  wb_all <- openxlsx::createWorkbook()
  for (st in unique(combined$stratum)) {
    for (ver in unique(combined$version)) {
      sheet_name <- paste0(st, "_", ver)
      sub_data   <- combined[combined$stratum == st & combined$version == ver, ]
      if (nrow(sub_data) > 0) {
        openxlsx::addWorksheet(wb_all, sheet_name)
        openxlsx::writeData(wb_all, sheet_name, sub_data, withFilter = TRUE)
      }
    }
  }
  openxlsx::saveWorkbook(wb_all, grand_xlsx, overwrite = TRUE)
  message(sprintf("Grand combined XLSX saved: %s", grand_xlsx))

  # --- Combined imputation logs ---
  imp_logs <- list()
  for (st in TP53_STRATA) {
    st_dir <- file.path(OUT_ROOT, st)
    if (!dir.exists(st_dir)) next
    imp_files <- list.files(st_dir, pattern = "_covariate_imputation_log\\.csv$",
                            recursive = TRUE, full.names = TRUE)
    for (f in imp_files) {
      imp_logs[[f]] <- data.table::fread(f)
    }
  }
  if (length(imp_logs) > 0) {
    combined_imp <- do.call(rbind, imp_logs)
    imp_combined_file <- file.path(OUT_ROOT,
      "ALL_datasets_ALL_strata_covariate_imputation_log.csv")
    data.table::fwrite(combined_imp, imp_combined_file)
    message(sprintf("Combined imputation log saved: %s", imp_combined_file))
  }
}


# ---- 10. Generate heatmaps from saved CSVs --------------------------------

message("\n============================================================")
message("  Generating mRNA expression correlation coefficient heatmaps")
message("  for all strata")
message("============================================================\n")

versions <- c("NoCovariate", "CovariateAdj")

for (ds in DATASETS) {

  # Determine which strata were run for this dataset
  if (ds %in% DS_SKIP_TP53_STRATIFICATION) {
    strata_to_plot <- "TP53_all"
  } else {
    strata_to_plot <- TP53_STRATA
  }

  for (st in strata_to_plot) {
    ds_out   <- file.path(OUT_ROOT, st, ds)
    csv_file <- file.path(ds_out,
      paste0(ds, "_", st, "_pairwise_correlations.csv"))

    if (!file.exists(csv_file)) {
      message(sprintf("[Heatmap] %s/%s: CSV not found, skip", ds, st))
      next
    }

    for (ver in versions) {
      grid_result <- tryCatch(
        .build_corr_grid(csv_file, ver, CSN_SUBUNITS),
        error = function(e) {
          message(sprintf("[Heatmap] %s/%s/%s: %s", ds, st, ver, conditionMessage(e)))
          NULL
        }
      )

      if (!is.null(grid_result)) {
        title_str <- sprintf("%s | CSN Subunits mRNA Expression\nPearson r (%s, %s)",
                             ds, st, ver)
        out_base  <- file.path(ds_out, sprintf("%s_%s_corrcoef_heatmap_%s",
                                                .safe_fs(ds), .safe_fs(st), .safe_fs(ver)))
        tryCatch(
          .plot_corr_heatmap_save(grid_result, title_str, out_base),
          error = function(e) {
            message(sprintf("[Heatmap] %s/%s/%s plot error: %s",
                            ds, st, ver, conditionMessage(e)))
          }
        )
      }
    }
    message(sprintf("[Heatmap] %s/%s: completed", ds, st))
  }
}

message("\n============================================================")
message("  All mRNA expression analyses completed!")
message(sprintf("  Output directory: %s", OUT_ROOT))
message("============================================================\n")
