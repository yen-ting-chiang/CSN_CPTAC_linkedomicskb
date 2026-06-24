## =========================================================================
##  CSN Subunits Protein Abundance Preranked GSEA (Meta-Analysis)
##  Data source: Meta Limma DEG results from CSN_subunits_protein_abundance_limma_DEG_meta
##
##  Purpose:
##    For each CSN subunit and TP53 stratum, use the meta Z-score
##    from Stouffer's method as a preranked gene list and perform
##    Gene Set Enrichment Analysis (GSEA) against MSigDB collections.
##
##  Statistical rationale:
##    Using a meta Z-score (like Stouffer's Z-score) as the ranking metric for
##    preranked GSEA is a robust and widely accepted approach in meta-analyses.
##    The meta Z-score inherently captures both the direction of the effect
##    and the combined statistical significance across multiple datasets,
##    making it ideal for identifying consistently enriched pathways.
##    References:
##      Subramanian A et al. (2005) PNAS 102(43):15545-15550.
##        https://doi.org/10.1073/pnas.0506580102
##      Willer, C. J. et al. (2010). Bioinformatics, 26(17), 2190-2191.
##        https://doi.org/10.1093/bioinformatics/btq340
## =========================================================================

# ---- 0. Load / install required packages ---------------------------------

required_cran <- c("data.table", "dplyr", "readr", "stringr", "msigdbr")
required_bioc <- c("fgsea")

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
  library(msigdbr)
  library(fgsea)
})


# ---- 1. Global configuration --------------------------------------------

BASE_DIR  <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb"
DEG_META_ROOT  <- file.path(BASE_DIR, "CSN_subunits_protein_abundance_limma_DEG_meta")
OUT_ROOT  <- file.path(BASE_DIR, "CSN_subunits_protein_abundance_limma_DEG_meta_GSEA")
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)


# ---- 1a. USER CONFIGURABLE: Strata, subunits ----------------------------
#
# By default, all available strata and subunits are analysed.
# To run only a subset, uncomment and modify the vectors below.
# Example:
#   STRATA_TO_RUN     <- c("TP53_all")
#   SUBUNITS_TO_RUN   <- c("GPS1", "COPS5")

STRATA_TO_RUN     <- NULL   # NULL = all available strata
SUBUNITS_TO_RUN   <- NULL   # NULL = all available subunits


# ---- 1b. USER CONFIGURABLE: MSigDB gene set collections -----------------
#
# Specify which MSigDB collections to use.
# Each element is a named list with:
#   category    : MSigDB category code (required)
#   subcategory : MSigDB subcategory code (optional, use NULL for all)
#   label       : short label for output file naming (required)
#
# Default: Hallmark only.
# To add more, append to the list below.

MSIGDB_COLLECTIONS <- list(
  list(category = "H", subcategory = NULL, label = "Hallmark")
)


# ---- 1c. GSEA parameters ------------------------------------------------

GSEA_MIN_SIZE  <- 15L    # minimum gene set size
GSEA_MAX_SIZE  <- 500L   # maximum gene set size
GSEA_NPERM     <- 10000L # number of permutations for fgsea
GSEA_SEED      <- 42L    # random seed for reproducibility
GSEA_EPS       <- 0      # boundary for calculating p-values (0 = exact)


# ---- 2. Discover available Meta DEG result files -------------------------

message("\n============================================================")
message("  CSN Subunits Protein Abundance Preranked GSEA (Meta-Analysis)")
message("  Ranking metric: meta_Z_score")
message("  Gene set database: MSigDB via msigdbr")
message("============================================================\n")

## Discover all per-subunit meta DEG files
deg_files <- list.files(
  DEG_META_ROOT,
  pattern    = "_meta_limma_DEG_predictor_.*\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)

if (length(deg_files) == 0) {
  stop("No Meta DEG result files found in: ", DEG_META_ROOT,
       "\nPlease run CSN_subunits_protein_abundance_limma_DEG_meta.R first.")
}

message(sprintf("[Discovery] Found %d per-subunit Meta DEG result files", length(deg_files)))

## Parse file paths to extract stratum and subunit
parse_deg_path <- function(filepath) {
  fname <- tools::file_path_sans_ext(basename(filepath))
  # Pattern: {STRATUM}_meta_limma_DEG_predictor_{SUBUNIT}
  m <- regmatches(fname, regexec("^(.+?)_meta_limma_DEG_predictor_(.+)$", fname))[[1]]
  if (length(m) == 3) {
    return(data.frame(
      filepath = filepath,
      stratum  = m[2],
      subunit  = m[3],
      stringsAsFactors = FALSE
    ))
  }
  return(NULL)
}

file_info_list <- lapply(deg_files, parse_deg_path)
file_info_list <- file_info_list[!sapply(file_info_list, is.null)]
file_info <- do.call(rbind, file_info_list)
rownames(file_info) <- NULL

message(sprintf("[Discovery] Parsed %d valid Meta DEG files across:", nrow(file_info)))
message(sprintf("  Strata  : %s", paste(sort(unique(file_info$stratum)), collapse = ", ")))
message(sprintf("  Subunits: %s", paste(sort(unique(file_info$subunit)), collapse = ", ")))

## Apply user filters
if (!is.null(STRATA_TO_RUN)) {
  file_info <- file_info[file_info$stratum %in% STRATA_TO_RUN, ]
  message(sprintf("[Filter] Strata restricted to: %s", paste(STRATA_TO_RUN, collapse = ", ")))
}
if (!is.null(SUBUNITS_TO_RUN)) {
  file_info <- file_info[file_info$subunit %in% SUBUNITS_TO_RUN, ]
  message(sprintf("[Filter] Subunits restricted to: %s", paste(SUBUNITS_TO_RUN, collapse = ", ")))
}

if (nrow(file_info) == 0) {
  stop("No Meta DEG files remain after applying user filters. Check your configuration.")
}

message(sprintf("\n[After filters] Will process %d Meta DEG files", nrow(file_info)))


# ---- 3. Retrieve MSigDB gene sets via msigdbr ----------------------------

message("\n[MSigDB] Retrieving gene sets ...")

all_gene_sets <- list()  # named list: label -> list of gene set vectors

for (coll in MSIGDB_COLLECTIONS) {
  cat_label <- coll$label
  message(sprintf("  Fetching collection: %s (category = %s, subcategory = %s)",
                  cat_label, coll$category,
                  ifelse(is.null(coll$subcategory), "ALL", coll$subcategory)))

  if (is.null(coll$subcategory)) {
    msig_df <- msigdbr(species = "Homo sapiens", category = coll$category)
  } else {
    msig_df <- msigdbr(species = "Homo sapiens",
                       category    = coll$category,
                       subcategory = coll$subcategory)
  }

  if (nrow(msig_df) == 0) {
    warning(sprintf("No gene sets found for collection '%s'. Skipping.", cat_label))
    next
  }

  ## Build named list: gene_set_name -> character vector of gene symbols
  gs_list <- split(msig_df$gene_symbol, msig_df$gs_name)
  gs_list <- lapply(gs_list, unique)  # deduplicate within sets

  message(sprintf("  Collection '%s': %d gene sets, %d unique genes",
                  cat_label, length(gs_list), length(unique(msig_df$gene_symbol))))

  all_gene_sets[[cat_label]] <- gs_list
}

if (length(all_gene_sets) == 0) {
  stop("No MSigDB gene sets were successfully retrieved. Check MSIGDB_COLLECTIONS.")
}


# ---- 4. Helper: run fgsea for one DEG file x one collection --------------

run_fgsea_one <- function(deg_filepath, stratum, subunit,
                          gs_list, coll_label,
                          min_size = GSEA_MIN_SIZE,
                          max_size = GSEA_MAX_SIZE,
                          n_perm   = GSEA_NPERM,
                          seed     = GSEA_SEED,
                          eps      = GSEA_EPS) {

  ## Read Meta DEG results
  deg_df <- data.table::fread(deg_filepath, header = TRUE)

  ## Filter out the self-predictor row (the CSN subunit predicting itself)
  if ("is_predictor" %in% names(deg_df)) {
    deg_df <- deg_df[deg_df$is_predictor == FALSE, ]
  }

  ## Check required columns
  if (!all(c("gene_symbol", "meta_Z_score") %in% names(deg_df))) {
    warning(sprintf("[%s|%s] Missing 'gene_symbol' or 'meta_Z_score' column. Skip.",
                    stratum, subunit))
    return(NULL)
  }

  ## Remove rows with NA gene_symbol or NA meta_Z_score
  deg_df <- deg_df[!is.na(deg_df$gene_symbol) & !is.na(deg_df$meta_Z_score), ]

  ## Remove genes where gene_symbol starts with "ENSG" (unmapped Ensembl IDs)
  deg_df <- deg_df[!grepl("^ENSG", deg_df$gene_symbol), ]

  ## Handle duplicate gene symbols: keep the one with the largest absolute meta_Z_score
  deg_df <- deg_df[order(-abs(deg_df$meta_Z_score)), ]
  deg_df <- deg_df[!duplicated(deg_df$gene_symbol), ]

  if (nrow(deg_df) < 50) {
    warning(sprintf("[%s|%s] Only %d genes after filtering. Skip.",
                    stratum, subunit, nrow(deg_df)))
    return(NULL)
  }

  ## Create named vector of meta Z-scores
  ranks <- setNames(deg_df$meta_Z_score, deg_df$gene_symbol)

  ## Sort by meta_Z_score (descending) as required by fgsea
  ranks <- sort(ranks, decreasing = TRUE)

  ## Run fgsea
  set.seed(seed)
  fgsea_res <- fgsea::fgsea(
    pathways  = gs_list,
    stats     = ranks,
    minSize   = min_size,
    maxSize   = max_size,
    nPermSimple = n_perm,
    eps       = eps
  )

  if (is.null(fgsea_res) || nrow(fgsea_res) == 0) {
    return(NULL)
  }

  ## Convert leadingEdge list column to a semicolon-delimited string
  fgsea_res$leadingEdge <- sapply(fgsea_res$leadingEdge, function(x) {
    paste(x, collapse = ";")
  })

  ## Add metadata columns
  fgsea_res$stratum       <- stratum
  fgsea_res$csn_subunit   <- subunit
  fgsea_res$collection    <- coll_label
  fgsea_res$n_ranked_genes <- length(ranks)

  ## Reorder columns: metadata first
  meta_cols <- c("stratum", "csn_subunit", "collection", "n_ranked_genes")
  result_cols <- c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge")
  col_order <- c(meta_cols, result_cols)
  col_order <- col_order[col_order %in% names(fgsea_res)]
  fgsea_res <- fgsea_res[, ..col_order]

  ## Sort by NES (descending) as requested
  fgsea_res <- fgsea_res[order(-fgsea_res$NES), ]

  return(as.data.frame(fgsea_res))
}


# ---- 5. Main loop: run GSEA across all combinations ---------------------

message("\n============================================================")
message("  Running preranked GSEA")
message(sprintf("  MSigDB collections: %s",
                paste(names(all_gene_sets), collapse = ", ")))
message(sprintf("  Gene set size filter: [%d, %d]", GSEA_MIN_SIZE, GSEA_MAX_SIZE))
message(sprintf("  Permutations: %d | Seed: %d", GSEA_NPERM, GSEA_SEED))
message("============================================================\n")

all_gsea_results <- list()
result_counter <- 0L

total_jobs <- nrow(file_info) * length(all_gene_sets)
job_counter <- 0L

for (i in seq_len(nrow(file_info))) {

  fi       <- file_info[i, ]
  stratum  <- fi$stratum
  subunit  <- fi$subunit
  deg_file <- fi$filepath

  for (coll_label in names(all_gene_sets)) {

    job_counter <- job_counter + 1L
    gs_list <- all_gene_sets[[coll_label]]

    message(sprintf("[%d/%d] %s | %s | Collection: %s",
                    job_counter, total_jobs, stratum, subunit, coll_label))

    res <- tryCatch(
      run_fgsea_one(
        deg_filepath = deg_file,
        stratum      = stratum,
        subunit      = subunit,
        gs_list      = gs_list,
        coll_label   = coll_label
      ),
      error = function(e) {
        message(sprintf("  ERROR: %s", conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(res) && nrow(res) > 0) {

      result_counter <- result_counter + 1L

      ## Summary
      n_sig_005 <- sum(res$padj < 0.05, na.rm = TRUE)
      n_sig_025 <- sum(res$padj < 0.25, na.rm = TRUE)
      message(sprintf("  Results: %d pathways tested | FDR < 0.05: %d | FDR < 0.25: %d",
                      nrow(res), n_sig_005, n_sig_025))

      ## Save per-subunit per-collection CSV
      out_dir <- file.path(OUT_ROOT, stratum)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      out_file <- file.path(out_dir,
        paste0(stratum, "_meta_GSEA_", coll_label,
               "_predictor_", subunit, ".csv"))
      data.table::fwrite(res, out_file)
      message(sprintf("  Saved: %s", basename(out_file)))

      ## Accumulate
      key <- paste(stratum, subunit, coll_label, sep = "|")
      all_gsea_results[[key]] <- res

    } else {
      message("  No results (skipped or empty)")
    }
  }
}


# ---- 6. Combined and summary outputs ------------------------------------

if (length(all_gsea_results) > 0) {

  grand_combined <- do.call(rbind, all_gsea_results)
  rownames(grand_combined) <- NULL

  ## --- Per-collection combined files within each stratum ---
  strata_present <- unique(grand_combined$stratum)
  colls_present  <- unique(grand_combined$collection)

  for (st in strata_present) {
    st_data <- grand_combined[grand_combined$stratum == st, ]
    if (nrow(st_data) == 0) next

    st_dir <- file.path(OUT_ROOT, st)
    dir.create(st_dir, showWarnings = FALSE, recursive = TRUE)

    for (cl in colls_present) {
      cl_data <- st_data[st_data$collection == cl, ]
      if (nrow(cl_data) == 0) next

      ## Per-stratum per-collection combined file
      combined_file <- file.path(st_dir,
        paste0(st, "_meta_GSEA_", cl, "_all_CSN_subunits.csv"))
      data.table::fwrite(cl_data, combined_file)
      message(sprintf("[%s|%s] Combined CSV: %s", st, cl, basename(combined_file)))
    }

    ## Per-stratum all-collections combined file
    st_all_file <- file.path(st_dir,
      paste0(st, "_meta_GSEA_all_collections_all_CSN_subunits.csv"))
    data.table::fwrite(st_data, st_all_file)
    message(sprintf("[%s] All-collections combined CSV: %s", st, basename(st_all_file)))
  }

  ## --- Grand combined file (all strata, all collections) ---
  grand_file <- file.path(OUT_ROOT,
    "ALL_strata_meta_GSEA_all_collections_all_CSN_subunits.csv")
  data.table::fwrite(grand_combined, grand_file)
  message(sprintf("\nGrand combined CSV: %s", grand_file))

  ## --- Summary table ---
  summary_df <- grand_combined %>%
    dplyr::group_by(stratum, csn_subunit, collection) %>%
    dplyr::summarise(
      n_ranked_genes    = dplyr::first(n_ranked_genes),
      n_pathways_tested = dplyr::n(),
      n_sig_FDR_0.05    = sum(padj < 0.05, na.rm = TRUE),
      n_sig_FDR_0.10    = sum(padj < 0.10, na.rm = TRUE),
      n_sig_FDR_0.25    = sum(padj < 0.25, na.rm = TRUE),
      n_up_FDR_0.25     = sum(padj < 0.25 & NES > 0, na.rm = TRUE),
      n_down_FDR_0.25   = sum(padj < 0.25 & NES < 0, na.rm = TRUE),
      top_pathway       = pathway[which.min(pval)],
      top_pathway_NES   = NES[which.min(pval)],
      top_pathway_padj  = padj[which.min(pval)],
      .groups = "drop"
    )

  summary_file <- file.path(OUT_ROOT,
    "ALL_strata_meta_GSEA_summary.csv")
  data.table::fwrite(summary_df, summary_file)
  message(sprintf("Summary CSV: %s", summary_file))

  ## Per-stratum summary files
  for (st in strata_present) {
    st_summary <- summary_df[summary_df$stratum == st, ]
    if (nrow(st_summary) == 0) next
    st_dir <- file.path(OUT_ROOT, st)
    st_summary_file <- file.path(st_dir,
      paste0(st, "_meta_GSEA_summary.csv"))
    data.table::fwrite(st_summary, st_summary_file)
    message(sprintf("[%s] Summary CSV: %s", st, basename(st_summary_file)))
  }

  ## Print summary table to console
  message("\n========== Meta GSEA Summary ==========")
  print(as.data.frame(summary_df), row.names = FALSE)

} else {
  message("\nNo Meta GSEA results were generated.")
}


message("\n============================================================")
message("  All Meta GSEA analyses completed!")
message(sprintf("  Output directory: %s", OUT_ROOT))
message("============================================================\n")
