## =========================================================================
##  CSN Subunits Phosphoprotein Abundance Preranked PTM-SEA (Meta-Analysis)
##  Data source: Meta Limma DPS results from
##               CSN_subunits_phosphoprotein_abundance_limma_DPS_meta
##
##  Purpose:
##    For each CSN subunit and TP53 stratum, use the meta Z-score
##    from Stouffer's method as a preranked phosphosite list and perform
##    PTM-SEA (Phosphorylation-site-level Set Enrichment Analysis) against
##    PTMsigDB v2.0.0 bi-directional perturbation signatures.
##
##  Statistical rationale:
##    Using a meta Z-score (like Stouffer's Z-score) as the ranking metric for
##    preranked PTM-SEA is a robust and widely accepted approach in meta-analyses.
##    The meta Z-score inherently captures both the direction of the effect
##    and the combined statistical significance across multiple datasets,
##    making it ideal for identifying consistently enriched phospho-signatures.
##    References:
##      Subramanian A et al. (2005) PNAS 102(43):15545-15550.
##        https://doi.org/10.1073/pnas.0506580102
##      Willer, C. J. et al. (2010). Bioinformatics, 26(17), 2190-2191.
##        https://doi.org/10.1093/bioinformatics/btq340
##
##  PTM-SEA bi-directional scoring:
##    PTMsigDB perturbation signatures contain BOTH up- and down-regulated
##    phosphosites (indicated by site.direction = 'u' or 'd'). To properly
##    score concordance, we implement bi-directional enrichment:
##      - For 'u'-tagged sites: original meta_Z_score is used as ranking metric
##      - For 'd'-tagged sites: sign-flipped meta_Z_score (-Z) is used
##    This ensures that when data matches the expected perturbation pattern
##    (u-sites up AND d-sites down), the NES is strongly positive.
##    Reference: Krug et al., Mol Cell Proteomics, 2019
##               (DOI: 10.1074/mcp.TIR118.000943)
##
##  ID Mapping Strategy:
##    DPS phosphosite ID format:
##      ENSG00000067840.12|ENSP00000164640.4|T150|GLMVCYRTDDEEDLG|1
##      (ENSG | ENSP | residue+position | 15-mer flanking sequence | multiplicity)
##
##    PTMsigDB site.annotation format:
##      PPP1R12A_T696:15226371;20801872
##      (GENE_SYMBOL_RESIDUE+POSITION:PMIDs)
##
##    Matching approach: Use the 15-mer flanking sequence (field 4 of DPS ID)
##    to look up the corresponding PTMsigDB gene_site ID (from site.annotation).
## =========================================================================

# ---- 0. Load / install required packages ---------------------------------

required_cran <- c("data.table", "dplyr", "readr", "stringr", "readxl")
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
  library(readxl)
  library(fgsea)
})


# ---- 1. Global configuration --------------------------------------------

BASE_DIR  <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb"
DPS_META_ROOT  <- file.path(BASE_DIR, "CSN_subunits_phosphoprotein_abundance_limma_DPS_meta")
OUT_ROOT  <- file.path(BASE_DIR, "CSN_subunits_phosphoprotein_abundance_limma_DPS_meta_GSEA")
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


# ---- 1b. PTMsigDB configuration -----------------------------------------
#
# PTMsigDB v2.0.0 is used as the phosphosite-level gene set collection.
# The xlsx file contains bi-directional perturbation signatures with
# site.direction = 'u' (up-regulated) or 'd' (down-regulated).

PTMSIGDB_FILE <- file.path(BASE_DIR, "data_PTMsigDB_all_sites_v2.0.0.xlsx")


# ---- 1c. GSEA parameters ------------------------------------------------

GSEA_MIN_SIZE  <- 5L      # minimum gene set size (smaller for PTMsigDB)
GSEA_MAX_SIZE  <- 500L    # maximum gene set size
GSEA_NPERM     <- 10000L  # number of permutations for fgsea
GSEA_SEED      <- 42L     # random seed for reproducibility
GSEA_EPS       <- 0       # boundary for calculating p-values (0 = exact)


# ---- 2. Load PTMsigDB bi-directional gene sets --------------------------

message("\n============================================================")
message("  CSN Subunits Phosphoprotein Abundance Preranked PTM-SEA (Meta-Analysis)")
message("  Ranking metric: meta_Z_score (bi-directional)")
message("  Phosphosite set database: PTMsigDB v2.0.0")
message("============================================================\n")

message("[PTMsigDB] Loading bi-directional phosphosite sets ...")

#' Read PTMsigDB and build bi-directional pathway lists
#'
#' Each PTMsigDB perturbation signature contains phosphosites tagged with
#' direction ('u' = up-regulated, 'd' = down-regulated after perturbation).
#' To enable bi-directional scoring with fgsea, this function creates
#' direction-aware gene_site IDs:
#'   - 'u'-tagged sites: stored as "GENE_SITE" (no suffix)
#'   - 'd'-tagged sites: stored as "GENE_SITE;d"
#'
#' At runtime, the stats vector is augmented with sign-flipped entries
#' ("GENE_SITE;d" -> -Z), so fgsea correctly scores concordance.
#'
#' @param fp Path to PTMsigDB xlsx file
#' @return A list with two elements:
#'   - pathways: named list of character vectors (signature -> direction-aware IDs)
#'   - flank_lookup: named character vector (flanking_seq -> base gene_site)
read_ptmsigdb <- function(fp) {
  if (!file.exists(fp)) {
    stop("[PTMsigDB] File does not exist: ", fp)
  }

  df <- readxl::read_xlsx(fp)
  req <- c("signature", "site.annotation", "site.flanking", "site.direction")
  if (!all(req %in% names(df))) {
    stop("[PTMsigDB] xlsx is missing necessary fields: ",
         paste(setdiff(req, names(df)), collapse = ", "))
  }

  # Extract gene_site from site.annotation (e.g., "PPP1R12A_T696:15226371" -> "PPP1R12A_T696")
  gene_site <- toupper(sub(":.*$", "", trimws(df$site.annotation)))

  # Extract flanking sequence (15-mer)
  flanking <- toupper(trimws(df$site.flanking))

  # Extract site direction (u = up-regulated, d = down-regulated in perturbation)
  direction <- tolower(trimws(df$site.direction))

  # Build direction-aware gene_site IDs for bi-directional scoring:
  #   'u'-tagged sites -> "GENE_SITE"   (original Z-score will be used by fgsea)
  #   'd'-tagged sites -> "GENE_SITE;d" (sign-flipped Z-score will be used)
  dir_gene_site <- ifelse(direction == "d",
                          paste0(gene_site, ";d"),
                          gene_site)

  # Build pathway lists with direction-aware IDs
  by_sig <- split(dir_gene_site, df$signature)
  pathways <- lapply(by_sig, function(v) unique(v[nzchar(v)]))
  pathways[lengths(pathways) == 0] <- NULL

  # Build flanking-to-gene_site lookup (maps to BASE gene_site without direction)
  # Direction is signature-specific, so the flanking lookup stores only base IDs
  valid <- nzchar(flanking) & nzchar(gene_site) & nchar(flanking) == 15
  flank_lookup <- setNames(gene_site[valid], flanking[valid])
  flank_lookup <- flank_lookup[!duplicated(names(flank_lookup))]

  # Report direction statistics
  n_up <- sum(direction == "u", na.rm = TRUE)
  n_dn <- sum(direction == "d", na.rm = TRUE)
  sigs_with_u <- unique(df$signature[direction == "u"])
  sigs_with_d <- unique(df$signature[direction == "d"])
  n_bidir <- length(intersect(sigs_with_u, sigs_with_d))

  message("  [PTMsigDB] Direction breakdown:")
  message("    Total rows: ", nrow(df))
  message("    Up-tagged sites (u): ", n_up)
  message("    Down-tagged sites (d): ", n_dn)
  message("    Bi-directional signatures (contain both u and d): ",
          n_bidir, " of ", length(pathways))

  list(pathways = pathways, flank_lookup = flank_lookup)
}

ptmsigdb_data <- read_ptmsigdb(PTMSIGDB_FILE)
pathways_ptm  <- ptmsigdb_data$pathways
flank_lookup  <- ptmsigdb_data$flank_lookup

message(sprintf("  PTMsigDB: %d phosphosite sets (bi-directional)", length(pathways_ptm)))
message(sprintf("  Flanking lookup table: %d unique 15-mer entries\n", length(flank_lookup)))

## Store in a named list for consistent iteration (matching template structure)
all_gene_sets <- list(PTMsigDB = pathways_ptm)


# ---- 3. Discover available Meta DPS result files -------------------------

## Discover all per-subunit meta DPS files
dps_files <- list.files(
  DPS_META_ROOT,
  pattern    = "_meta_limma_DPS_predictor_.*\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)

if (length(dps_files) == 0) {
  stop("No Meta DPS result files found in: ", DPS_META_ROOT,
       "\nPlease run CSN_subunits_phosphoprotein_abundance_limma_DPS_meta.R first.")
}

message(sprintf("[Discovery] Found %d per-subunit Meta DPS result files", length(dps_files)))

## Parse file paths to extract stratum and subunit
parse_dps_path <- function(filepath) {
  fname <- tools::file_path_sans_ext(basename(filepath))
  # Pattern: {STRATUM}_meta_limma_DPS_predictor_{SUBUNIT}
  m <- regmatches(fname, regexec("^(.+?)_meta_limma_DPS_predictor_(.+)$", fname))[[1]]
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

file_info_list <- lapply(dps_files, parse_dps_path)
file_info_list <- file_info_list[!sapply(file_info_list, is.null)]
file_info <- do.call(rbind, file_info_list)
rownames(file_info) <- NULL

message(sprintf("[Discovery] Parsed %d valid Meta DPS files across:", nrow(file_info)))
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
  stop("No Meta DPS files remain after applying user filters. Check your configuration.")
}

message(sprintf("\n[After filters] Will process %d Meta DPS files", nrow(file_info)))


# ---- 4. Helper: map DPS phosphosites to PTMsigDB IDs --------------------

#' Map DPS phosphosite IDs to PTMsigDB gene_site IDs via flanking sequence
#'
#' DPS ID format: ENSG00000067840.12|ENSP00000164640.4|T150|GLMVCYRTDDEEDLG|1
#' This function extracts the 15-mer flanking sequence (field 4) and uses the
#' prebuilt lookup table to find the corresponding PTMsigDB gene_site ID
#' (e.g., "PPP1R12A_T696").
#'
#' @param dps_ids Character vector of DPS phosphosite IDs
#' @param lookup Named character vector (flanking -> gene_site)
#' @return Character vector of PTMsigDB gene_site IDs (NA if no match)
map_dps_to_ptmsigdb <- function(dps_ids, lookup) {
  parts <- strsplit(dps_ids, "|", fixed = TRUE)

  # Extract flanking sequence (field 4) and convert to uppercase
  flanking <- vapply(parts, function(x) {
    if (length(x) >= 4) toupper(x[4]) else NA_character_
  }, character(1))

  # Look up gene_site via flanking sequence
  mapped <- lookup[flanking]

  # Report mapping statistics
  n_total <- length(dps_ids)
  n_mapped <- sum(!is.na(mapped))
  message(sprintf("    ID mapping: %d of %d phosphosites mapped to PTMsigDB (%.1f%%)",
                  n_mapped, n_total, 100 * n_mapped / n_total))

  mapped
}


# ---- 5. Helper: run PTM-SEA for one DPS file x one collection -----------

#' Run PTM-SEA bi-directional scoring from a Meta DPS table
#'
#' Creates an augmented stats vector with sign-flipped entries for 'd'-tagged
#' sites, enabling bi-directional concordance scoring:
#'   - "GENE_SITE"   -> original meta_Z_score (for 'u'-tagged pathway members)
#'   - "GENE_SITE;d" -> -meta_Z_score (sign-flipped, for 'd'-tagged pathway members)
#'
#' Interpretation of results:
#'   - Positive NES: data is CONCORDANT with the perturbation signature
#'     (u-sites tend to be up AND d-sites tend to be down in data)
#'   - Negative NES: data is DISCORDANT with the perturbation signature
#'
#' @param dps_filepath Path to Meta DPS CSV file
#' @param stratum     Stratum label
#' @param subunit     CSN subunit name
#' @param gs_list     Named list of direction-aware phosphosite sets (from PTMsigDB)
#' @param coll_label  Collection label for output
#' @param lookup      Flanking-to-gene_site lookup table (base IDs)
#' @param min_size    Minimum gene set size
#' @param max_size    Maximum gene set size
#' @param n_perm      Number of permutations for fgsea
#' @param seed        Random seed
#' @param eps         Boundary for calculating p-values (0 = exact)
#' @return fgsea result data.frame or NULL
run_ptmsea_one <- function(dps_filepath, stratum, subunit,
                           gs_list, coll_label, lookup,
                           min_size = GSEA_MIN_SIZE,
                           max_size = GSEA_MAX_SIZE,
                           n_perm   = GSEA_NPERM,
                           seed     = GSEA_SEED,
                           eps      = GSEA_EPS) {

  ## Read Meta DPS results
  dps_df <- data.table::fread(dps_filepath, header = TRUE)

  ## Filter out the self-predictor row (the CSN subunit predicting itself)
  if ("is_predictor" %in% names(dps_df)) {
    dps_df <- dps_df[dps_df$is_predictor == FALSE, ]
  }

  ## Check required columns
  if (!all(c("phospho_site_id", "meta_Z_score") %in% names(dps_df))) {
    warning(sprintf("[%s|%s] Missing 'phospho_site_id' or 'meta_Z_score' column. Skip.",
                    stratum, subunit))
    return(NULL)
  }

  ## Remove rows with NA phospho_site_id or NA meta_Z_score
  dps_df <- dps_df[!is.na(dps_df$phospho_site_id) & !is.na(dps_df$meta_Z_score), ]

  if (nrow(dps_df) < 50) {
    warning(sprintf("[%s|%s] Only %d phosphosites after filtering. Skip.",
                    stratum, subunit, nrow(dps_df)))
    return(NULL)
  }

  ## Map DPS phosphosite IDs to PTMsigDB gene_site IDs via flanking sequence
  gene_site_ids <- map_dps_to_ptmsigdb(dps_df$phospho_site_id, lookup)

  ## Keep only successfully mapped phosphosites
  mapped_mask <- !is.na(gene_site_ids)
  if (sum(mapped_mask) < 50) {
    warning(sprintf("[%s|%s] Too few mapped phosphosites: %d. Skip.",
                    stratum, subunit, sum(mapped_mask)))
    return(NULL)
  }

  mapped_ids <- gene_site_ids[mapped_mask]
  mapped_z   <- dps_df$meta_Z_score[mapped_mask]

  ## Build named vector of meta Z-scores (base gene_site IDs)
  stats_vec <- setNames(mapped_z, mapped_ids)

  ## Handle duplicates: keep the one with the largest absolute meta_Z_score
  stats_vec <- stats_vec[order(-abs(stats_vec))]
  stats_vec <- stats_vec[!duplicated(names(stats_vec))]

  ## Remove non-finite values
  stats_vec <- stats_vec[is.finite(stats_vec)]

  if (length(stats_vec) < 50) {
    warning(sprintf("[%s|%s] Too few unique mapped phosphosites after dedup: %d. Skip.",
                    stratum, subunit, length(stats_vec)))
    return(NULL)
  }

  ## ---- PTM-SEA bi-directional augmentation ----
  ## Create sign-flipped entries for bi-directional scoring:
  ##   Original:  "GENE_SITE"   -> Z   (used by 'u'-tagged pathway members)
  ##   Flipped:   "GENE_SITE;d" -> -Z  (used by 'd'-tagged pathway members)
  ##
  ## Logic: when a 'd'-tagged site has negative Z in the data (i.e., it went
  ## down as expected by the perturbation), the flipped value (-(-Z) = +Z) is
  ## positive. This makes it rank high, contributing to positive enrichment.
  ## Conversely, if a 'd'-tagged site went up (unexpected), the flipped value
  ## is negative, contributing to negative enrichment.
  flipped_vec <- setNames(-stats_vec, paste0(names(stats_vec), ";d"))
  augmented_stats <- c(stats_vec, flipped_vec)

  ## Sort descending (required by fgsea)
  augmented_stats <- sort(augmented_stats, decreasing = TRUE)

  message(sprintf("    Augmented stats: %d base + %d flipped = %d total entries",
                  length(stats_vec), length(flipped_vec), length(augmented_stats)))

  ## Run fgsea with augmented bi-directional stats
  set.seed(seed)
  fgsea_res <- tryCatch({
    suppressWarnings(fgsea::fgseaMultilevel(
      pathways = gs_list,
      stats    = augmented_stats,
      minSize  = min_size,
      maxSize  = max_size,
      eps      = eps
    ))
  }, error = function(e) {
    message(sprintf("    [ERROR] fgsea failed: %s", conditionMessage(e)))
    NULL
  })

  if (is.null(fgsea_res) || nrow(fgsea_res) == 0) {
    return(NULL)
  }

  ## Annotate leadingEdge with direction labels for interpretability
  ## Sites ending in ";d" are down-regulated members of the perturbation signature
  fgsea_res$leadingEdge <- sapply(fgsea_res$leadingEdge, function(x) {
    labels <- ifelse(grepl(";d$", x),
                     paste0(sub(";d$", "", x), "(dn)"),
                     paste0(x, "(up)"))
    paste(labels, collapse = ";")
  })

  ## Add metadata columns
  fgsea_res$stratum         <- stratum
  fgsea_res$csn_subunit     <- subunit
  fgsea_res$collection      <- coll_label
  fgsea_res$n_ranked_sites  <- length(stats_vec)
  fgsea_res$n_augmented     <- length(augmented_stats)

  ## Reorder columns: metadata first
  meta_cols <- c("stratum", "csn_subunit", "collection", "n_ranked_sites", "n_augmented")
  result_cols <- c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge")
  col_order <- c(meta_cols, result_cols)
  col_order <- col_order[col_order %in% names(fgsea_res)]
  fgsea_res <- fgsea_res[, ..col_order]

  ## Sort by NES (descending) as requested
  fgsea_res <- fgsea_res[order(-fgsea_res$NES), ]

  return(as.data.frame(fgsea_res))
}


# ---- 6. Main loop: run PTM-SEA across all combinations ------------------

message("\n============================================================")
message("  Running preranked PTM-SEA (bi-directional)")
message(sprintf("  PTMsigDB collections: %s",
                paste(names(all_gene_sets), collapse = ", ")))
message(sprintf("  Phosphosite set size filter: [%d, %d]", GSEA_MIN_SIZE, GSEA_MAX_SIZE))
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
  dps_file <- fi$filepath

  for (coll_label in names(all_gene_sets)) {

    job_counter <- job_counter + 1L
    gs_list <- all_gene_sets[[coll_label]]

    message(sprintf("[%d/%d] %s | %s | Collection: %s",
                    job_counter, total_jobs, stratum, subunit, coll_label))

    res <- tryCatch(
      run_ptmsea_one(
        dps_filepath = dps_file,
        stratum      = stratum,
        subunit      = subunit,
        gs_list      = gs_list,
        coll_label   = coll_label,
        lookup       = flank_lookup
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
      n_concord <- sum(res$padj < 0.05 & res$NES > 0, na.rm = TRUE)
      n_discord <- sum(res$padj < 0.05 & res$NES < 0, na.rm = TRUE)
      message(sprintf("  Results: %d signatures tested | FDR < 0.05: %d (%d concordant, %d discordant) | FDR < 0.25: %d",
                      nrow(res), n_sig_005, n_concord, n_discord, n_sig_025))

      ## Save per-subunit per-collection CSV (sorted by NES descending)
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


# ---- 7. Combined and summary outputs ------------------------------------

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

      ## Sort by NES descending
      cl_data <- cl_data[order(-cl_data$NES), ]

      ## Per-stratum per-collection combined file
      combined_file <- file.path(st_dir,
        paste0(st, "_meta_GSEA_", cl, "_all_CSN_subunits.csv"))
      data.table::fwrite(cl_data, combined_file)
      message(sprintf("[%s|%s] Combined CSV: %s", st, cl, basename(combined_file)))
    }

    ## Sort by NES descending
    st_data <- st_data[order(-st_data$NES), ]

    ## Per-stratum all-collections combined file
    st_all_file <- file.path(st_dir,
      paste0(st, "_meta_GSEA_all_collections_all_CSN_subunits.csv"))
    data.table::fwrite(st_data, st_all_file)
    message(sprintf("[%s] All-collections combined CSV: %s", st, basename(st_all_file)))
  }

  ## --- Grand combined file (all strata, all collections) ---
  ## Sort by NES descending
  grand_combined <- grand_combined[order(-grand_combined$NES), ]

  grand_file <- file.path(OUT_ROOT,
    "ALL_strata_meta_GSEA_all_collections_all_CSN_subunits.csv")
  data.table::fwrite(grand_combined, grand_file)
  message(sprintf("\nGrand combined CSV: %s", grand_file))

  ## --- Summary table ---
  summary_df <- grand_combined %>%
    dplyr::group_by(stratum, csn_subunit, collection) %>%
    dplyr::summarise(
      n_ranked_sites       = dplyr::first(n_ranked_sites),
      n_augmented_entries  = dplyr::first(n_augmented),
      n_signatures_tested  = dplyr::n(),
      n_sig_FDR_0.05       = sum(padj < 0.05, na.rm = TRUE),
      n_sig_FDR_0.10       = sum(padj < 0.10, na.rm = TRUE),
      n_sig_FDR_0.25       = sum(padj < 0.25, na.rm = TRUE),
      n_concordant_FDR_0.25 = sum(padj < 0.25 & NES > 0, na.rm = TRUE),
      n_discordant_FDR_0.25 = sum(padj < 0.25 & NES < 0, na.rm = TRUE),
      top_signature        = pathway[which.min(pval)],
      top_signature_NES    = NES[which.min(pval)],
      top_signature_padj   = padj[which.min(pval)],
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
  message("\n========== Meta PTM-SEA Summary ==========")
  print(as.data.frame(summary_df), row.names = FALSE)

} else {
  message("\nNo Meta PTM-SEA results were generated.")
}


message("\n============================================================")
message("  All Meta PTM-SEA analyses completed!")
message(sprintf("  Output directory: %s", OUT_ROOT))
message("============================================================\n")
