# CSN_subunits_phosphoprotein_abundance_DPS_meta_geneset_heatmap.R
# =========================================================================
# This script generates phosphosite set heatmaps based on DPS meta-analysis
# results from CSN_subunits_phosphoprotein_abundance_limma_DPS_meta.R.
#
# For each specified PTMsigDB phosphosite set, it:
#   1. Identifies phosphosites in the DPS meta results that belong to the set
#      (via flanking-sequence mapping to PTMsigDB gene_site IDs).
#   2. Ranks phosphosites by the CSN_SCORE predictor's meta_Z_score.
#   3. Selects the top or bottom N phosphosites (default 20).
#   4. Plots a heatmap of meta_Z_scores across all CSN subunit predictors,
#      with dots indicating meta_BH_FDR < 0.05.
#
# Output: .tiff, .pdf, and .csv data files for each heatmap.
# =========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(readxl)
})

# ---- 1. Global configuration --------------------------------------------

BASE_DIR   <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb"
INPUT_DIR  <- file.path(BASE_DIR, "CSN_subunits_phosphoprotein_abundance_limma_DPS_meta")
OUTPUT_DIR <- file.path(BASE_DIR, "CSN_subunits_phosphoprotein_abundance_DPS_meta_phosphosite_set_heatmap")

# PTMsigDB v2.0.0 phosphosite set database
PTMSIGDB_FILE <- file.path(BASE_DIR, "data_PTMsigDB_all_sites_v2.0.0.xlsx")


# ---- 2. Define target phosphosite sets -----------------------------------

# Phosphosite sets for which to plot the TOP 20 phosphosites
# (ranked by CSN_SCORE meta_Z_score in descending order)
phosphosite_sets_top20 <- c(
  "KINASE-PSP_Akt1/AKT1",
  "KINASE-PSP_PKACA/PRKACA",
  "PATH-NP_EGFR1_PATHWAY",
  "KINASE-PSP_AurB/AURKB",
  "KINASE-PSP_CAMK2A",
  "PATH-BI_ISCHEMIA",
  "KINASE-iKiP_AKT2"
)

# Phosphosite sets for which to plot the BOTTOM 20 phosphosites
# (ranked by CSN_SCORE meta_Z_score in ascending order)
phosphosite_sets_bottom20 <- c(
  "KINASE-PSP_CDK2",
  "KINASE-PSP_CDK1",
  "KINASE-PSP_CDC7"
)

# Default number of phosphosites to select
N_SITES <- 20


# ---- 3. Predictor order -------------------------------------------------

pred_order <- c(
  "Z_CSN_SCORE", "Z_GPS1",
  "Z_COPS2", "Z_COPS3", "Z_COPS4", "Z_COPS5", "Z_COPS6",
  "Z_COPS7A", "Z_COPS7B", "Z_COPS8", "Z_COPS9",
  "Z_RESIDUAL_GPS1",
  "Z_RESIDUAL_COPS2", "Z_RESIDUAL_COPS3", "Z_RESIDUAL_COPS4", "Z_RESIDUAL_COPS5", "Z_RESIDUAL_COPS6",
  "Z_RESIDUAL_COPS7A", "Z_RESIDUAL_COPS7B", "Z_RESIDUAL_COPS8", "Z_RESIDUAL_COPS9"
)


# ---- 4. Load PTMsigDB ---------------------------------------------------

message("\n============================================================")
message("  CSN Subunits Phosphoprotein Abundance DPS Meta Phosphosite Set Heatmaps")
message("============================================================\n")

message("[PTMsigDB] Loading phosphosite sets from: ", PTMSIGDB_FILE)

if (!file.exists(PTMSIGDB_FILE)) {
  stop("[PTMsigDB] File does not exist: ", PTMSIGDB_FILE)
}

ptmsigdb_raw <- read_xlsx(PTMSIGDB_FILE)
req_cols <- c("signature", "site.annotation", "site.flanking", "site.direction")
if (!all(req_cols %in% names(ptmsigdb_raw))) {
  stop("[PTMsigDB] xlsx is missing necessary fields: ",
       paste(setdiff(req_cols, names(ptmsigdb_raw)), collapse = ", "))
}

# Extract gene_site from site.annotation (e.g., "PPP1R12A_T696:15226371" -> "PPP1R12A_T696")
ptmsigdb_raw$gene_site <- toupper(sub(":.*$", "", trimws(ptmsigdb_raw$site.annotation)))

# Extract flanking sequence (15-mer) for ID mapping
ptmsigdb_raw$flanking <- toupper(trimws(ptmsigdb_raw$site.flanking))

# Build flanking-to-gene_site lookup table
valid_flank <- nzchar(ptmsigdb_raw$flanking) & nzchar(ptmsigdb_raw$gene_site) & nchar(ptmsigdb_raw$flanking) == 15
flank_lookup <- setNames(ptmsigdb_raw$gene_site[valid_flank], ptmsigdb_raw$flanking[valid_flank])
flank_lookup <- flank_lookup[!duplicated(names(flank_lookup))]

# Build phosphosite set membership lists (gene_site -> set membership)
# Each entry: signature name -> vector of unique gene_site IDs (ignoring direction for membership)
ptmsigdb_sets <- split(ptmsigdb_raw$gene_site, ptmsigdb_raw$signature)
ptmsigdb_sets <- lapply(ptmsigdb_sets, function(v) unique(v[nzchar(v)]))
ptmsigdb_sets[lengths(ptmsigdb_sets) == 0] <- NULL

message(sprintf("  PTMsigDB: %d phosphosite sets loaded", length(ptmsigdb_sets)))
message(sprintf("  Flanking lookup table: %d unique 15-mer entries", length(flank_lookup)))


# ---- 5. Helper: map DPS phosphosite IDs to PTMsigDB IDs ------------------

#' Map DPS phosphosite IDs to PTMsigDB gene_site IDs via flanking sequence
#'
#' DPS ID format: ENSG...|ENSP...|S27|GLMVCYRTDDEEDLG|1
#' This function extracts the 15-mer flanking sequence (field 4) and uses the
#' prebuilt lookup table to find the corresponding PTMsigDB gene_site ID.
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


# ---- 6. Helper: extract predictor name from file path --------------------

get_predictor_name <- function(file_path, stratum) {
  fname <- basename(file_path)
  # Remove the prefix corresponding to the stratum
  prefix <- paste0(stratum, "_meta_limma_DPS_predictor_")
  p_name <- str_remove(fname, prefix)
  p_name <- str_remove(p_name, "\\.csv$")

  # Format residual names
  if (str_starts(p_name, "resid_")) {
    p_name <- str_replace(p_name, "resid_", "RESIDUAL_")
  }

  paste0("Z_", p_name)
}


# ---- 7. Read all predictor CSVs for a specific stratum -------------------

read_stratum_data <- function(stratum) {
  stratum_dir <- file.path(INPUT_DIR, stratum)
  if (!dir.exists(stratum_dir)) {
    warning("Stratum directory does not exist: ", stratum_dir)
    return(NULL)
  }

  csv_files <- list.files(stratum_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) {
    warning("No CSV files found in: ", stratum_dir)
    return(NULL)
  }

  all_data <- lapply(csv_files, function(file) {
    pred_name <- get_predictor_name(file, stratum)
    df <- read_csv(file, show_col_types = FALSE)
    if (all(c("phospho_site_id", "meta_Z_score", "meta_BH_FDR") %in% names(df))) {
      df %>%
        select(phospho_site_id, gene_symbol, phospho_site, meta_Z_score, meta_BH_FDR) %>%
        mutate(predictor = pred_name)
    } else {
      NULL
    }
  })

  bind_rows(all_data)
}


# ---- 8. Map phosphosite IDs and create display labels --------------------

#' Create a human-readable label for a phosphosite
#' Combines gene_symbol and phospho_site (e.g., "AKT1_S473")
make_phosphosite_label <- function(gene_symbol, phospho_site) {
  paste0(gene_symbol, "_", phospho_site)
}


# ---- 9. Main function to plot phosphosite set heatmap --------------------

plot_phosphosite_set_heatmap <- function(df_long, ps_set_name, stratum,
                                         direction = "top", n_sites = N_SITES) {
  if (!(ps_set_name %in% names(ptmsigdb_sets))) {
    warning("Phosphosite set not found in PTMsigDB: ", ps_set_name)
    return(NULL)
  }

  ps_set_members <- ptmsigdb_sets[[ps_set_name]]

  # Map DPS phosphosite IDs to PTMsigDB gene_site IDs
  unique_phosphosites <- df_long %>%
    filter(predictor == "Z_CSN_SCORE") %>%
    distinct(phospho_site_id)

  ptmsigdb_ids <- map_dps_to_ptmsigdb(unique_phosphosites$phospho_site_id, flank_lookup)

  # Create mapping table: phospho_site_id -> ptmsigdb_gene_site
  id_map <- tibble(
    phospho_site_id = unique_phosphosites$phospho_site_id,
    ptmsigdb_gene_site = unname(ptmsigdb_ids)
  )

  # Filter to phosphosites that belong to this phosphosite set
  id_map_in_set <- id_map %>%
    filter(!is.na(ptmsigdb_gene_site), ptmsigdb_gene_site %in% ps_set_members)

  if (nrow(id_map_in_set) == 0) {
    message("No intersecting phosphosites found for ", ps_set_name, " in ", stratum)
    return(NULL)
  }

  # Merge mapping back to data
  df_sub <- df_long %>%
    inner_join(id_map_in_set, by = "phospho_site_id")

  if (nrow(df_sub) == 0) {
    message("No data after joining for ", ps_set_name, " in ", stratum)
    return(NULL)
  }

  # Create display label
  df_sub <- df_sub %>%
    mutate(display_label = make_phosphosite_label(gene_symbol, phospho_site))

  # Determine phosphosite ordering based on CSN_SCORE meta_Z_score
  df_csn <- df_sub %>% filter(predictor == "Z_CSN_SCORE")
  if (nrow(df_csn) == 0) {
    message("No CSN_SCORE data found to order phosphosites for ", ps_set_name, " in ", stratum)
    return(NULL)
  }

  if (direction == "top") {
    ordered_sites <- df_csn %>%
      arrange(desc(meta_Z_score)) %>%
      slice_head(n = n_sites) %>%
      pull(display_label)
  } else {
    ordered_sites <- df_csn %>%
      arrange(meta_Z_score) %>%
      slice_head(n = n_sites) %>%
      pull(display_label)
  }

  # Filter to the selected phosphosites
  df_plot <- df_sub %>% filter(display_label %in% ordered_sites)

  # Format factors for plotting
  df_plot$display_label <- factor(df_plot$display_label, levels = rev(ordered_sites))

  present_preds <- intersect(pred_order, unique(df_plot$predictor))

  # Setup x-axis positions and gaps (same layout as template)
  gap <- 0.4
  needs_gap1 <- all(c("Z_CSN_SCORE", "Z_GPS1") %in% present_preds)
  needs_gap2 <- "Z_RESIDUAL_GPS1" %in% present_preds

  pos_map <- list()
  pos <- 0
  for (pname in present_preds) {
    if (pname == "Z_GPS1" && needs_gap1) {
      pos <- pos + gap
    }
    if (pname == "Z_RESIDUAL_GPS1" && needs_gap2) {
      pos <- pos + gap
    }
    pos <- pos + 1
    pos_map[[pname]] <- pos
  }
  pos_map <- unlist(pos_map)

  df_plot$xpos <- unname(pos_map[df_plot$predictor])

  # Colors for heatmap
  palette <- c(low = "#053061", mid = "#FFFFFF", high = "#67001F")
  L <- ceiling(max(abs(df_plot$meta_Z_score[is.finite(df_plot$meta_Z_score)]), na.rm = TRUE))
  if (L == 0 || is.na(L) || is.infinite(L)) L <- 1

  # Build plot title
  dir_label <- if (direction == "top") paste0("top ", n_sites) else paste0("bottom ", n_sites)
  plot_title <- paste0(ps_set_name, " (", dir_label, " by CSN_SCORE)")

  p <- ggplot(df_plot, aes(x = xpos, y = display_label, fill = meta_Z_score)) +
    geom_tile(width = 1, height = 0.9, color = NA) +
    geom_point(
      data = filter(df_plot, is.finite(meta_BH_FDR), meta_BH_FDR < 0.05),
      aes(x = xpos, y = display_label),
      shape = 16, size = 1.6, color = "black", inherit.aes = FALSE
    ) +
    scale_fill_gradient2(
      low = palette[["low"]], mid = palette[["mid"]], high = palette[["high"]],
      limits = c(-L, L), midpoint = 0, oob = scales::squish, name = "Z Score"
    ) +
    scale_x_continuous(
      breaks = unname(pos_map[present_preds]),
      labels = present_preds, expand = expansion(mult = c(0.01, 0.01)), position = "top"
    ) +
    labs(x = NULL, y = NULL, title = plot_title) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 9, margin = margin(b = 8)),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      plot.margin = margin(6, 12, 6, 6),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    ) +
    coord_cartesian(clip = "off")

  # Prepare output data (wide format)
  df_wide <- df_plot %>%
    select(display_label, phospho_site_id, gene_symbol, phospho_site,
           ptmsigdb_gene_site, predictor, meta_Z_score, meta_BH_FDR) %>%
    pivot_wider(
      id_cols = c(display_label, phospho_site_id, gene_symbol, phospho_site, ptmsigdb_gene_site),
      names_from = predictor,
      values_from = c(meta_Z_score, meta_BH_FDR)
    ) %>%
    arrange(match(display_label, ordered_sites))

  n_path <- length(ordered_sites)
  W <- max(8, length(present_preds) * 0.45)
  H <- max(6, n_path * 0.22)

  list(plot = p, data = df_wide, width = W, height = H)
}


# ---- 10. Wrapper function to run the pipeline for selected strata --------

run_heatmap_pipeline <- function(strata = NULL) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

  # If strata is NULL, use all available strata
  if (is.null(strata)) {
    strata <- list.dirs(INPUT_DIR, full.names = FALSE, recursive = FALSE)
    strata <- strata[strata != ""]
  }

  message(sprintf("Will process %d strata: %s", length(strata), paste(strata, collapse = ", ")))

  for (stratum in strata) {
    message("\n========================================")
    message("Processing stratum: ", stratum)
    message("========================================")
    df_long <- read_stratum_data(stratum)

    if (is.null(df_long) || nrow(df_long) == 0) {
      message("Skipping stratum (no data): ", stratum)
      next
    }

    stratum_out_dir <- file.path(OUTPUT_DIR, stratum)
    dir.create(stratum_out_dir, recursive = TRUE, showWarnings = FALSE)

    # Process TOP 20 phosphosite sets
    for (ps_set in phosphosite_sets_top20) {
      message(sprintf("\n  [%s] Phosphosite set: %s (top %d)", stratum, ps_set, N_SITES))
      res <- plot_phosphosite_set_heatmap(df_long, ps_set, stratum, direction = "top", n_sites = N_SITES)
      if (!is.null(res)) {
        # Sanitize phosphosite set name for filename (replace "/" with "_")
        ps_set_safe <- gsub("/", "_", ps_set)
        base_name <- file.path(stratum_out_dir, paste0("heatmap_DPS_", stratum, "_", ps_set_safe, "_top", N_SITES))

        # Save plots
        ggsave(paste0(base_name, ".tiff"), res$plot, width = res$width, height = res$height,
               dpi = 600, bg = "white", compression = "lzw")
        ggsave(paste0(base_name, ".pdf"), res$plot, width = res$width, height = res$height,
               bg = "white")

        # Save data
        write_csv(res$data, paste0(base_name, "_data.csv"))
        message("    Saved outputs for: ", basename(base_name))
      }
    }

    # Process BOTTOM 20 phosphosite sets
    for (ps_set in phosphosite_sets_bottom20) {
      message(sprintf("\n  [%s] Phosphosite set: %s (bottom %d)", stratum, ps_set, N_SITES))
      res <- plot_phosphosite_set_heatmap(df_long, ps_set, stratum, direction = "bottom", n_sites = N_SITES)
      if (!is.null(res)) {
        # Sanitize phosphosite set name for filename (replace "/" with "_")
        ps_set_safe <- gsub("/", "_", ps_set)
        base_name <- file.path(stratum_out_dir, paste0("heatmap_DPS_", stratum, "_", ps_set_safe, "_bottom", N_SITES))

        # Save plots
        ggsave(paste0(base_name, ".tiff"), res$plot, width = res$width, height = res$height,
               dpi = 600, bg = "white", compression = "lzw")
        ggsave(paste0(base_name, ".pdf"), res$plot, width = res$width, height = res$height,
               bg = "white")

        # Save data
        write_csv(res$data, paste0(base_name, "_data.csv"))
        message("    Saved outputs for: ", basename(base_name))
      }
    }
  }

  message("\n============================================================")
  message("All heatmaps generated successfully in: ", OUTPUT_DIR)
  message("============================================================")
}


# ---- 11. Execute ---------------------------------------------------------

# Default execution: Run for all available strata.
# To run for specific strata only, pass a character vector, e.g.:
#   run_heatmap_pipeline(strata = c("TP53_all", "TP53_MUT"))
run_heatmap_pipeline()
