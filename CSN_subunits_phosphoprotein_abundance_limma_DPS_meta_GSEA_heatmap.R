## =========================================================================
##  CSN Subunits Phosphoprotein Abundance limma DPS Meta-Analysis GSEA Heatmaps
##  Data source: Meta-analysis output from CSN_subunits_phosphoprotein_abundance_limma_DPS_meta_GSEA
##
##  Purpose:
##    Generate heatmaps of GSEA NES, indicating padj < 0.05
##    with dots. Combines predictors for the same stratum and collection.
##    Outputs .tiff, .pdf and the combined data as .csv.
##
##  Note:
##    By default, only pathways whose names begin with "KINASE-PSP" or
##    "PATH-NP" are plotted. You can change this by editing the
##    PATHWAY_PREFIXES parameter below.
## =========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# ---- 1. Global configuration --------------------------------------------

BASE_DIR <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb"
INPUT_DIR <- file.path(BASE_DIR, "CSN_subunits_phosphoprotein_abundance_limma_DPS_meta_GSEA")
OUT_ROOT <- file.path(BASE_DIR, "CSN_subunits_phosphoprotein_abundance_limma_DPS_meta_GSEA_heatmap")

dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# ---- Pathway prefix filter -----------------------------------------------
# Only pathways whose names start with one of these prefixes will be plotted.
# Change this vector to include/exclude different pathway categories.
# Examples of available prefixes in PTMsigDB:
#   "KINASE-PSP"  : Kinase-substrate signatures from PhosphoSitePlus
#   "PATH-NP"     : Pathway signatures from NetPath
#   "PERT-PSP"    : Perturbation signatures from PhosphoSitePlus
#   "KINASE-iKiP" : Kinase-substrate signatures from iKiP
#   "DIS-PSP"     : Disease signatures from PhosphoSitePlus
# Set to NULL to plot ALL pathways (no filtering).
PATHWAY_PREFIXES <- c("KINASE-PSP", "PATH-NP")

.pred_order_all <- c(
  "CSN_SCORE", "GPS1", "COPS2", "COPS3", "COPS4", "COPS5", "COPS6",
  "COPS7A", "COPS7B", "COPS8", "COPS9",
  "resid_GPS1", "resid_COPS2", "resid_COPS3", "resid_COPS4", 
  "resid_COPS5", "resid_COPS6", "resid_COPS7A", "resid_COPS7B", 
  "resid_COPS8", "resid_COPS9"
)
.pred_order_meta <- paste0("NES_", .pred_order_all)

# ---- 2. Find and aggregate input files ----------------------------------

message("Scanning for meta-analysis result files in: ", INPUT_DIR)
all_csvs <- list.files(INPUT_DIR, pattern = "_meta_GSEA_.*_predictor_.*\\.csv$", recursive = TRUE, full.names = TRUE)

if (length(all_csvs) == 0) {
  stop("No GSEA meta-analysis result files found in: ", INPUT_DIR)
}

message(sprintf("Found %d CSV files. Reading and aggregating...", length(all_csvs)))
dt_list <- lapply(all_csvs, function(f) {
  dt <- try(data.table::fread(f), silent = TRUE)
  if (inherits(dt, "try-error") || nrow(dt) == 0) return(NULL)
  
  # Ensure the necessary columns exist
  req_cols <- c("pathway", "collection", "csn_subunit", "NES", "padj")
  if (all(req_cols %in% names(dt))) {
    dt_sub <- dt[, ..req_cols]
    dt_sub[, stratum := basename(dirname(f))]
    return(dt_sub)
  } else {
    return(NULL)
  }
})

all_dt <- data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE)

if (nrow(all_dt) == 0) {
  stop("No valid data could be read from the CSV files.")
}

# ---- 2b. Apply pathway prefix filter ------------------------------------

if (!is.null(PATHWAY_PREFIXES) && length(PATHWAY_PREFIXES) > 0) {
  prefix_pattern <- paste0("^(", paste(PATHWAY_PREFIXES, collapse = "|"), ")")
  n_before <- length(unique(all_dt$pathway))
  all_dt <- all_dt[grepl(prefix_pattern, pathway)]
  n_after <- length(unique(all_dt$pathway))
  message(sprintf("Pathway prefix filter applied: %d -> %d unique pathways (prefixes: %s)",
                  n_before, n_after, paste(PATHWAY_PREFIXES, collapse = ", ")))
  
  if (nrow(all_dt) == 0) {
    stop("No pathways remaining after prefix filter. Check PATHWAY_PREFIXES setting.")
  }
} else {
  message("No pathway prefix filter applied. Plotting all pathways.")
}

# Create predictor column and map padj to padj_meta to match heatmap logic
all_dt[, predictor := paste0("NES_", csn_subunit)]
all_dt[, padj_meta := padj]

# Group by stratum and collection
combinations <- unique(all_dt[, .(stratum, collection)])

# ---- 3. Heatmap functions -----------------------------------------------

.meta_compute_y_order <- function(df_long, coll_tok) {
  present <- intersect(.pred_order_meta, unique(df_long$predictor))
  if (length(present) == 0) stop("No available predictors.")
  
  df_long <- df_long[predictor %in% present]
  
  if ("NES_CSN_SCORE" %in% present) {
    nes_tbl <- unique(df_long[predictor == "NES_CSN_SCORE", .(pathway, NES)])
    .top_n <- 25
    .bot_n <- 25
    keep_up <- head(nes_tbl[order(-NES)], .top_n)$pathway
    keep_dn <- head(nes_tbl[order(NES)], .bot_n)$pathway
    keep <- unique(c(keep_up, keep_dn))
    if (length(keep) > 0) {
      df_long <- df_long[pathway %in% keep]
      nes_tbl <- nes_tbl[pathway %in% keep]
    }
  }
  
  if ("NES_CSN_SCORE" %in% present) {
    ord <- unique(df_long[predictor == "NES_CSN_SCORE", .(pathway, NES)])[order(-NES)]$pathway
  } else {
    mean_NES <- df_long[, .(m = mean(NES, na.rm = TRUE)), by = pathway]
    ord <- mean_NES[order(-m)]$pathway
  }
  return(ord)
}

.meta_make_heatmap_plot <- function(df_long, coll_tok, y_order = NULL,
                                    palette = c(low = "#053061", mid = "#FFFFFF", high = "#67001F")) {
  present <- intersect(.pred_order_meta, unique(df_long$predictor))
  if (length(present) == 0) return(NULL)
  
  df_long <- df_long[predictor %in% present]
  
  if (is.null(y_order)) {
    if ("NES_CSN_SCORE" %in% present) {
      .top_n <- 25
      .bot_n <- 25
      csn_tbl <- df_long[predictor == "NES_CSN_SCORE"]
      keep_up <- head(csn_tbl[order(-NES)], .top_n)$pathway
      keep_dn <- head(csn_tbl[order(NES)], .bot_n)$pathway
      keep <- unique(c(keep_up, keep_dn))
      if (length(keep) > 0) {
        df_long <- df_long[pathway %in% keep]
      }
    }
    y_order <- .meta_compute_y_order(df_long, coll_tok)
  }
  
  df_long <- df_long[pathway %in% y_order]
  y_levels <- y_order[y_order %in% df_long$pathway]
  if (length(y_levels) == 0) return(NULL)
  
  df_long[, pathway := factor(pathway, levels = rev(y_levels))]
  
  gap <- 0.4
  needs_gap1 <- all(c("NES_CSN_SCORE", "NES_GPS1") %in% present)
  needs_gap2 <- "NES_resid_GPS1" %in% present
  
  pos_map <- numeric(length(present))
  names(pos_map) <- present
  pos <- 0
  for (p in present) {
    if (p == "NES_GPS1" && needs_gap1) pos <- pos + gap
    if (p == "NES_resid_GPS1" && needs_gap2) pos <- pos + gap
    pos <- pos + 1
    pos_map[p] <- pos
  }
  
  df_long[, xpos := pos_map[predictor]]
  
  nes_fin <- df_long$NES[is.finite(df_long$NES)]
  L <- if (length(nes_fin) > 0) max(abs(nes_fin), na.rm = TRUE) else 0
  if (L == 0) L <- 1
  
  sig_data <- df_long[is.finite(padj_meta) & padj_meta < 0.05]
  
  p <- ggplot(df_long, aes(x = xpos, y = pathway, fill = NES)) +
    geom_tile(width = 1, height = 0.9, color = NA) +
    geom_point(
      data = sig_data,
      aes(x = xpos, y = pathway),
      shape = 16, size = 1.6, color = "black", inherit.aes = FALSE
    ) +
    scale_fill_gradient2(
      low = palette[["low"]], mid = palette[["mid"]], high = palette[["high"]],
      limits = c(-L, L), midpoint = 0, oob = scales::squish, name = "NES"
    ) +
    scale_x_continuous(
      breaks = pos_map[present],
      labels = present,
      expand = expansion(mult = c(0.01, 0.01)),
      position = "top"
    ) +
    labs(x = NULL, y = NULL) +
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
  
  n_path <- nlevels(df_long$pathway)
  W <- max(8, length(present) * 0.45)
  H <- max(6, n_path * 0.22)
  
  list(plot = p, width = W, height = H, data = df_long)
}

# ---- 4. Execute plotting ------------------------------------------------

for (i in seq_len(nrow(combinations))) {
  strat <- combinations$stratum[i]
  coll <- combinations$collection[i]
  
  message(sprintf("Generating heatmap for Stratum: %s | Collection: %s", strat, coll))
  
  sub_dt <- all_dt[stratum == strat & collection == coll]
  
  # Set up palette based on stratum
  pal_cell <- c(low = "#053061", mid = "#FFFFFF", high = "#67001F") # Blue-White-Red
  if (grepl("interaction", tolower(strat))) {
    # Use Green-Purple palette for interaction models
    pal_cell <- c(low = "#2E7D32", mid = "#FFFFFF", high = "#FB8C00") 
  }
  
  # Create wide data for output CSV
  wide_nes <- dcast(sub_dt, pathway ~ predictor, value.var = "NES")
  wide_p <- dcast(sub_dt, pathway ~ predictor, value.var = "padj_meta")
  
  p_cols <- setdiff(names(wide_p), "pathway")
  setnames(wide_p, p_cols, gsub("^NES_", "padj_meta_", p_cols))
  
  wide_out <- merge(wide_nes, wide_p, by = "pathway", all = TRUE)
  
  # Sort wide_out rows by NES_CSN_SCORE if available
  if ("NES_CSN_SCORE" %in% names(wide_out)) {
    wide_out <- wide_out[order(-NES_CSN_SCORE)]
  }
  
  strat_out_dir <- file.path(OUT_ROOT, strat)
  dir.create(strat_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  base_name <- paste0(strat, "_meta_GSEA_heatmap_", coll)
  
  # Generate heatmap
  res <- try(.meta_make_heatmap_plot(sub_dt, coll_tok = coll, palette = pal_cell), silent = TRUE)
  
  if (inherits(res, "try-error") || is.null(res)) {
    message("  Plotting failed or no data for: ", base_name)
    next
  }
  
  # Filter the wide output to only include pathways that were actually plotted
  plotted_pathways <- levels(res$data$pathway)
  wide_out_plotted <- wide_out[pathway %in% plotted_pathways]
  
  # Reorder to match plot from top to bottom (y_levels order, which is reversed factor levels)
  wide_out_plotted <- wide_out_plotted[match(rev(plotted_pathways), wide_out_plotted$pathway)]
  
  # Save CSV data corresponding to the heatmap
  out_csv <- file.path(strat_out_dir, paste0(base_name, ".csv"))
  fwrite(wide_out_plotted, out_csv)
  
  # Save plots
  out_tiff <- file.path(strat_out_dir, paste0(base_name, ".tiff"))
  out_pdf <- file.path(strat_out_dir, paste0(base_name, ".pdf"))
  
  ggsave(out_tiff, res$plot, width = res$width, height = res$height, units = "in", dpi = 600, bg = "white", compression = "lzw")
  ggsave(out_pdf, res$plot, width = res$width, height = res$height, units = "in", bg = "white")
  
  message("  Saved ", base_name)
  
  # Generate ordered version for specific strata based on TP53_all
  if (strat %in% c("TP53_interaction", "TP53_MUT", "TP53_WT")) {
    all_strat <- "TP53_all"
    sub_dt_all <- all_dt[stratum == all_strat & collection == coll]
    
    if (nrow(sub_dt_all) > 0) {
      yo_full <- try(.meta_compute_y_order(sub_dt_all, coll_tok = coll), silent = TRUE)
      if (!inherits(yo_full, "try-error") && length(yo_full) > 0) {
        res_ord <- try(.meta_make_heatmap_plot(sub_dt, coll_tok = coll, y_order = yo_full, palette = pal_cell), silent = TRUE)
        
        if (!inherits(res_ord, "try-error") && !is.null(res_ord)) {
          base_name_ord <- paste0(base_name, "_ordered")
          
          plotted_pathways_ord <- levels(res_ord$data$pathway)
          wide_out_plotted_ord <- wide_out[pathway %in% plotted_pathways_ord]
          wide_out_plotted_ord <- wide_out_plotted_ord[match(rev(plotted_pathways_ord), wide_out_plotted_ord$pathway)]
          
          out_csv_ord <- file.path(strat_out_dir, paste0(base_name_ord, ".csv"))
          fwrite(wide_out_plotted_ord, out_csv_ord)
          
          out_tiff_ord <- file.path(strat_out_dir, paste0(base_name_ord, ".tiff"))
          out_pdf_ord <- file.path(strat_out_dir, paste0(base_name_ord, ".pdf"))
          
          ggsave(out_tiff_ord, res_ord$plot, width = res_ord$width, height = res_ord$height, units = "in", dpi = 600, bg = "white", compression = "lzw")
          ggsave(out_pdf_ord, res_ord$plot, width = res_ord$width, height = res_ord$height, units = "in", bg = "white")
          
          message("  Saved ", base_name_ord)
        } else {
          message("  Ordered plotting failed or no data for: ", paste0(base_name, "_ordered"))
        }
      } else {
        message("  Failed to compute y_order from TP53_all for: ", base_name)
      }
    } else {
      message("  TP53_all data not found to order: ", base_name)
    }
  }
}

message("\nAll heatmaps generated successfully in: ", OUT_ROOT)
