## =========================================================================
##  CSN Subunits Protein Abundance GSEA Per-Dataset Heatmap
##  Data source: GSEA results from CSN_subunits_protein_abundance_GSEA
##
##  Purpose:
##    Plot GSEA heatmaps per dataset, stratum, and collection based on the
##    reference CSN_CPTAC heatmap plotting code.
##    Outputs .tiff, .pdf, and .csv data files.
## =========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(scales)
})

# ---- Configuration ----
BASE_DIR <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb"
INPUT_FILE <- file.path(BASE_DIR, "CSN_subunits_protein_abundance_GSEA", "ALL_datasets_ALL_strata_GSEA_all_collections_all_CSN_subunits.csv")
OUT_ROOT <- file.path(BASE_DIR, "CSN_subunits_protein_abundance_GSEA_per_dataset_heatmap")

dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# Expected X-axis order for predictors
pred_order_all <- c(
  "CSN_SCORE", "GPS1",
  "COPS2", "COPS3", "COPS4", "COPS5", "COPS6",
  "COPS7A", "COPS7B", "COPS8", "COPS9",
  "resid_GPS1",
  "resid_COPS2", "resid_COPS3", "resid_COPS4", "resid_COPS5", "resid_COPS6",
  "resid_COPS7A", "resid_COPS7B", "resid_COPS8", "resid_COPS9"
)

# ---- Helper Functions ----

.compute_y_order <- function(df_long, coll_tok) {
  present <- intersect(pred_order_all, unique(df_long$csn_subunit))
  if (length(present) == 0) stop("No available predictors.")
  
  df_long <- df_long %>% dplyr::filter(csn_subunit %in% present)
  
  is_hallmark <- toupper(coll_tok) %in% c("HALLMARK", "H")
  if (!is_hallmark && "CSN_SCORE" %in% present) {
    csn_tbl <- df_long %>% dplyr::filter(csn_subunit == "CSN_SCORE")
    keep_up <- csn_tbl %>%
      dplyr::arrange(dplyr::desc(NES)) %>%
      dplyr::slice_head(n = 25) %>%
      dplyr::pull(pathway)
    keep_dn <- csn_tbl %>%
      dplyr::arrange(NES) %>%
      dplyr::slice_head(n = 25) %>%
      dplyr::pull(pathway)
    keep <- unique(c(keep_up, keep_dn))
    if (length(keep) > 0) {
      df_long <- df_long %>% dplyr::filter(pathway %in% keep)
    }
  }
  
  if ("CSN_SCORE" %in% present) {
    nes_csn <- df_long %>%
      dplyr::filter(csn_subunit == "CSN_SCORE") %>%
      dplyr::select(pathway, NES) %>%
      dplyr::distinct()
    ord <- nes_csn %>%
      dplyr::arrange(dplyr::desc(NES)) %>%
      dplyr::pull(pathway)
  } else {
    ord <- df_long %>%
      dplyr::group_by(pathway) %>%
      dplyr::summarise(m = mean(NES, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(m)) %>%
      dplyr::pull(pathway)
  }
  return(ord)
}

.make_heatmap_plot <- function(df_long, coll_tok, y_order = NULL,
                               palette = c(low = "#053061", mid = "#FFFFFF", high = "#67001F")) {
  present_preds <- intersect(pred_order_all, unique(df_long$csn_subunit))
  if (length(present_preds) == 0) return(NULL)
  
  df_long <- df_long %>%
    dplyr::filter(csn_subunit %in% present_preds) %>%
    dplyr::mutate(csn_subunit = factor(csn_subunit, levels = present_preds))
  
  if (is.null(y_order)) {
    is_hallmark <- toupper(coll_tok) %in% c("HALLMARK", "H")
    if (!is_hallmark && "CSN_SCORE" %in% present_preds) {
      csn_tbl <- df_long %>% dplyr::filter(csn_subunit == "CSN_SCORE")
      keep_up <- csn_tbl %>%
        dplyr::arrange(dplyr::desc(NES)) %>%
        dplyr::slice_head(n = 25) %>%
        dplyr::pull(pathway)
      keep_dn <- csn_tbl %>%
        dplyr::arrange(NES) %>%
        dplyr::slice_head(n = 25) %>%
        dplyr::pull(pathway)
      keep <- unique(c(keep_up, keep_dn))
      if (length(keep) > 0) {
        df_long <- df_long %>% dplyr::filter(pathway %in% keep)
      }
    }
    y_order <- .compute_y_order(df_long, coll_tok)
  }
  
  df_long <- df_long %>% dplyr::filter(pathway %in% y_order)
  y_levels <- y_order[y_order %in% df_long$pathway]
  if (length(y_levels) == 0) return(NULL)
  
  df_long <- df_long %>% dplyr::mutate(pathway = factor(pathway, levels = rev(y_levels)))
  
  gap <- 0.4
  needs_gap1 <- all(c("CSN_SCORE", "GPS1") %in% present_preds)
  needs_gap2 <- "resid_GPS1" %in% present_preds
  
  pos_map <- list()
  pos <- 0
  for (p in present_preds) {
    if (p == "GPS1" && needs_gap1) pos <- pos + gap
    if (p == "resid_GPS1" && needs_gap2) pos <- pos + gap
    pos <- pos + 1
    pos_map[[p]] <- pos
  }
  pos_map <- unlist(pos_map)
  df_long <- df_long %>% dplyr::mutate(xpos = unname(pos_map[as.character(csn_subunit)]))
  
  L <- max(abs(df_long$NES), na.rm = TRUE)
  if (!is.finite(L) || L == 0) L <- 1
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = xpos, y = pathway, fill = NES)) +
    ggplot2::geom_tile(width = 1, height = 0.9, color = NA) +
    ggplot2::geom_point(
      data = df_long %>% dplyr::filter(is.finite(padj) & padj < 0.05),
      ggplot2::aes(x = xpos, y = pathway),
      shape = 16, size = 1.6, color = "black", inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      low = palette[["low"]], mid = palette[["mid"]], high = palette[["high"]],
      limits = c(-L, L), oob = scales::squish, name = "NES"
    ) +
    ggplot2::scale_x_continuous(
      breaks = unname(pos_map[present_preds]),
      labels = present_preds,
      expand = ggplot2::expansion(mult = c(0.01, 0.01)),
      position = "top"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0, size = 9, margin = ggplot2::margin(b = 8)),
      axis.text.y = ggplot2::element_text(size = 9),
      legend.position = "right",
      plot.margin = ggplot2::margin(6, 12, 6, 6),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    ) +
    ggplot2::coord_cartesian(clip = "off")
  
  n_path <- nlevels(df_long$pathway)
  W <- max(8, length(present_preds) * 0.45)
  H <- max(6, n_path * 0.22)
  
  list(plot = p, width = W, height = H, data = df_long)
}

# ---- Load Data ----
if (!file.exists(INPUT_FILE)) {
  stop("Input file not found: ", INPUT_FILE, "\nPlease ensure CSN_subunits_protein_abundance_GSEA.R has been executed to generate this file.")
}
message("Loading grand GSEA result file: ", INPUT_FILE)
df_all <- data.table::fread(INPUT_FILE)

req_cols <- c("dataset", "stratum", "collection", "csn_subunit", "pathway", "NES", "padj")
missing_cols <- setdiff(req_cols, names(df_all))
if (length(missing_cols) > 0) {
  stop("Missing required columns in input file: ", paste(missing_cols, collapse = ", "))
}

combinations <- unique(df_all[, .(dataset, stratum, collection)])

# ---- Plotting Loop ----
for (i in seq_len(nrow(combinations))) {
  ds <- combinations$dataset[i]
  strat <- combinations$stratum[i]
  coll <- combinations$collection[i]
  
  message(sprintf("Processing: Dataset=%s | Stratum=%s | Collection=%s", ds, strat, coll))
  
  sub_dt <- df_all[dataset == ds & stratum == strat & collection == coll]
  
  # Palette selection
  pal_cell <- c(low = "#053061", mid = "#FFFFFF", high = "#67001F") # Blue-White-Red
  if (grepl("interaction", tolower(strat))) {
    pal_cell <- c(low = "#2E7D32", mid = "#FFFFFF", high = "#FB8C00") # Green-White-Orange
  }
  
  out_dir <- file.path(OUT_ROOT, ds)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  safe_strat <- gsub("[^A-Za-z0-9._-]+", "_", strat)
  safe_coll <- gsub("[^A-Za-z0-9._-]+", "_", coll)
  base_name <- sprintf("heatmap_%s_%s_%s", ds, safe_strat, safe_coll)
  out_prefix <- file.path(out_dir, base_name)
  
  # Generate heatmap
  res <- try(.make_heatmap_plot(sub_dt, coll_tok = coll, palette = pal_cell), silent = TRUE)
  
  if (inherits(res, "try-error") || is.null(res)) {
    message("  Plotting failed or no data for: ", base_name)
    next
  }
  
  plotted_pathways <- levels(res$data$pathway)
  
  # Save plots
  ggplot2::ggsave(paste0(out_prefix, ".tiff"), res$plot, width = res$width, height = res$height, units = "in", dpi = 600, bg = "white", compression = "lzw")
  ggplot2::ggsave(paste0(out_prefix, ".pdf"), res$plot, width = res$width, height = res$height, units = "in", bg = "white")
  
  # Save Corresponding Data as CSV
  df_out <- sub_dt %>%
    dplyr::filter(pathway %in% plotted_pathways & csn_subunit %in% pred_order_all) %>%
    dplyr::arrange(factor(pathway, levels = rev(plotted_pathways))) %>%
    dplyr::select(pathway, csn_subunit, NES, padj) %>%
    tidyr::pivot_wider(
      names_from = csn_subunit,
      values_from = c(NES, padj),
      names_glue = "{.value}_{csn_subunit}"
    )
  data.table::fwrite(df_out, paste0(out_prefix, "_data.csv"))
  message("  Saved ", base_name)
  
  # Generate ordered version for specific strata based on TP53_all
  if (strat %in% c("TP53_interaction", "TP53_MUT", "TP53_WT")) {
    all_strat <- "TP53_all"
    sub_dt_all <- df_all[dataset == ds & stratum == all_strat & collection == coll]
    
    if (nrow(sub_dt_all) > 0) {
      yo_full <- try(.compute_y_order(sub_dt_all, coll_tok = coll), silent = TRUE)
      if (!inherits(yo_full, "try-error") && length(yo_full) > 0) {
        res_ord <- try(.make_heatmap_plot(sub_dt, coll_tok = coll, y_order = yo_full, palette = pal_cell), silent = TRUE)
        
        if (!inherits(res_ord, "try-error") && !is.null(res_ord)) {
          base_name_ord <- paste0(base_name, "_ordered")
          out_prefix_ord <- file.path(out_dir, base_name_ord)
          
          plotted_pathways_ord <- levels(res_ord$data$pathway)
          
          ggplot2::ggsave(paste0(out_prefix_ord, ".tiff"), res_ord$plot, width = res_ord$width, height = res_ord$height, units = "in", dpi = 600, bg = "white", compression = "lzw")
          ggplot2::ggsave(paste0(out_prefix_ord, ".pdf"), res_ord$plot, width = res_ord$width, height = res_ord$height, units = "in", bg = "white")
          
          df_out_ord <- sub_dt %>%
            dplyr::filter(pathway %in% plotted_pathways_ord & csn_subunit %in% pred_order_all) %>%
            dplyr::arrange(factor(pathway, levels = rev(plotted_pathways_ord))) %>%
            dplyr::select(pathway, csn_subunit, NES, padj) %>%
            tidyr::pivot_wider(
              names_from = csn_subunit,
              values_from = c(NES, padj),
              names_glue = "{.value}_{csn_subunit}"
            )
          data.table::fwrite(df_out_ord, paste0(out_prefix_ord, "_data.csv"))
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

message("All per-dataset GSEA heatmaps and corresponding CSV data have been successfully generated!")
