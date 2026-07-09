## =========================================================================
##  CSN Subunits Protein Abundance vs Cell Type Correlation
##  Single-Dataset Heatmap
##  Data source: Per-dataset correlation results from
##               CSN_subunits_protein_abundance_cell_type_correlation.R
##
##  Purpose:
##    Produce per-prefix, per-stratum, per-dataset heatmaps of single-dataset
##    correlation results.
##    x-axis: predictor (CSN_SCORE, individual subunits, resid_*)
##    y-axis: cell_type (ordered top-to-bottom by CSN_SCORE's
##            pearson_r from largest to smallest)
##    Tile fill color: pearson_r
##    Significance marker: "*" where BH_FDR < 0.05
##
##  x-axis predictor ordering and gap rules follow the same convention
##  as the companion meta heatmap script
##  (CSN_subunits_protein_abundance_cell_type_correlation_meta_heatmap.R).
##
##  Output structure:
##    CSN_subunits_protein_abundance_cell_type_correlation_single_dataset_heatmap/
##      {TP53_all, TP53_MUT, TP53_WT}/
##        {CIBERSORT, ESTIMATE, xCell}/
##          {dataset}_{stratum}_{prefix}_single_dataset_correlation_heatmap.tiff
##          {dataset}_{stratum}_{prefix}_single_dataset_correlation_heatmap_data.csv
## =========================================================================

# ---- 0. Load / install required packages ---------------------------------

required_cran <- c("data.table", "dplyr", "ggplot2", "scales", "ragg")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(ragg)
})

# ---- 1. Global configuration --------------------------------------------

BASE_DIR  <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb"
INPUT_DIR <- file.path(BASE_DIR, "CSN_subunits_protein_abundance_cell_type_correlation")
OUT_ROOT  <- file.path(BASE_DIR, "CSN_subunits_protein_abundance_cell_type_correlation_single_dataset_heatmap")
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Heatmap color palette options ------------------------------------

HM_COLORSET <- "BLUE_RED_CRIMSON"

HM_PALETTES <- list(
  GSEA_DEFAULT      = c(neg = "#053061", mid = "#FFFFFF", pos = "#67001F"),
  BLUE_RED          = c(neg = "#2166AC", mid = "#F7F7F7", pos = "#B2182B"),
  BLUE_ORANGE       = c(neg = "#2B8CBE", mid = "#F7F7F7", pos = "#E34A33"),
  GREEN_MAGENTA     = c(neg = "#1B7837", mid = "#F7F7F7", pos = "#762A83"),
  BLUE_RED_DEEP     = c(neg = "#08306B", mid = "#F7F7F7", pos = "#7F0000"),
  BLUE_RED_BRIGHT   = c(neg = "#1F78B4", mid = "#FFFFFF", pos = "#E31A1C"),
  BLUE_RED_LIGHT    = c(neg = "#6BAED6", mid = "#F7F7F7", pos = "#FB6A4A"),
  BLUE_RED_TWILIGHT = c(neg = "#2C7FB8", mid = "#EEEEEE", pos = "#D7301F"),
  BLUE_RED_CRIMSON  = c(neg = "#0B3C5D", mid = "#FAFAFA", pos = "#B80C09")
)

# ---- 3. x-axis predictor ordering (same as meta heatmaps) ---------------

PRED_ORDER_ALL <- c(
  "CSN_SCORE", "GPS1",
  "COPS2", "COPS3", "COPS4", "COPS5", "COPS6",
  "COPS7A", "COPS7B", "COPS8", "COPS9",
  "resid_GPS1",
  "resid_COPS2", "resid_COPS3", "resid_COPS4", "resid_COPS5", "resid_COPS6",
  "resid_COPS7A", "resid_COPS7B", "resid_COPS8", "resid_COPS9"
)

# Output settings
OUT_FORMATS <- c("TIFF")   # Can be c("TIFF","PNG","PDF")
OUT_DPI     <- 600

# ---- 4. Scan and load single-dataset correlation result files ------------

message("Scanning for single-dataset correlation result files in: ", INPUT_DIR)
all_files <- list.files(
  INPUT_DIR,
  pattern    = "_correlation_predictor_.*\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)

if (length(all_files) == 0) {
  stop("No single-dataset correlation result files found in: ", INPUT_DIR)
}

message(sprintf("Found %d result files.", length(all_files)))

dt_list <- lapply(all_files, data.table::fread)
full_dt <- data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# Ensure required columns
req_cols <- c("dataset", "stratum", "prefix", "predictor",
              "cell_type", "pearson_r", "BH_FDR")
missing_cols <- setdiff(req_cols, names(full_dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

message(sprintf("Total records loaded: %d", nrow(full_dt)))
message(sprintf("Datasets found: %s", paste(sort(unique(full_dt$dataset)), collapse = ", ")))
message(sprintf("Strata found: %s", paste(unique(full_dt$stratum), collapse = ", ")))
message(sprintf("Prefixes found: %s", paste(unique(full_dt$prefix), collapse = ", ")))

# ---- 5. Helper: build x-axis positions with gaps -------------------------
#      Gap between CSN_SCORE and GPS1, and between COPS9 and resid_GPS1.

build_x_positions <- function(present_preds) {
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
  unlist(pos_map)
}

# ---- 6. Heatmap plotting function (single-dataset) -----------------------

plot_correlation_single_dataset_heatmap <- function(
    plot_dt,
    dataset_name,
    stratum,
    prefix,
    out_dir,
    palette     = HM_PALETTES[[HM_COLORSET]],
    pred_order  = PRED_ORDER_ALL,
    out_formats = OUT_FORMATS,
    dpi         = OUT_DPI
) {

  if (nrow(plot_dt) == 0) {
    message(sprintf("  [%s|%s|%s] No data, skip heatmap",
                    dataset_name, stratum, prefix))
    return(invisible(NULL))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # --- x-axis: predictor ordering ---
  present_preds <- intersect(pred_order, unique(plot_dt$predictor))
  if (length(present_preds) == 0) {
    message(sprintf("  [%s|%s|%s] No matching predictors, skip",
                    dataset_name, stratum, prefix))
    return(invisible(NULL))
  }
  plot_dt <- plot_dt[predictor %in% present_preds]

  # --- y-axis: cell_type ordered by CSN_SCORE's pearson_r (desc) ---
  csn_score_dt <- plot_dt[predictor == "CSN_SCORE"]
  if (nrow(csn_score_dt) > 0) {
    y_order <- csn_score_dt[order(-pearson_r)]$cell_type
  } else {
    # Fallback: use overall mean pearson_r across predictors
    mean_r <- plot_dt[, .(mean_r = mean(pearson_r, na.rm = TRUE)),
                      by = cell_type]
    y_order <- mean_r[order(-mean_r)]$cell_type
    message(sprintf("  [%s|%s|%s] CSN_SCORE not found, using mean pearson_r for y-order",
                    dataset_name, stratum, prefix))
  }

  # Keep only cell types that appear in y_order
  plot_dt <- plot_dt[cell_type %in% y_order]
  plot_dt[, cell_type := factor(cell_type, levels = rev(y_order))]

  # Build x-axis positions with gaps
  pos_map <- build_x_positions(present_preds)
  plot_dt[, xpos := pos_map[predictor]]

  # Clean predictor labels for x-axis display
  pred_labels <- present_preds
  names(pred_labels) <- present_preds

  # --- Symmetric color scale ---
  L <- max(abs(plot_dt$pearson_r), na.rm = TRUE)
  if (!is.finite(L) || L == 0) L <- 1

  # --- Build heatmap ---
  g <- ggplot2::ggplot(
    plot_dt,
    ggplot2::aes(x = xpos, y = cell_type, fill = pearson_r)
  ) +
    ggplot2::geom_tile(width = 1, height = 0.9, color = NA) +
    # Significance marker: "*" for BH_FDR < 0.05
    ggplot2::geom_text(
      data = plot_dt[is.finite(BH_FDR) & BH_FDR < 0.05],
      ggplot2::aes(x = xpos, y = cell_type, label = "*"),
      size = 4, color = "black", fontface = "bold",
      vjust = 0.75, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      low = palette[["neg"]], mid = palette[["mid"]], high = palette[["pos"]],
      limits = c(-L, L), midpoint = 0, oob = scales::squish,
      name = "Pearson r"
    ) +
    ggplot2::scale_x_continuous(
      breaks = unname(pos_map[present_preds]),
      labels = pred_labels,
      expand = ggplot2::expansion(mult = c(0.01, 0.01)),
      position = "top"
    ) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = sprintf("CSN Subunit Protein vs %s (%s)", prefix, dataset_name),
      subtitle = sprintf("Stratum: %s | Color: Pearson r | * BH FDR < 0.05",
                          stratum)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      axis.text.x.top  = ggplot2::element_text(
        angle = 45, hjust = 0, vjust = 0, size = 9,
        margin = ggplot2::margin(b = 8)
      ),
      axis.text.y      = ggplot2::element_text(size = 9),
      legend.position  = "right",
      plot.margin      = ggplot2::margin(6, 12, 6, 6),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.title       = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = 10, color = "grey40")
    ) +
    ggplot2::coord_cartesian(clip = "off")

  # --- Dynamic figure dimensions ---
  n_ct <- length(y_order)
  W <- max(8, length(present_preds) * 0.50 + 2)
  H <- max(4, n_ct * 0.28 + 1.5)

  # --- Save heatmap ---
  base_bn <- paste0(dataset_name, "_", stratum, "_", prefix,
                    "_single_dataset_correlation_heatmap")
  fmts <- toupper(out_formats)
  any_ok <- FALSE

  for (f in fmts) {
    target <- file.path(out_dir, paste0(base_bn, ".", tolower(f)))
    if (f == "TIFF") {
      tryCatch({
        ragg::agg_tiff(
          filename = target, width = W, height = H,
          units = "in", res = dpi, compression = "lzw"
        )
        invisible(print(g))
        grDevices::dev.off()
        message(sprintf("  [%s|%s|%s] Saved: %s",
                        dataset_name, stratum, prefix, basename(target)))
        any_ok <- TRUE
      }, error = function(e) {
        message(sprintf("  [%s|%s|%s] TIFF save failed: %s",
                        dataset_name, stratum, prefix, conditionMessage(e)))
        try(grDevices::dev.off(), silent = TRUE)
      })
    } else if (f %in% c("PNG", "JPG", "JPEG", "PDF")) {
      tryCatch({
        ggplot2::ggsave(
          filename = target, plot = g,
          width = W, height = H, units = "in", dpi = dpi,
          bg = "white", limitsize = FALSE
        )
        message(sprintf("  [%s|%s|%s] Saved: %s",
                        dataset_name, stratum, prefix, basename(target)))
        any_ok <- TRUE
      }, error = function(e) {
        message(sprintf("  [%s|%s|%s] %s save failed: %s",
                        dataset_name, stratum, prefix, f, conditionMessage(e)))
      })
    }
  }

  if (!any_ok) {
    message(sprintf("  [%s|%s|%s] WARNING: No output files generated",
                    dataset_name, stratum, prefix))
  }

  # --- Save corresponding data CSV ---
  csv_out <- file.path(out_dir, paste0(base_bn, "_data.csv"))
  # Export the plotted data with cell_type restored as character (for readability)
  export_dt <- copy(plot_dt)
  export_dt[, cell_type := as.character(cell_type)]
  # Re-sort: predictor order x cell_type y_order
  export_dt[, predictor := factor(predictor, levels = present_preds)]
  export_dt[, cell_type_f := factor(cell_type, levels = y_order)]
  export_dt <- export_dt[order(predictor, cell_type_f)]
  export_dt[, c("xpos", "cell_type_f") := NULL]
  data.table::fwrite(export_dt, csv_out)
  message(sprintf("  [%s|%s|%s] Data CSV saved: %s",
                  dataset_name, stratum, prefix, basename(csv_out)))

  invisible(list(plot = g, width = W, height = H))
}

# ---- 7. Main loop: iterate over strata, prefixes, and datasets -----------

message("\n============================================================")
message("  CSN Subunits Protein Abundance vs Cell Type")
message("  Correlation Single-Dataset Heatmaps")
message("  Color: Pearson r | * : BH FDR < 0.05")
message("============================================================\n")

strata   <- unique(full_dt$stratum)
prefixes <- unique(full_dt$prefix)
datasets <- sort(unique(full_dt$dataset))

for (strat in strata) {
  for (pfx in prefixes) {
    for (ds in datasets) {

      message(sprintf("\n---------- %s | %s | %s ----------", ds, strat, pfx))

      sub_dt <- full_dt[stratum == strat & prefix == pfx & dataset == ds]
      if (nrow(sub_dt) == 0) {
        message(sprintf("  [%s|%s|%s] No data found, skip", ds, strat, pfx))
        next
      }

      pfx_out_dir <- file.path(OUT_ROOT, strat, pfx)

      tryCatch(
        plot_correlation_single_dataset_heatmap(
          plot_dt      = sub_dt,
          dataset_name = ds,
          stratum      = strat,
          prefix       = pfx,
          out_dir      = pfx_out_dir
        ),
        error = function(e) {
          message(sprintf("  [%s|%s|%s] ERROR: %s",
                          ds, strat, pfx, conditionMessage(e)))
        }
      )
    }
  }
}


message("\n============================================================")
message("  All single-dataset heatmaps completed!")
message(sprintf("  Output directory: %s", OUT_ROOT))
message("============================================================\n")
