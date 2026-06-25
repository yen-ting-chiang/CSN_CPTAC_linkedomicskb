## =========================================================================
## ===== GSEA Bubble Plot Generation =====
## =========================================================================
## Creates publication-quality cross-dataset bubble plots for GSEA results 
## using ggplot2.
## Output formats: .tiff, .pdf, and .csv

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")

## -------------------------------------------------------------------------
## Configuration
## -------------------------------------------------------------------------
GSEA_PREFIX <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb/CSN_subunits_protein_abundance_GSEA"
OUTPUT_DIR <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb/CSN_subunits_protein_abundance_GSEA_cross_dataset_bubble_plot"

DEFAULT_DATASETS <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC")

## Define combinations of pathways
## You can configure your own combinations here.
PATHWAY_COMBOS <- list(
  combo1 = c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
    "HALLMARK_GLYCOLYSIS", 
    "HALLMARK_FATTY_ACID_METABOLISM"
  ),
  combo2 = c(
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT"
  )
)

## Select which combinations to draw
## You can specify which combination to run by modifying this vector.
COMBOS_TO_RUN <- c("combo1", "combo2")

## Define parameters
DEFAULT_STRATUM <- "TP53_all"
DEFAULT_PREDICTOR <- "CSN_SCORE"
DEFAULT_COLLECTION <- "Hallmark"

## -------------------------------------------------------------------------
## Function Definition
## -------------------------------------------------------------------------
generate_gsea_bubble_plot <- function(
  gsea_prefix = GSEA_PREFIX,
  output_dir = OUTPUT_DIR,
  stratum = DEFAULT_STRATUM,
  predictor = DEFAULT_PREDICTOR,
  collection = DEFAULT_COLLECTION,
  genesets,
  datasets = DEFAULT_DATASETS,
  color_low = "#2166AC",
  color_mid = "white",
  color_high = "#B2182B",
  output_prefix = "GSEA_bubble_plot",
  width = 8,
  height = 6,
  dpi = 300
) {
  cat(sprintf("===== Generating GSEA Bubble Plot =====\n"))
  cat(sprintf("  Stratum: %s | Predictor: %s\n", stratum, predictor))
  cat(sprintf("  Gene sets: %s\n", paste(genesets, collapse = ", ")))
  
  ## Collect data from all datasets
  plot_data <- data.frame()
  
  for (ds in datasets) {
    # Expected filename format: BRCA_TP53_all_GSEA_Hallmark_predictor_CSN_SCORE.csv
    filename <- sprintf("%s_%s_GSEA_%s_predictor_%s.csv", ds, stratum, collection, predictor)
    csv_path <- file.path(gsea_prefix, stratum, ds, filename)
    
    if (file.exists(csv_path)) {
      dt <- tryCatch(
        data.table::fread(csv_path, na.strings = c("NA", "NaN", "")),
        error = function(e) NULL
      )
      
      if (!is.null(dt) && nrow(dt) > 0) {
        ## Filter for selected genesets
        if ("pathway" %in% names(dt)) {
          dt_sub <- dt[dt$pathway %in% genesets, ]
          if (nrow(dt_sub) > 0) {
            df_sub <- data.frame(
              dataset = ds,
              pathway = dt_sub$pathway,
              NES = as.numeric(dt_sub$NES),
              padj = as.numeric(dt_sub$padj),
              stringsAsFactors = FALSE
            )
            plot_data <- rbind(plot_data, df_sub)
          }
        }
      }
    } else {
      cat(sprintf("  [Warning] File not found: %s\n", csv_path))
    }
  }
  
  if (nrow(plot_data) == 0) {
    cat(sprintf("[ERROR] No data found for bubble plot. Please check parameters.\n"))
    return(invisible(NULL))
  }
  
  ## Calculate -log10(padj) for bubble size
  plot_data$neg_log10_padj <- -log10(plot_data$padj + 1e-300) # Add small value to avoid Inf
  
  ## Create factor levels for proper ordering
  plot_data$dataset <- factor(plot_data$dataset, levels = rev(datasets)) # Reverse for y-axis
  plot_data$pathway <- factor(plot_data$pathway, levels = genesets)
  
  ## Create shortened labels for x-axis
  short_labels <- gsub("^HALLMARK_", "", genesets)
  short_labels <- gsub("_", " ", short_labels)
  
  ## Determine NES range for symmetric color scale
  nes_max <- max(abs(plot_data$NES), na.rm = TRUE)
  if (is.infinite(nes_max) || is.na(nes_max)) {
    nes_max <- 1
  }
  
  ## Create bubble plot using ggplot2
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = pathway, y = dataset)) +
    ggplot2::geom_point(ggplot2::aes(size = neg_log10_padj, fill = NES),
      shape = 21, color = "black", stroke = 0.5
    ) +
    ggplot2::scale_fill_gradient2(
      low = color_low,
      mid = color_mid,
      high = color_high,
      midpoint = 0,
      limits = c(-nes_max, nes_max),
      name = "NES"
    ) +
    ggplot2::scale_size_continuous(
      range = c(3, 12),
      name = expression(-log[10](padj))
    ) +
    ggplot2::scale_x_discrete(labels = short_labels) +
    ggplot2::labs(
      x = "Hallmark Gene Sets",
      y = "Dataset",
      title = sprintf("GSEA Results: %s | %s", stratum, predictor)
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
  
  ## Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## Generate output filename
  out_base <- file.path(output_dir, sprintf("%s_%s_%s", output_prefix, stratum, predictor))
  
  ## Save as CSV
  csv_file <- paste0(out_base, ".csv")
  write.csv(plot_data, csv_file, row.names = FALSE)
  cat(sprintf("  [Saved] CSV: %s\n", csv_file))
  
  ## Save as TIFF (publication quality, 300 DPI)
  tiff_file <- paste0(out_base, ".tiff")
  ggplot2::ggsave(
    filename = tiff_file,
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    compression = "lzw"
  )
  cat(sprintf("  [Saved] TIFF: %s\n", tiff_file))
  
  ## Save as PDF
  pdf_file <- paste0(out_base, ".pdf")
  ggplot2::ggsave(
    filename = pdf_file,
    plot = p,
    width = width,
    height = height
  )
  cat(sprintf("  [Saved] PDF: %s\n", pdf_file))
  
  cat(sprintf("===== GSEA Bubble Plot Generation Complete =====\n\n"))
  invisible(list(plot = p, tiff = tiff_file, pdf = pdf_file, csv = csv_file))
}

## -------------------------------------------------------------------------
## Execution
## -------------------------------------------------------------------------

## Loop through selected combinations to plot
for (combo_name in COMBOS_TO_RUN) {
  if (combo_name %in% names(PATHWAY_COMBOS)) {
    genesets <- PATHWAY_COMBOS[[combo_name]]
    
    ## Adjust width based on the number of genesets (optional tuning)
    plot_width <- max(6, length(genesets) * 0.8 + 3)
    
    generate_gsea_bubble_plot(
      genesets = genesets,
      output_prefix = paste0("GSEA_bubble_plot_", combo_name),
      width = plot_width,
      height = 6
    )
  } else {
    cat(sprintf("[Warning] Combination '%s' not found in PATHWAY_COMBOS.\n", combo_name))
  }
}

cat("All selected bubble plots have been generated.\n")
