## =========================================================================
## ===== Phosphoprotein Abundance GSEA Cross-Dataset Bubble Plot =====
## =========================================================================
## Creates publication-quality cross-dataset bubble plots for phosphoprotein
## GSEA results (PTMsigDB collection) using ggplot2.
## Input: GSEA results from CSN_subunits_phosphoprotein_abundance_GSEA.R
## Output formats: .tiff, .pdf, and .csv

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")

## -------------------------------------------------------------------------
## Configuration
## -------------------------------------------------------------------------
GSEA_PREFIX <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb/CSN_subunits_phosphoprotein_abundance_GSEA"
OUTPUT_DIR  <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb/CSN_subunits_phosphoprotein_abundance_GSEA_cross_dataset_bubble_plot"

DEFAULT_DATASETS <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC")

## Define combinations of PTMsigDB pathways
## You can configure your own combinations here.
## Each combination defines a set of PTMsigDB Hallmark collection pathways
## to display on the x-axis from left to right.
PATHWAY_COMBOS <- list(
  combo1 = c(
    "KINASE-PSP_PKACA/PRKACA",
    "KINASE-PSP_Akt1/AKT1",
    "PATH-NP_EGFR1_PATHWAY",
    "KINASE-PSP_CDK2"
  ),
  combo2 = c(
    "KINASE-PSP_Akt1/AKT1",
    "KINASE-PSP_PKACA/PRKACA",
    "PATH-NP_EGFR1_PATHWAY",
    "KINASE-PSP_AurB/AURKB",
    "KINASE-PSP_CAMK2A",
    "PATH-BI_ISCHEMIA",
    "KINASE-iKiP_AKT2",
    "KINASE-PSP_CDK2",
    "KINASE-PSP_CDK1",
    "KINASE-PSP_CDC7"
  )
)

## Select which combinations to draw
## You can specify which combinations to run by modifying this vector.
COMBOS_TO_RUN <- c("combo1", "combo2")

## Select which stratums to draw
## You can specify which stratums to run by modifying this vector.
STRATA_TO_RUN <- c("TP53_all", "TP53_WT", "TP53_MUT", "TP53_interaction")

## Select which predictors to draw
## You can specify which predictors to run by modifying this vector.
PREDICTORS_TO_RUN <- c("CSN_SCORE", "COPS7A", "COPS7B")

## Define collection parameter
DEFAULT_COLLECTION <- "PTMsigDB"

## Whether to strip PTMsigDB prefixes (KINASE-PSP_, KINASE-iKiP_, PATH-NP_,
## PATH-BI_) from x-axis labels. Set to TRUE to remove prefixes.
STRIP_PREFIX <- FALSE

## -------------------------------------------------------------------------
## Function Definition
## -------------------------------------------------------------------------
generate_gsea_bubble_plot <- function(
  gsea_prefix = GSEA_PREFIX,
  output_dir = OUTPUT_DIR,
  stratum,
  predictor,
  collection = DEFAULT_COLLECTION,
  strip_prefix = STRIP_PREFIX,
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
  cat(sprintf("===== Generating Phosphoprotein GSEA Bubble Plot =====\n"))
  cat(sprintf("  Stratum: %s | Predictor: %s\n", stratum, predictor))
  cat(sprintf("  Gene sets: %s\n", paste(genesets, collapse = ", ")))
  
  ## Collect data from all datasets
  plot_data <- data.frame()
  
  for (ds in datasets) {
    # Expected filename format: BRCA_TP53_all_GSEA_PTMsigDB_predictor_CSN_SCORE.csv
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
  
  ## Determine border thickness based on padj
  plot_data$border_thickness <- ifelse(plot_data$padj < 0.05, 1.5, 0.5)
  
  ## Create factor levels for proper ordering
  plot_data$dataset <- factor(plot_data$dataset, levels = rev(datasets)) # Reverse for y-axis
  plot_data$pathway <- factor(plot_data$pathway, levels = genesets)
  
  ## Create x-axis labels
  ## Optionally remove PTMsigDB prefixes: KINASE-PSP_, KINASE-iKiP_, PATH-NP_, PATH-BI_
  if (strip_prefix) {
    short_labels <- gsub("^KINASE-PSP_", "", genesets)
    short_labels <- gsub("^KINASE-iKiP_", "", short_labels)
    short_labels <- gsub("^PATH-NP_", "", short_labels)
    short_labels <- gsub("^PATH-BI_", "", short_labels)
  } else {
    short_labels <- genesets
  }
  
  ## Determine NES range for symmetric color scale
  nes_max <- max(abs(plot_data$NES), na.rm = TRUE)
  if (is.infinite(nes_max) || is.na(nes_max)) {
    nes_max <- 1
  }
  
  ## Create bubble plot using ggplot2
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = pathway, y = dataset)) +
    ggplot2::geom_point(ggplot2::aes(size = neg_log10_padj, fill = NES, stroke = border_thickness),
      shape = 21, color = "black"
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
    ggplot2::scale_discrete_identity(aesthetics = "stroke") +
    ggplot2::guides(stroke = "none") +
    ggplot2::scale_x_discrete(labels = short_labels) +
    ggplot2::labs(
      x = "PTMsigDB Hallmark Collection",
      y = "Dataset",
      title = sprintf("Phosphoprotein GSEA Results: %s | %s", stratum, predictor)
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
  
  cat(sprintf("===== Phosphoprotein GSEA Bubble Plot Generation Complete =====\n\n"))
  invisible(list(plot = p, tiff = tiff_file, pdf = pdf_file, csv = csv_file))
}

## -------------------------------------------------------------------------
## Execution
## -------------------------------------------------------------------------

## Loop through selected combinations, stratums, and predictors to plot
for (stratum_name in STRATA_TO_RUN) {
  for (predictor_name in PREDICTORS_TO_RUN) {
    for (combo_name in COMBOS_TO_RUN) {
      if (combo_name %in% names(PATHWAY_COMBOS)) {
        genesets <- PATHWAY_COMBOS[[combo_name]]
        
        ## Adjust width based on the number of genesets (optional tuning)
        plot_width <- max(6, length(genesets) * 0.8 + 3)
        
        generate_gsea_bubble_plot(
          stratum = stratum_name,
          predictor = predictor_name,
          genesets = genesets,
          output_prefix = paste0("GSEA_bubble_plot_", combo_name),
          width = plot_width,
          height = 6
        )
      } else {
        cat(sprintf("[Warning] Combination '%s' not found in PATHWAY_COMBOS.\n", combo_name))
      }
    }
  }
}

cat("All selected phosphoprotein GSEA bubble plots have been generated.\n")
