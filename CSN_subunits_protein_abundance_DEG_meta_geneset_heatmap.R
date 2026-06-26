# CSN_subunits_protein_abundance_DEG_meta_geneset_heatmap.R
# This script generates gene set heatmaps based on DEG meta-analysis results
# Output: .tiff, .pdf, and .csv data files for each heatmap

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(msigdbr)
  library(tidyr)
})

# Define Paths
INPUT_DIR <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb/CSN_subunits_protein_abundance_limma_DEG_meta"
OUTPUT_DIR <- "C:/Users/danny/Documents/R_project/CSN_CPTAC_linkedomicskb/CSN_subunits_protein_abundance_DEG_meta_geneset_heatmap"

# Define Target Gene Sets
genesets_top15 <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
  "HALLMARK_GLYCOLYSIS", 
  "HALLMARK_FATTY_ACID_METABOLISM", 
  "HALLMARK_MYC_TARGETS_V1"
)

genesets_bottom15 <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
  "HALLMARK_MYC_TARGETS_V1", 
  "HALLMARK_E2F_TARGETS", 
  "HALLMARK_G2M_CHECKPOINT"
)

# Fetch Hallmark gene sets from MSigDB
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

# Predictor order matching the reference script
pred_order <- c(
  "Z_CSN_SCORE", "Z_GPS1",
  "Z_COPS2", "Z_COPS3", "Z_COPS4", "Z_COPS5", "Z_COPS6",
  "Z_COPS7A", "Z_COPS7B", "Z_COPS8", "Z_COPS9",
  "Z_RESIDUAL_GPS1",
  "Z_RESIDUAL_COPS2", "Z_RESIDUAL_COPS3", "Z_RESIDUAL_COPS4", "Z_RESIDUAL_COPS5", "Z_RESIDUAL_COPS6",
  "Z_RESIDUAL_COPS7A", "Z_RESIDUAL_COPS7B", "Z_RESIDUAL_COPS8", "Z_RESIDUAL_COPS9"
)

# Helper function to extract and format predictor name from file path
get_predictor_name <- function(file_path, stratum) {
  fname <- basename(file_path)
  # Remove the prefix corresponding to the stratum
  prefix <- paste0(stratum, "_meta_limma_DEG_predictor_")
  p_name <- str_remove(fname, prefix)
  p_name <- str_remove(p_name, "\\.csv$")
  
  # Format residual names
  if (str_starts(p_name, "resid_")) {
    p_name <- str_replace(p_name, "resid_", "RESIDUAL_")
  }
  
  paste0("Z_", p_name)
}

# Function to read all predictor CSVs for a specific stratum
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
    if ("gene_symbol" %in% names(df) && "meta_Z_score" %in% names(df) && "meta_BH_FDR" %in% names(df)) {
      df %>% 
        select(gene_symbol, meta_Z_score, meta_BH_FDR) %>%
        mutate(predictor = pred_name)
    } else {
      NULL
    }
  })
  
  bind_rows(all_data)
}

# Main function to plot heatmaps
plot_geneset_heatmap <- function(df_long, gs_name, stratum, direction = "top", n_genes = 15) {
  if (!(gs_name %in% names(hallmark_list))) {
    warning("Gene set not found in MSigDB Hallmark collection: ", gs_name)
    return(NULL)
  }
  
  gs_genes <- hallmark_list[[gs_name]]
  
  # Filter data for genes in the gene set
  df_sub <- df_long %>% filter(gene_symbol %in% gs_genes)
  
  if (nrow(df_sub) == 0) {
    message("No intersecting genes found for ", gs_name, " in ", stratum)
    return(NULL)
  }
  
  # Determine gene ordering based on CSN_SCORE meta_Z_score
  df_csn <- df_sub %>% filter(predictor == "Z_CSN_SCORE")
  if (nrow(df_csn) == 0) {
    message("No CSN_SCORE data found to order genes for ", gs_name, " in ", stratum)
    return(NULL)
  }
  
  if (direction == "top") {
    ordered_genes <- df_csn %>% 
      arrange(desc(meta_Z_score)) %>% 
      slice_head(n = n_genes) %>% 
      pull(gene_symbol)
  } else {
    ordered_genes <- df_csn %>% 
      arrange(meta_Z_score) %>% 
      slice_head(n = n_genes) %>% 
      pull(gene_symbol)
  }
  
  # Filter to the selected genes
  df_plot <- df_sub %>% filter(gene_symbol %in% ordered_genes)
  
  # Format factors for plotting
  df_plot$gene_symbol <- factor(df_plot$gene_symbol, levels = rev(ordered_genes))
  
  present_preds <- intersect(pred_order, unique(df_plot$predictor))
  
  # Setup x-axis positions and gaps similar to reference script
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
  
  p <- ggplot(df_plot, aes(x = xpos, y = gene_symbol, fill = meta_Z_score)) +
    geom_tile(width = 1, height = 0.9, color = NA) +
    geom_point(
      data = filter(df_plot, is.finite(meta_BH_FDR), meta_BH_FDR < 0.05),
      aes(x = xpos, y = gene_symbol),
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
    labs(x = NULL, y = NULL, title = paste0(gs_name, " (", direction, " 15)")) +
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
    select(gene_symbol, predictor, meta_Z_score, meta_BH_FDR) %>%
    pivot_wider(
      names_from = predictor,
      values_from = c(meta_Z_score, meta_BH_FDR)
    ) %>%
    arrange(match(gene_symbol, ordered_genes))
  
  n_path <- length(ordered_genes)
  W <- max(8, length(present_preds) * 0.45)
  H <- max(6, n_path * 0.22)
  
  list(plot = p, data = df_wide, width = W, height = H)
}

# Wrapper function to run the pipeline for selected strata
run_heatmap_pipeline <- function(strata = c("TP53_all", "TP53_mutant", "TP53_wild_type")) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  for (stratum in strata) {
    message("Processing stratum: ", stratum)
    df_long <- read_stratum_data(stratum)
    
    if (is.null(df_long)) next
    
    # Process Top 15 genesets
    for (gs in genesets_top15) {
      res <- plot_geneset_heatmap(df_long, gs, stratum, direction = "top", n_genes = 15)
      if (!is.null(res)) {
        base_name <- file.path(OUTPUT_DIR, paste0("heatmap_DEG_", stratum, "_", gs, "_top15"))
        
        # Save plots
        ggsave(paste0(base_name, ".tiff"), res$plot, width = res$width, height = res$height, dpi = 600, bg = "white", compression = "lzw")
        ggsave(paste0(base_name, ".pdf"), res$plot, width = res$width, height = res$height, bg = "white")
        
        # Save data
        write_csv(res$data, paste0(base_name, "_data.csv"))
        message("Saved outputs for: ", base_name)
      }
    }
    
    # Process Bottom 15 genesets
    for (gs in genesets_bottom15) {
      res <- plot_geneset_heatmap(df_long, gs, stratum, direction = "bottom", n_genes = 15)
      if (!is.null(res)) {
        base_name <- file.path(OUTPUT_DIR, paste0("heatmap_DEG_", stratum, "_", gs, "_bottom15"))
        
        # Save plots
        ggsave(paste0(base_name, ".tiff"), res$plot, width = res$width, height = res$height, dpi = 600, bg = "white", compression = "lzw")
        ggsave(paste0(base_name, ".pdf"), res$plot, width = res$width, height = res$height, bg = "white")
        
        # Save data
        write_csv(res$data, paste0(base_name, "_data.csv"))
        message("Saved outputs for: ", base_name)
      }
    }
  }
}

# Define available stratums based on folders in INPUT_DIR
available_strata <- list.dirs(INPUT_DIR, full.names = FALSE, recursive = FALSE)
available_strata <- available_strata[available_strata != ""]

# Default execution: Run for all available strata
# (Uncomment and modify the strata argument to run for specific subsets)
run_heatmap_pipeline(strata = available_strata)
