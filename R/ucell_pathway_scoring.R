# Load required packages
library(Seurat)
library(ggplot2)
library(UCell)
library(viridis)
library(ggpubr) # dependency of stat_compare_means

RunPathwayAnalysis <- function(seu_object, pathway_set, GSEA_set, cell_color, 
                               group_col = "Celltype", 
                               cell_levels = NULL, 
                               my_comparisons = NULL) {
  
  cat(paste0(">>> Processing pathway: ", pathway_set, "\n"))
  
  # 1. Extract the gene list of the specified metabolic pathway
  pathway <- subset(GSEA_set, term == pathway_set)
  if(nrow(pathway) == 0) stop("Pathway not found in GSEA_set, please check the spelling!")
  
  pathway_list <- list(pathway$gene)
  names(pathway_list)[1] <- pathway_set
  
  # 2. Compute UCell scores (AddModuleScore_UCell automatically returns a Seurat object containing the scores)
  seu_object <- AddModuleScore_UCell(seu_object, features = pathway_list, name = "_score")
  score_name <- paste0(pathway_set, "_score")
  
  # 3. Visualize and save the UMAP
  umap_pdf <- paste0(pathway_set, "_UMAP.pdf")
  pdf(file = umap_pdf, width = 7, height = 5.5)
  # Note: within a function the plot must be wrapped in print() to be written to the PDF successfully
  p_umap <- FeaturePlot(seu_object, features = score_name, order = TRUE, cols = plasma(256))
  print(p_umap)
  dev.off()
  cat(paste0("UMAP plot saved to: ", umap_pdf, "\n"))
  
  # 4. Summarize and print the cell type distribution
  cat("Current cell distribution:\n")
  print(table(seu_object[[group_col]]))
  
  # 5. Prepare violin plot data
  df <- FetchData(seu_object, vars = c(group_col, score_name))
  colnames(df) <- c("Celltype", "score")
  
  # Dynamically set the cell factor order (if cell_levels is provided)
  if (!is.null(cell_levels)) {
    df$Celltype <- factor(df$Celltype, levels = cell_levels)
    # Remove cells not in levels (e.g. when focusing on a few subpopulations)
    df <- df[!is.na(df$Celltype), ] 
  } else {
    df$Celltype <- as.factor(df$Celltype)
  }
  
  # 6. Draw the violin plot
  vln_pdf <- paste0(pathway_set, "_violin.pdf")
  pdf(file = vln_pdf, width = 5, height = 5)
  
  p_vln <- ggplot(df, aes(x = Celltype, y = score, fill = Celltype)) +
    geom_violin(color = "black", linewidth = 1) + # ggplot2 recommends linewidth instead of size
    theme_classic() +
    theme(text = element_text(size = 10, colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 15),
          axis.text.x = element_text(colour = "black", size = 12),
          axis.text.y = element_text(colour = "black", size = 10),
          axis.title.y = element_text(color = "black", size = 12),
          axis.line = element_line(linewidth = 1),
          legend.position = "none") +
    scale_fill_manual(values = cell_color) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA)
  
  # Dynamically add statistical tests (if my_comparisons is provided)
  if (!is.null(my_comparisons)) {
    p_vln <- p_vln + stat_compare_means(method = "t.test", hide.ns = FALSE,
                                        comparisons = my_comparisons,
                                        label = "p.signif", bracket.size = 0.8, size = 6)
  }
  
  print(p_vln)
  dev.off()
  cat(paste0("Violin plot saved to: ", vln_pdf, "\n\n"))
  
  # Return the new Seurat object with UCell scores for subsequent use
  return(seu_object)
}



RunPathwayAnalysis_KEGG <- function(seu_object, pathway_name, gene_set_list, cell_color, 
                                    group_col = "Celltype", 
                                    cell_levels = NULL, 
                                    my_comparisons = NULL) {
  
  cat(paste0(">>> Processing pathway: ", pathway_name, "\n"))
  
  # 1. Minimal gene extraction logic (adapted to the msigdbr List format)
  if (!(pathway_name %in% names(gene_set_list))) {
    stop("Pathway not found in the provided list, please check the spelling (use names(gene_set_list) to inspect)!")
  }
  
  # Extract the corresponding gene set directly from the list, keeping the single-element list format for UCell
  pathway_list <- gene_set_list[pathway_name] 
  
  # 2. Compute UCell scores
  seu_object <- AddModuleScore_UCell(seu_object, features = pathway_list, name = "_score")
  score_name <- paste0(pathway_name, "_score")
  
  # 3. Visualize and save the UMAP
  umap_pdf <- paste0(pathway_name, "_UMAP.pdf")
  pdf(file = umap_pdf, width = 7, height = 5.5)
  p_umap <- FeaturePlot(seu_object, features = score_name, order = TRUE, cols = plasma(256))
  print(p_umap)
  dev.off()
  cat(paste0("UMAP plot saved to: ", umap_pdf, "\n"))
  
  # 4. Summarize the cell distribution
  cat("Current cell distribution:\n")
  print(table(seu_object[[group_col]]))
  
  # 5. Prepare violin plot data
  df <- FetchData(seu_object, vars = c(group_col, score_name))
  colnames(df) <- c("Celltype", "score")
  
  if (!is.null(cell_levels)) {
    df$Celltype <- factor(df$Celltype, levels = cell_levels)
    df <- df[!is.na(df$Celltype), ] 
  } else {
    df$Celltype <- as.factor(df$Celltype)
  }
  
  # 6. Draw the violin plot
  vln_pdf <- paste0(pathway_name, "_violin.pdf")
  pdf(file = vln_pdf, width = 5, height = 5)
  
  p_vln <- ggplot(df, aes(x = Celltype, y = score, fill = Celltype)) +
    geom_violin(color = "black", linewidth = 1) + 
    theme_classic() +
    theme(text = element_text(size = 10, colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 15),
          axis.text.x = element_text(colour = "black", size = 12),
          axis.text.y = element_text(colour = "black", size = 10),
          axis.title.y = element_text(color = "black", size = 12),
          axis.line = element_line(linewidth = 1),
          legend.position = "none") +
    scale_fill_manual(values = cell_color) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA)
  
  if (!is.null(my_comparisons)) {
    p_vln <- p_vln + stat_compare_means(method = "t.test", hide.ns = FALSE,
                                        comparisons = my_comparisons,
                                        label = "p.signif", bracket.size = 0.8, size = 6)
  }
  
  print(p_vln)
  dev.off()
  cat(paste0("Violin plot saved to: ", vln_pdf, "\n\n"))
  
  return(seu_object)
}




RunSinglePathwayUCell <- function(
    seu_object,
    pathway_genes,
    pathway_name = "Pathway",
    cell_color,
    group_col = "Celltype",
    cell_levels = NULL,
    my_comparisons = NULL
) {
  
  cat(paste0(">>> Processing pathway: ", pathway_name, "\n"))
  
  # 1. Build the UCell input format
  pathway_list <- list(pathway_genes)
  names(pathway_list)[1] <- pathway_name
  
  # 2. Compute UCell scores
  seu_object <- AddModuleScore_UCell(
    seu_object,
    features = pathway_list,
    name = "_score"
  )
  
  score_name <- paste0(pathway_name, "_score")
  
  # 3. UMAP visualization
  umap_pdf <- paste0(pathway_name, "_UMAP.pdf")
  pdf(umap_pdf, width = 7, height = 5.5)
  print(FeaturePlot(
    seu_object,
    features = score_name,
    order = TRUE,
    cols = viridis::plasma(256)
  ))
  dev.off()
  cat("UMAP plot saved:", umap_pdf, "\n")
  
  # 4. Cell type distribution
  cat("Current cell distribution:\n")
  print(table(seu_object[[group_col]]))
  
  # 5. Prepare violin data
  df <- FetchData(seu_object, vars = c(group_col, score_name))
  colnames(df) <- c("Celltype", "score")
  
  if (!is.null(cell_levels)) {
    df$Celltype <- factor(df$Celltype, levels = cell_levels)
    df <- df[!is.na(df$Celltype), ]
  } else {
    df$Celltype <- as.factor(df$Celltype)
  }
  
  # 6. Violin plot
  vln_pdf <- paste0(pathway_name, "_violin.pdf")
  pdf(vln_pdf, width = 5, height = 5)
  
  p <- ggplot(df, aes(x = Celltype, y = score, fill = Celltype)) +
    geom_violin(color = "black", linewidth = 1) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = cell_color) +
    theme_classic() +
    theme(
      text = element_text(size = 10, colour = "black"),
      plot.title = element_text(hjust = 0.5, size = 15),
      axis.text.x = element_text(colour = "black", size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 10),
      axis.title.y = element_text(size = 12),
      axis.line = element_line(linewidth = 1),
      legend.position = "none"
    ) +
    labs(title = pathway_name, y = "UCell score", x = "")
  
  if (!is.null(my_comparisons)) {
    p <- p + ggpubr::stat_compare_means(
      method = "t.test",
      comparisons = my_comparisons,
      label = "p.signif",
      bracket.size = 0.8,
      size = 6
    )
  }
  
  print(p)
  dev.off()
  cat("Violin plot saved:", vln_pdf, "\n\n")
  
  return(seu_object)
}







