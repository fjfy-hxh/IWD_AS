add.flag <- function(pheatmap, kept.labels, repel.degree = 1) {
  
  library(grid) # Ensure the grid package is loaded
  heatmap <- pheatmap$gtable
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # Keep only labels in kept.labels; set the rest to empty
  new.label$label <- ifelse(new.label$label %in% kept.labels, new.label$label, "")
  
  # Compute evenly distributed new y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # Use a native grid function to safely convert complex units to plain numeric values (npc)
    strip.npc.safe <- function(dd) {
      as.numeric(grid::convertY(dd, "npc", valueOnly = TRUE))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc.safe(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc.safe(d[d.select][i]))
    
    # Recompute the repelled coordinates
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  
  new.y.positions <- repelled.y(new.label$y, d.select = new.label$label != "")
  
  # Draw the flag segments
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # Shift the positions of the retained labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # Add the flag segments to the gtable
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4)
  
  # Replace the original labels
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # Plot
  grid.newpage()
  grid.draw(heatmap)
  
  invisible(heatmap)
}






###################################
#' Plot text-box panels for multi-module enrichment analysis
# #########################
plot_module_text_boxes <- function(gene_list, 
                                   species = "mouse", 
                                   rowcolor = NULL, 
                                   total_terms = 20,
                                   wrap_width = 38,
                                   text_size = 3.8) {
  
  library(clusterProfiler)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  
  # 1. Configure the species database
  if (species == "mouse") {
    library(org.Mm.eg.db)
    org_db <- "org.Mm.eg.db"
  } else if (species == "human") {
    library(org.Hs.eg.db)
    org_db <- "org.Hs.eg.db"
  } else {
    stop("Species must be 'human' or 'mouse'")
  }
  
  # Unify column names
  if("cluster" %in% colnames(gene_list)) gene_list <- rename(gene_list, Module = cluster)
  
  # 2. Loop over GO enrichment analysis
  message("Running Enrichment Analysis...")
  Module_GO <- data.frame()
  unique_modules <- unique(gene_list$Module)
  
  for (i in unique_modules) {
    curr_genes <- gene_list %>% filter(Module == i) %>% pull(gene)
    
    # ID conversion
    df_id <- bitr(curr_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)
    
    if(nrow(df_id) > 0) {
      go <- enrichGO(gene          = unique(df_id$ENTREZID),
                     OrgDb         = org_db,
                     keyType       = 'ENTREZID',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05, # 👈 If errors persist, try relaxing this to 0.1
                     qvalueCutoff  = 0.05, # 👈 Or relax this to 0.2
                     readable      = TRUE)
      
      # Ensure the enrichment result exists and is not empty
      if (!is.null(go) && nrow(go@result) > 0) {
        # Manually filter significant results
        res <- go@result %>% filter(p.adjust < 0.05)
        if (nrow(res) > 0) {
          res <- res %>% mutate(cluster = as.character(i)) # Convert to character to avoid errors
          Module_GO <- rbind(Module_GO, res)
        }
      }
    }
  }
  
  # ✨ Key fix: check whether any enrichment result actually exists
  if (nrow(Module_GO) == 0) {
    stop("🚨 Analysis aborted: no significantly enriched GO terms found!\nPossible reasons: too few genes, no coherent function, or thresholds that are too strict. Try relaxing pvalueCutoff in enrichGO.")
  }
  
  # 3. Allocate term quotas by gene proportion
  message("Allocating Term quotas...")
  module_counts <- as.data.frame(table(gene_list$Module))
  colnames(module_counts) <- c("cluster", "gene_count")
  module_counts$cluster <- as.character(module_counts$cluster) # Force to character for safe matching
  
  module_counts <- module_counts %>%
    mutate(ratio = gene_count / sum(gene_count),
           term_limit = pmax(1, round(ratio * total_terms)))
  
  final_selected_terms <- data.frame()
  for (i in unique(Module_GO$cluster)) {
    limit <- module_counts$term_limit[module_counts$cluster == i]
    temp_df <- Module_GO %>%
      filter(cluster == i) %>%
      arrange(qvalue) %>%
      head(limit)
    final_selected_terms <- rbind(final_selected_terms, temp_df)
  }
  
  # 4. Prepare plotting data
  if (is.null(rowcolor)) {
    rowcolor <- hue_pal()(length(unique_modules))
  }
  
  cluster_df <- final_selected_terms %>%
    mutate(Description = str_wrap(Description, width = wrap_width)) %>% 
    group_by(cluster) %>%
    summarise(
      combined_text = paste(Description, collapse = "\n"),
      n_lines = sum(str_count(Description, "\n") + 1)
    ) %>%
    ungroup() %>%
    arrange(cluster)
  
  # Coordinate allocation logic
  cluster_df <- cluster_df %>%
    mutate(
      spacing = (n_lines * 2.5) + 15, 
      y_offset = lag(cumsum(spacing), default = 0),
      y_pos = -y_offset 
    )
  
  cluster_df$border_color <- rowcolor[1:nrow(cluster_df)]
  cluster_df$fill_color <- alpha(cluster_df$border_color, 0.05)
  
  # 5. Plot
  message("Generating Plot...")
  p <- ggplot(cluster_df, aes(x = 0, y = y_pos)) +
    geom_label(
      aes(label = combined_text, color = border_color, fill = fill_color),
      family = "sans", # Use a common sans-serif font for better system compatibility
      size = text_size,
      hjust = 0, vjust = 1,
      lineheight = 0.9,
      label.padding = unit(0.8, "lines"), 
      label.size = 1.2,
      label.r = unit(0.2, "lines")
    ) +
    scale_color_identity() +
    scale_fill_identity() +
    xlim(-0.05, 1) +
    expand_limits(y = min(cluster_df$y_pos) - max(cluster_df$n_lines)*3 - 20) +
    theme_void() +
    theme(plot.margin = margin(10, 10, 10, 10))
  
  return(p)
}

######################################################
######################################################
##Compute matrix
######################################################
GetDistanceSmoothMatrix <- function(seu_obj, distance_vals, n_bins = 100, 
                                    smooth_df = 3, logfc_thresh = 0.25, 
                                    min_pct = 0.25, assay = NULL) {
  
  if (is.null(assay)) {
    assay <- DefaultAssay(seu_obj)
    message(sprintf("No assay specified, using the current default assay: '%s'", assay))
  }
  
  # 1. Sort and reorder the Seurat object by distance
  sorted_indices <- order(distance_vals)
  seu_sorted <- seu_obj[, sorted_indices]
  distances_sorted <- distance_vals[sorted_indices]
  
  # 2. Split into n_bins equal-width segments
  distance_range <- range(distances_sorted, na.rm = TRUE)
  breaks <- seq(distance_range[1], distance_range[2], length.out = n_bins + 1)
  
  # ✨ Key fix: change "Bin_" to "Bin-" to avoid match failures caused by Seurat forcing renames
  group_id_bins <- cut(distances_sorted, 
                       breaks = breaks, 
                       labels = paste0("Bin-", sprintf("%03d", 1:n_bins)),
                       include.lowest = TRUE)
  seu_sorted$distance_bin <- group_id_bins
  
  # ✨ Matching fix: also use a hyphen for the binary group names
  median_distance <- median(distances_sorted, na.rm = TRUE)
  group_id_binary <- ifelse(distances_sorted <= median_distance, "Low-Distance", "High-Distance")
  seu_sorted$distance_binary <- group_id_binary
  
  # 3. Compute differentially expressed genes
  message("Finding differentially expressed genes between Low-Distance and High-Distance groups...")
  Idents(seu_sorted) <- "distance_binary"
  markers <- FindMarkers(seu_sorted,
                         ident.1 = "Low-Distance",
                         ident.2 = "High-Distance",
                         min.pct = min_pct,
                         logfc.threshold = logfc_thresh,
                         assay = assay)
  
  valid_genes <- rownames(markers)
  if (length(valid_genes) == 0) stop("No significant differential genes found! Try lowering logfc_thresh.")
  
  # 4. Compute mean expression per bin
  message("Computing mean expression per bin...")
  Idents(seu_sorted) <- "distance_bin"
  avg_expr_list <- AverageExpression(seu_sorted, features = valid_genes, 
                                     assays = assay, slot = "data")
  avg_expr <- avg_expr_list[[assay]]
  
  # Extract the Bin columns that exist
  valid_bins <- levels(group_id_bins)[levels(group_id_bins) %in% colnames(avg_expr)]
  
  # ✨ Add a safeguard: if matching fails, throw a plain, candid error immediately
  if (length(valid_bins) == 0) {
    stop("Fatal error: cannot match bin groups in the matrix; check the grouping naming behaviour of your Seurat version.")
  }
  
  plot_matrix <- as.matrix(avg_expr[, valid_bins])
  
  # Remove zero-variance genes (genes with no variation across all Bins, another cause of smooth.spline errors)
  gene_vars <- apply(plot_matrix, 1, var)
  plot_matrix <- plot_matrix[gene_vars > 0, , drop = FALSE]
  
  if (nrow(plot_matrix) < 2) {
    stop("Too few genes left after removing zero-variance genes; cannot plot.")
  }
  
  # 5. Data smoothing and Z-score
  message("Smoothing curves and applying Z-score normalization...")
  plot_matrix <- t(apply(plot_matrix, 1, function(x) {
    smooth.spline(x, df = smooth_df)$y
  }))
  plot_matrix <- t(apply(plot_matrix, 1, function(x) {
    (x - mean(x)) / sd(x)
  }))
  
  plot_matrix[is.na(plot_matrix)] <- 0
  message(sprintf("✅ Matrix processing done! Output dimensions: %d genes x %d bins", nrow(plot_matrix), ncol(plot_matrix)))
  
  return(plot_matrix)
}


#########################################################
#' Plot the heatmap of spatial-distance fitted trajectories
#'
#' @param plot_matrix Expression matrix after smoothing and normalization
#' @param k_clusters Number of gene clusters to split into (default 4)
#' @param show_rownames Whether to show row names (gene names); set to FALSE when there are many genes
#' @param title Heatmap title
#' @param cluster_colors Color vector for the cluster sidebar
#' @param heatmap_colors Colors for the heatmap body
#' @param order_colors Color gradient for the top distance bar
#' @return A list containing the heatmap object (heatmap) and gene cluster information (gene_clusters)
PlotDistanceHeatmap <- function(plot_matrix, 
                                k_clusters = 4,
                                show_rownames = FALSE,
                                title = "Distance (near to far)",
                                cluster_colors = c('#F3AE63','#73558B','#93CBAE','#E2BECB', '#82B0D2', '#FFBE7A'),
                                heatmap_colors = colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                                order_colors = colorRampPalette(c('#313536', '#57147D','#FB8760','#FCE9AC'))(100)) {
  
  # 1. Define the ordering callback function (Callback)
  callback <- function(hc, mat){
    sv <- svd(t(mat))$v[,1]
    dend <- reorder(as.dendrogram(hc), wts = sv)
    as.hclust(dend)
  }
  
  # 2. Generate a pre-clustering heatmap (silent, not drawn) only to obtain the clustering tree (tree_row)
  message("Performing hierarchical clustering of genes...")
  p_temp <- pheatmap::pheatmap(plot_matrix, 
                               silent = TRUE,           # Key: run silently, do not draw this layer
                               cluster_cols = FALSE, 
                               cluster_rows = TRUE, 
                               clustering_method = "ward.D2",
                               cutree_rows = k_clusters,
                               clustering_callback = callback)
  
  # 3. Extract row annotation (Cluster)
  annotation_row <- data.frame(Cluster = factor(cutree(p_temp$tree_row, k_clusters)))
  rownames(annotation_row) <- rownames(plot_matrix)
  
  # Prepare row annotation colors
  if (k_clusters > length(cluster_colors)) {
    warning("Fewer cluster_colors than k_clusters provided; colors will be extended automatically.")
    cluster_colors <- colorRampPalette(cluster_colors)(k_clusters)
  }
  rowcolor <- cluster_colors[1:k_clusters]
  names(rowcolor) <- levels(annotation_row$Cluster)
  
  # 4. Extract column annotation (Order - distance from near to far)
  if (is.null(colnames(plot_matrix))) {
    colnames(plot_matrix) <- paste0("C", 1:ncol(plot_matrix))
  }
  col_ann <- data.frame(Order = 1:ncol(plot_matrix))
  rownames(col_ann) <- colnames(plot_matrix)
  
  # 5. Combine the color dictionary
  ann_colors <- list(
    Cluster = rowcolor,
    Order = order_colors
  )
  
  # 6. Draw the final heatmap
  message("Drawing the final heatmap...")
  p_final <- pheatmap::pheatmap(plot_matrix, 
                                cluster_cols = FALSE, 
                                cluster_rows = TRUE, 
                                show_rownames = show_rownames, 
                                show_colnames = FALSE, 
                                clustering_method = "ward.D2",
                                cutree_rows = k_clusters,
                                border_color = NA,
                                fontsize_row = 8,
                                fontsize_main = 22, 
                                color = heatmap_colors,
                                annotation_row = annotation_row,
                                annotation_col = col_ann,
                                annotation_colors = ann_colors,
                                annotation_names_row = FALSE,
                                annotation_names_col = FALSE,
                                clustering_callback = callback,
                                main = title)
  
  # 7. Return the heatmap object and gene cluster data
  return(list(
    heatmap = p_final,
    gene_clusters = annotation_row
  ))
}
