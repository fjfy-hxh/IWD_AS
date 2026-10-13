
plot_umap_with_unchull <- function(
    seurat_obj,
    reduction = "umap",
    cluster_col = "Type",
    colors = NULL,
    point_size = 0.01,
    unchull_alpha = 0.00,
    unchull_size = 0.25,
    unchull_delta = 0.3,
    axis_trunc = 2,
    aspect_ratio = 1,
    show_legend = TRUE,
    arrow_length = 0.2
) {
  if (!requireNamespace("ggunchull", quietly = TRUE)) {
    stop("Please install the ggunchull package first: devtools::install_github('const-ae/ggunchull')")
  }
  if (!requireNamespace("ggh4x", quietly = TRUE)) {
    stop("Please install the ggh4x package first: install.packages('ggh4x')")
  }
  
  # Check the input
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop("Specified reduction not found in the Seurat object: ", reduction)
  }
  if (!cluster_col %in% names(seurat_obj@meta.data)) {
    stop("Specified cluster column not found in metadata: ", cluster_col)
  }
  axis <- ggh4x::guide_axis_truncated(
    trunc_lower = unit(0, "npc"),
    trunc_upper = unit(3, "cm")
  )
  library(ggunchull)
  library(ggh4x)

  plotData <- as.data.frame(seurat_obj[[reduction]]@cell.embeddings)
  Cluster <- seurat_obj@meta.data[[cluster_col]]

  if (is.null(colors)) {
    n_clusters <- length(unique(Cluster))
    colors <- scales::hue_pal()(n_clusters)
    warning("No color vector provided, using the default ggplot2 palette")
  }

  p <- ggplot(plotData, aes(x = .data[[paste0(reduction, "_1")]], 
                            y = .data[[paste0(reduction, "_2")]], 
                            fill = Cluster, 
                            color = Cluster)) +
    stat_unchull(
      alpha = unchull_alpha, 
      size = unchull_size, 
      lty = 0, 
      delta = unchull_delta
    ) +
    geom_point(
      size = point_size, 
      show.legend = show_legend
    ) +
    theme(
      aspect.ratio = aspect_ratio,
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(
        arrow = arrow(type = "closed", length = unit(arrow_length, "inches"))
      ),
      axis.title = element_text(hjust = 0.05, face = "italic")
    ) +
    guides(
      color = guide_legend(override.aes = list(alpha = 1)),
      fill = guide_legend(override.aes = list(alpha = 1))
    ) +
    scale_x_continuous(
      breaks = NULL,
      guide = guide_axis_truncated(
        trunc_lower = unit(0, "npc"), 
        trunc_upper = unit(axis_trunc, "cm")
      )
    ) +
    scale_y_continuous(
      breaks = NULL,
      guide = guide_axis_truncated(
        trunc_lower = unit(0, "npc"), 
        trunc_upper = unit(axis_trunc, "cm")
      )
    )

  if (!is.null(colors)) {
    p <- p + 
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors)
  }

  p <- p + 
    xlab(paste0(toupper(reduction), "_1")) +
    ylab(paste0(toupper(reduction), "_2"))
  
  return(p)
}