# =============================================================================
# 08 - UMAP aesthetics
# -----------------------------------------------------------------------------
# Description : Draw outline (unchull) UMAPs and cloud-style UMAPs with customized legends and axes.
# Inputs      : A Seurat object with a UMAP embedding.
# Outputs     : UMAP PDF / PNG figures.
# Depends on  : ./R/umap_with_unchull.R
# Note        : All paths are relative to the repository root. Open the .Rproj
#               file (or set the working directory to the repository root)
#               before running this script.
# =============================================================================

#http://www.360doc.com/content/24/0430/20/82995563_1122020662.shtml


#############################
source("./R/umap_with_unchull.R")


###small arrow styled plot####
#plotData$UMAP_1
p <- plot_umap_with_unchull(
  seurat_obj = af,
  cluster_col = "Celltype",
  colors = colors3,
  point_size = 0.03#0.2-0.3
)

print(p)

pdf(file = " UMAP.pdf",width =6,height = 4)
p
dev.off()

###cloud-like UMAP
p1 <- plotUMAP(
  seurat_obj = af.epi,
  reduction = "umap",
  cluster_col = "Type",
  color_palette = colors3,
  axis_line_size = 0.3,      # 
  legend_point_size = 3,     #
  legend_text_size = 5,      # 
  legend_ncol = 1 ,          # 
  legend_alpha = 0.7 ,
  point_size = 0.2
)
p1
ggsave(filename ="Episub UMAP .png",height = 3,width = 4,plot = p1,dpi = 300,bg="white")





