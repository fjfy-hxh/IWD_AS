# =============================================================================
# 07 - Spatial visualization and data export
# -----------------------------------------------------------------------------
# Description : Export cell-type tables, plot gene expression across samples, spatial DimPlots and density plots.
# Inputs      : A spatial Seurat object.
# Outputs     : Cell-type CSV, per-sample expression plots, spatial PNG/PDF figures.
# Depends on  : ./R/STdistance.R (external helper library, not bundled - see README)
# Note        : All paths are relative to the repository root. Open the .Rproj
#               file (or set the working directory to the repository root)
#               before running this script.
# =============================================================================


af.output<-  subset(af.stroma, CN == "CN_5")

# 
export_data <- data.frame(
  cell_id = colnames(af.ST),    #
  group = af.ST$Celltype_Epi_Fib  # 
)

# 
write.csv(export_data, 
          file = "Celltype_Fibsub.csv", # 
          row.names = FALSE)         # 

###########################################################################
source("STdistance.R")


#####dimplot plotting##########
table(af.ST$group)


af.niche <- subset(af.ST, niche %in% c('niche_5','niche_3'))
af.ST$vertical_distance
# 2. add centroid coordinates into a reduction layer
ks_plot_gene_across_samples(
  seurat_obj = af.ST,
  feature    = "vertical_distance",
  group_col = "group",
  group_ids  = c("A6"),
  ncol       = 1,   # lay them out in a single row
  pt_size    = 0.1,
  colors=c('#313536', '#61BECD','#89BF69','#E1D432','#C3400C','#821407')

)

####################
 ggsave(   filename = "vertical distance .png",   width = 7, height   = 5, units    = "in", dpi  = 300 )
 

p<-DimPlot(af.treat_filter, reduction = "spatial", label = F, label.size = 2.5,
        group.by = "Signal_Status",
        pt.size = 0.01,cols= colors3)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
p <- p + coord_fixed(ratio = 1)
p
# save as PNG (handy for sharing / preview)
ggsave(
  filename = "A6 celltype.png",
  plot     = p,
  width    = 7,      # 
  height   = 7,      # 
  units    = "in",
  dpi      = 300
)

plot_density(af.plot, features = c("FOXE1"),reduction = 'spatial',slot='data',
             pal = 'magma',method='ks',raster = T, size = 0.8) &
  theme_blank()&#
  theme(legend.frame = element_rect(colour = "black"),#
        legend.ticks = element_line(colour = "black", linewidth  = 0),#
        legend.key.width = unit(0.3, "cm"),#
        legend.key.height = unit(0.8, "cm"),#
        legend.title = element_text(color = 'black', face = "bold", size=8))#
