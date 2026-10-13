# =============================================================================
# 06 - Spatial distance-composition analysis
# -----------------------------------------------------------------------------
# Description : Merge distance metadata from CSV files, bin cells by distance percentile and plot domain composition per distance bin and cell composition along a distance axis.
# Inputs      : A spatial Seurat object and pre-computed distance CSV files.
# Outputs     : Stacked domain percentage plots, distance composition curves.
# Depends on  : ./R/ST_distance_source.R
# Note        : All paths are relative to the repository root. Open the .Rproj
#               file (or set the working directory to the repository root)
#               before running this script.
# =============================================================================


source("./R/ST_distance_source.R")
library(ggplot2)
library(dplyr)
library(viridis)

#########################################################################
##
#####################################################
##
##################################################
af.ST_Fib1 <- add_metadata_from_csv(
  csv_file = "./neighborhood_number.csv",  # 
  seurat_obj = af.ST_Fib1,         # 
  metadata_cols = c("Myeloid_cells" )  # 
)
af.ST_Fib1$Celltype_Epi_Fib
seu_radiusrange<- subset(af.plaque_radiusrange, Type %in% c('CWD'))
seu_radiusrange<- subset(seu_radiusrange, subset = domain != "Peri vascular")

seu_radiusrange <- add_distance_percentile(seu_radiusrange,
                                           dist_col = "min_distance", 
                                           group_col = "group", 
                                           probs = seq(0, 1, 
                                                       0.2),#
                                           new_col_name = "distance_percentile")


###
custom_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")

p2 <- plot_domain_percentage(
  seurat_obj = seu_radiusrange,
  radius_range_col = "distance_percentile",##the x axis must be a categorical variable
  domain_col = "domain",##y axis proportion
  x_labels = c("0-20%", "20-40%", "40-60%", "60-80%",'80-100%'),#adjust the x-axis categories
  color_palette = custom_colors,
  xlab = "Radius Range",
  ylab = "Percentage (%)",
  title = " ",
  legend_title = "Domain Types"
)
print(p2)
###
#####################################################
##

###################################
#
#######################################
af.ST_vertical_distance <- add_metadata_from_csv(
  csv_file = "./EPI_vertical_distance.csv",  # 
  seurat_obj = af.ST,         # 
  metadata_cols = c("vertical_distance",'dist')  # 
)



# 
p_prop <- plot_cell_composition_by_distance(
  seurat_obj   = af.niche6,
  distance_col = "vertical_distance",
  cell_col     = "Type",
  bin_num      = 60,               
  plot_type    = "proportion",
  colors=colors_Type
)
print(p_prop)
table(af.ST$Cell_Anno)

# 

af.niche6<- subset(af.ST, Cell_Anno %in% c('Plasma cells'))
colors_Type
# 
my_colors <- c("A" = "#F3AE63", "M" = "#73558B","Y" = "#93CBAE")

# 

plot_fast_curve(af.ST, 
                x_var = "vertical_distance",
                y_var = "factor5", 
                group_var = "Type",#use_relative = T,
                colors = my_colors)

