# =============================================================================
# 09 - UCell pathway scoring
# -----------------------------------------------------------------------------
# Description : Score a REACTOME gene signature per cell with UCell and visualize the score across groups.
# Inputs      : A Seurat object and a GMT gene set file.
# Outputs     : Per-cell UCell scores and group comparison plots.
# Depends on  : ./R/ucell_pathway_scoring.R
# Note        : All paths are relative to the repository root. Open the .Rproj
#               file (or set the working directory to the repository root)
#               before running this script.
# =============================================================================

#########Ucell########


#
source('./R/ucell_pathway_scoring.R')
library(UCell)
library(stringr)
library(viridis)
library(clusterProfiler)
library(msigdbr)
#




#####################

senescence_gs <- read.gmt("REACTOME_CELLULAR_SENESCENCE.v2026.1.Hs.gmt")
senescence_genes <- senescence_gs$gene

cell_color <- c(
  "A" = "#F3AE63",
  "M" = "#73558B",
  "Y"   = "#93CBAE"
)
# specify the comparison groups
my_comps <- list(
  c("A", "M"),
  c("A", "Y"),
  c('M','Y')
)

# run
af.ST <- RunSinglePathwayUCell(
  seu_object     = afd,
  pathway_genes  = senescence_genes,
  pathway_name   = "REACTOME_CELLULAR_SENESCENCE",
  cell_color     = cell_color,
  group_col      = "Type",
  cell_levels    = c("A", "M", "Y"),
  my_comparisons = my_comps
)


table(af.ST$dist)
