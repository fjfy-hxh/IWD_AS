# =============================================================================
# 05 - Spatial transcriptomics main analysis (Banksy + RCTD)
# -----------------------------------------------------------------------------
# Description : Read and normalize spatial data, attach Banksy niche labels, run RCTD deconvolution to transfer cell-type annotations, and build PCA / UMAP embeddings.
# Inputs      : A spatial Seurat object (Xenium / VisiumHD) and a single-cell reference for RCTD.
# Outputs     : Annotated spatial Seurat objects (.qs), cell-type CSV tables.
# Depends on  : ./R/spatial_io_functions.R (the ST_Rprocessing.R referenced in the script; adjust the source() call accordingly)
# Note        : All paths are relative to the repository root. Open the .Rproj
#               file (or set the working directory to the repository root)
#               before running this script.
# =============================================================================

df.domain<-read.csv( "stroma_banksy_niche.csv")

af.ST$domain <- df.domain$domain[match(colnames(af.ST),
                                       df.domain$cell_id)]


library(Seurat)
library(argparse)
library(dplyr)
library(ggplot2)
library(Banksy)
library(SeuratWrappers)
library(harmony)
library(qs)
library(reticulate)
library(hdf5r)
library(sceasy)
library(BiocParallel)
library(dplyr)
library(spacexr)
library(qs)
source("ST_Rprocessing.R")


af.stroma<-qread("xenium_stroma.qs")




af.ST$cytoplasm_area <- af.ST$cell_area - af.ST$nucleus_area
af.ST$nc_ratio <- af.ST$nucleus_area / af.ST$cytoplasm_area
af.ST$nc_ratio[is.infinite(af.ST$nc_ratio)] <- NA

# remove cells with 0 counts
af.ST <- subset(af.ST, subset = nCount_Spatial > 20)
af.ST
af.ST <- NormalizeData(af.ST, normalization.method = "LogNormalize", scale.factor = 10000)
af.ST <- FindVariableFeatures(af.ST, selection.method = "vst", nfeatures = 2000)
scale.genes<-VariableFeatures(af.ST)
af.ST <-ScaleData(af.ST,features = scale.genes)
af.ST <- Seurat::RunPCA(af.ST, features = VariableFeatures(object = af.ST))

afd=FindNeighbors(afd, dims = 1:20, reduction = "pca")
afd <- RunUMAP(afd, reduction = "pca", dims = 1:20)


#############################################################
#############################################################
###########RCTD####################################################
#############################################################
table(af.ST$group)
# Extract the counts matrix (must be a sparse or regular matrix)
af
table(af$Celltype)
Idents(af)<-af$celltype
af_RCTD<- af
#
cell_counts <- table(af_RCTD$Celltype)
table(af_RCTD$Celltype)
#celltypes_to_remove <- names(cell_counts[cell_counts < 50])
#af_RCTD <- subset(af_RCTD, subset = Celltype != 'NK cell')

counts <- GetAssayData(af_RCTD, assay = "RNA", slot = "counts")

# 
cell_types <- Idents(af_RCTD)  # 

#
reference_RCTD <- Reference(counts, cell_types)


#af.ST_RCTD<-subset(af.ST_both, Cell_Anno %in% c("FBs"))
af.ST_RCTD_1<-subset(af.ST, Type %in% c("A"))
af.ST_RCTD_2<-subset(af.ST, Type %in% c("M",'Y'))


#############
af.ST_RCTD_1 <- RunRCTD_Annotations(
  seurat_obj = af.ST_RCTD_1,
  assay_name = "Spatial",
  data_type = "pcf",#visiumHD #pcf
  reference_RCTD = reference_RCTD,
  UMI_min = 20,  # 
  max_cores = 25
)

af.ST_RCTD_1$Cell_Anno[is.na(af.ST_RCTD_1$Cell_Anno)] <- "low_Qs"


############
af.ST_RCTD_1 <- subset(af.ST_RCTD_1, subset = Cell_Anno != "low_Qs")

af.ST_RCTD_2 <- RunRCTD_Annotations(
  seurat_obj = af.ST_RCTD_2,
  assay_name = "Spatial",
  data_type = "pcf",#visiumHD #pcf
  reference_RCTD = reference_RCTD,
  UMI_min = 20,  # 
  max_cores = 25
)

af.ST_RCTD_2$Cell_Anno[is.na(af.ST_RCTD_2$Cell_Anno)] <- "low_Qs"


###########
af.ST_RCTD_2 <- subset(af.ST_RCTD_2, subset = Cell_Anno != "low_Qs")


af.ST<-merge(af.ST_RCTD_1,af.ST_RCTD_2)

qsave(af.ST, file = "af_ckd_spatial.qs")
#########
af.ST_both<-sce_subset
barcodes <- colnames(af.ST)##
cell_anno_values <- af.ST$Cell_Anno

#
output_df <- data.frame(
  cell_id = barcodes,
  group = cell_anno_values,
  row.names = NULL # 
)
head(output_df)

write.csv(output_df, "af_ST_celltype.csv", row.names = FALSE)



############
af_list <- SplitObject(af.ST, split.by = "group")


#
processed_af_list <- list()

# 
for (i in 1:length(af_list)) {
  subset_name <- names(af_list)[i]
  cat("Processing subset: ", subset_name, "\n")
  
  # 
  seurat_subset <- af_list[[i]]
  
  # 
  seurat_annotated <- RunRCTD_Annotations(
    seurat_obj = seurat_subset,
    assay_name = "RNA",
    reference_RCTD = reference_RCTD,
    UMI_min = 50,
    max_cores = 50
  )
  
  # 
  processed_af_list[[subset_name]] <- seurat_annotated
}
# 
af_combined <- processed_af_list[[1]]
# 
for (j in 2:length(processed_af_list)) {
  af_combined <- merge(af_combined, processed_af_list[[j]])
}


qsave(af_combined, file = "slide1_RCTD.qs")



