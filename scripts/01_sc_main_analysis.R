# =============================================================================
# 01 - Single-cell main analysis pipeline
# -----------------------------------------------------------------------------
# Description : Quality control and filtering, SCTransform normalization, PCA and Harmony integration, graph-based clustering, marker gene detection, manual cell-type annotation, cell proportion and metabolic pathway analysis.
# Inputs      : A multi-sample Seurat object for single-cell RNA-seq.
# Outputs     : QC violin / correlation plots, elbow plot, cluster and annotated UMAPs, marker gene CSV, dot plots, cell proportion and metabolic pathway plots.
# Depends on  : ./R/ks_scAverExp.R, ./R/ks_scDotplot.R, ./R/Singlecellratio_plotstat.R, ./R/ks_cellration_plot.R, ./R/Mouse_sc_metabolism.R
# Note        : All paths are relative to the repository root. Open the .Rproj
#               file (or set the working directory to the repository root)
#               before running this script.
# =============================================================================

rm(list = ls())
options(stringsAsFactors = F)
gc()

af
af[["RNA"]] <- as(af[["RNA"]], "Assay")#

af[["RNA"]] <- as(object = af[["RNA"]], Class = "Assay5")#


colors_type<-c("royalblue",  "indianred")

nmf_color<-c('#9C0642','#5FB6A9','#614DA1','#CD3B4E','#EA6547',
             '#F79C5A','#FACF7E','#FCF2A8','#F1FBAE','#D0ED9D','#97D5A5'
             )

colors_all<-c('#B7B663','#9d9ac6',"#BBE379",'#960937',"#0F7240","#edc951")

colors_domain<-c('#00f3d8','#ffabc3','#8c86ff','#bf9800')
colors_niche<-c('#87C164',"#ffabc3",'#ea7400','#b4a100','#43bce7',"#00bea8",
                '#cdb6ff','#6c94ff',"#ff6464",'#f237ff' )

colors_Celltype<-c('#d51f26','#f3758a','#ff87bbff',"#fbaeb5",#EC
               '#0F7240', #T
              '#B7B663',#B
              '#d5a339',"#f8ff00",'#fccc3a',#Mac
              "#9d9ac6", #SMC
              "#bbe379",  "#00a7fa",'#b9e6fb','#354c85',
              "#b1ddcc", "#C9E1AF", "#31aaaa"#FMC
           )
colors_FMC<-c("#bbe379",  "#00a7fa",'#b9e6fb','#354c85',
              "#b1ddcc", "#C9E1AF", "#31aaaa"#FMC
)
colors_meta<-c( "#00a7fa",'#354c85','#31aaaa',"#9d9ac6"#FMC
)

colors_carotid<-c("#0F7240","#BBE379","#edc951",'#960937',"#79B99D","#95A6DA"  )


library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(patchwork)
library(ggplot2)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)
library(ggpubr)
library(Biobase)
library(Matrix)
library(RcppML) 
library(magick)
library(GeneNMF)#R>=4.3
library(scRNAtoolVis)
library(Nebulosa)
library(reshape2)
library(ggnetwork)
library(ggimage)
library(mascarade)
library(UCell)
library(ggrastr)
library(qs)

af.bei<-af
# 
pdf(file = "01.vlnplot.pdf",width = 8,height = 5)
VlnPlot(af, features = c("nFeature_RNA", "nCount_RNA"), ncol = 4)
dev.off()

# 
plot1 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "percent.mt")+ RotatedAxis()
plot2 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "percent.rb")+ RotatedAxis()
plot3 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ RotatedAxis()
#
pdf(file = "01.corqc.pdf",width =12,height = 5)
plot1+plot2+plot3+plot_layout(ncol = 3)      #plot_layout, patchwork
dev.off()

plot3
# 
exp_m <- GetAssayData(af, layer = 'counts')
gene.all <- rownames(exp_m)
mt.gene  = grep('^MT-', gene.all, value = TRUE)#
ribosome.gene = grep('^RPL|^RPS|^MRPS|^MRPL', gene.all, value = TRUE)#
features=setdiff(gene.all,c(mt.gene,ribosome.gene))#
filtered_counts <- GetAssayData(af, layer = 'counts')[features, ]# 
af <- SetAssayData(af, layer = 'counts', new.data = filtered_counts)# 


table(af$Celltype)

options(future.globals.maxSize = 4 * 1024^3)
af <- SCTransform(af, assay = "RNA", verbose = FALSE)

af <- Seurat::RunPCA(af)

pct <- af[["pca"]]@stdev / sum(af[["pca"]]@stdev) * 100 ; cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
ElbowPlot(af)$data %>% ggplot() +
  geom_point(aes(x = dims,y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred") +
  theme_bw() + labs(title = "Elbow plot: quantitative approach")

af <- Seurat::RunPCA(af,npcs = 17)
af <- RunHarmony(af, group.by.vars = "Type")

rownames(af.ST)


table(afe$leiden)
stage_vector <- rep(NA, length(afe$leiden))

af <- FindNeighbors(af, dims = 1:17,reduction = "harmony")

af <- FindClusters(af, resolution = 0.6)

# 
head(Idents(af), 5)

# 
head(af@meta.data)
table(af@meta.data$seurat_clusters)
# 

af.markers <- FindAllMarkers(afd, only.pos = T, min.pct = 0.25, 
                             logfc.threshold = 0.1)
write.csv(af.markers,file = "FMC cluster_markers.csv")
af.markers_cluster <- af.markers[af.markers$cluster == 8, ]

myMarkers=c("COL4A1",'TIMP2','APP','CCN1','PECAM1','SYNM','SPTBN1')
pdf(file = "cluster.pdf",width =4.5,height = 4)
jjVolcano(diffData = af.markers, aesCol = c('#42AFFC','tomato'),topGeneN=2,
          #myMarkers=myMarkers,
          flip = F,legend.position = c(0.9,0.2),#polar=T,
          tile.col = c('#9d9ac6','#00a7fa','#354c85','#31aaaa'))##cluster
dev.off()

## 
af <- RunUMAP(af, dims = 1:17, reduction = "harmony",min.dist = 0.7)

# 
pdf(file = "05-cluster.UMAP.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "umap", label = F, label.size = 2.5,pt.size = 1,
       group.by = "seurat_clusters",
        cols= colors2)+#ggsci::scale_color_igv()+
  theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()

af$seurat_clusters<-Idents(af)

####################

pdf(file = "08.ann_cluster_marker.pdf",width =20,height = 7)
do_DotPlot(sample = af,features = genes,dot.scale = 12,
           legend.length = 50,
           legend.framewidth = 2, font.size =12)
dev.off()

#

af <- FindClusters(af, resolution = 0.6)
table(af@active.ident)
afe$`MP1`
ann.ids <- c("T cell",  #0Injured Proximal tubular cell  
             "FBs like",  #1Injured Proximal tubular cell
             "vSMC",  #2
             "FBs like", #3
             "Mono/Mac", #4
             "Mono/Mac",#5
             "EC",#6
             "Mono/Mac",#7
             "T cell",#8
             "Mono/Mac",#9
             "Plasma cell",#10
             "Plasma cell",#11
             "FBs like",#12
             "Plasma cell",#13
             "Mast"#14
)


####

Idents(af)

afidens=mapvalues(Idents(af), from = levels(Idents(af)), to = ann.ids)
Idents(af)=afidens

af$Celltype=Idents(af)
########
pdf(file = "anno.UMAP.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "umap", label = F, label.size = 2.5,pt.size = 0.0001,
        group.by = "Celltype",
        cols= colors2)+#ggsci::scale_color_igv()+
  theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()




####################

colors<-c("#FB8072","gray")

pdf(file = "controlneu.UMAP.pdf",width =7,height = 5.5)
DimPlot(af.neu ,reduction="umap", group.by = "Type",
        label =F)+
  scale_color_manual(values=colors)+
  NoLegend()
dev.off()

source(file = "./R/ks_scAverExp.R")
source(file = "./R/ks_scDotplot.R")

afgenes=c('Pecam1','Vwf','Acta2','Myh11','Col1a1','Dcn','Cd3d',
          'Nkg7','Cd79a','Cd68','Cd14')

afgenes=c('Tnf','Il1b','Il6','Spp1','Ly6C','Fcgr4','Itgam','S100a9','S100a8',
          'Trem2','Apoe','Fabp5','Lyve1','Folr2','Nampt')
afgenes=c('Icam1','Vcam1','Cdh5',
          'Klf2','Klf4','Nos3','Pecam1','Vwf','Acta2','Vim','Tagln',
          'Apoe')
table(af$Celltype)
featureSets <- list("Inf EC" = c('Icam1'),#, "Mzb1"
                    "EndMT EC" = c('Vcam1'),#Endothelial cell
                    'Angio EC'=c('Tmem252'),
                    'Res EC'=c("Pecam1"),
                    'T cell'=c('Cdh5'),
                    'B cell'=c('Cd3d'),
                    'Res Mac'=c('Cd79a','Cd14','CD68'),
                    'LAMs'=c('Lyve1'),
                    "Inf Mac"=c('Folr2'),#Smooth muscle cell
                    "vSMC"= c('Trem2'),#Macrophage
                    "FMC_1" = c('Apoe'),#neutrophil
                    "FMC_2" =c('Fabp5'),
                    "FMC_3" = c('Tnf'),
                    "FMC_4" = c('Il1b'),
                    "FMC_5" = c('Acta2','Myh11'),
                    "FMC_6" = c('Col1a1'),
                    "FMC_7" = c('Dcn')#nature killer
)
featureSets <- list(  'T cell'=c("CD8A",'CD3D'),
                      'FBs like'=c('COL1A1','DCN'),
                      "vSMC" = c('ACTA2','MYH11'),
                      'Mono/Mac'=c('CD68','CD14'),
                      "EC" = c('PECAM1','VWF'),#, "Mzb1"
                      'Plasma cell'=c('MZB1','XBP1'),
                      'Mast cell'=c('TPSAB1','MS4A2')
                   
)

sce<-af
sce$celltype<-sce$Celltype
Exp_scRNA <- ks_scAverExp(obj = sce,#seurat or ArchR proj
                          features = featureSets,
                          CellTypes = unique(sce$celltype))

scRNA_dotplot = ks_scDotplot(obj = sce,#seurat or ArchR proj
                             features = featureSets,
                             CellTypes = unique(sce$celltype),
                             avgPctMat= Exp_scRNA,
                             pal = colors_carotid )


x<-c( "#960937","#9d9ac6","#BBE379","#B7B663","#0F7240" , "#edc951" )
########################

# only.pos:
af.markers <- FindAllMarkers(afd, only.pos = T, min.pct = 0.25, logfc.threshold = 0.75)

# get top 10 genes
topaf.markers <- af.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)


# plot
pdf(file = "09-cell_marker.hetmap.pdf",width =7,height = 7)
AverageHeatmap(object=afd,
               markerGene = topaf.markers$gene,slot = "data",assays='RNA',
               clusterAnnoName =F,annoCol=T,myanCol= colors2[1:2],# htRange=c(-0.5, 0, 2),
               htCol =c('#435B95','white','#C3400C')
               )



afgenes=c('Spp1','Lgals3','Vcam1','Icam1','Thbs1','Col1a1','Col3a1','Fn1',
          'Postn','Mmp2','Mmp14','Mmp3',
          'Cx3cl1','Cxcl12','Mif')
afgenes=c('Acss2','Acly','Ldha','Hif1a')
afgenes=c('Pecam1','Vwf','Acta2','Myh11','Col1a1','Dcn','Cd8a',
          'Nkg7','Cd79a','Cd68','Cd14')

afgenes=c('Tnf','Il1b','Il6','Spp1','Ly6C','Fcgr4','Itgam','S100a9','S100a8',
          'Trem2','Apoe','Fabp5','Lyve1','Folr2','Nampt')
afgenes=c('Icam1','Vcam1','Cdh5',
          'Klf2','Klf4','Nos3','Pecam1','Vwf','Acta2','Vim','Tagln',
          'Apoe')

afgenes=c('MP1','MP2','MP3','MP4','MP5','MP6','MP7','MP8','MP9','MP10','MP11','MP12')
afgenes=c('Sox9','Fosl1','Creb3l1','Klf9')

afgenes=c('PECAM1','VWF','ACTA2','MYH11','COL1A1','DCN','CD8A','CD3D',
          'NKG7','CD79A','CD68','CD14','TPSAB1','MS4A2','XBP1','MZB1')
afgenes=c('FMC_2','FMC_4','FMC_7','LAMs')

Idents(afe)<-afe$oldType
Idents(afd)<-afd$Celltype
pdf(file = "spatial_marker_dotplot.pdf",width =7,height = 3)
DotPlot(af, features = afgenes)+theme_bw()+#

theme(panel.grid = element_blank(),
           
axis.text.x=element_text(angle=45,hjust = 1,vjust=1))+#
  
scale_color_gradientn(values = seq(0,1,0.2),
                      colours = c('#64A2CC',"white",'#E07B56',"#731929"))+#
 # c('#64A2CC',"white",'#E07B56',"#731929")
labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+coord_flip()
dev.off()
  ###################

################################
Idents(af)<-af$Celltype


#################
#
afgenes=c('Tnf','Il1b','Il6','Spp1','Ly6C','Fcgr4','Itgam','S100a9','S100a8',
          'Trem2','Apoe','Fabp5','Lyve1','Folr2','Nampt')
afgenes=c('Icam1','Vcam1','Cdh5',
          'Klf2','Klf4','Nos3','Pecam1','Vwf','Acta2','Vim','Tagln',
          'Apoe')
afgenes=c('Dcn','Col1a1','Acta2','Myh11','Pdgfra','Il6','Cxcl12')
afgenes=c('Acta2','Myh11','Tagln','Lum','Col1a1','Cd68')

afgenes=c('Icam1','Vcam1','Tmem252','Pecam1','Cdh5','Cd3d','Cd79a','Cd68','Cd14',
          'Lyve1','Folr2','Trem2','Apoe','Fabp5','Tnf','Il1b','Col1a1','Lum',
          'Tagln','Myh11')

#afgenes=read.table("stroke_geneFA.txt",header = F,sep = "\t")[,1]
pdf(file = "Mac_subcluster_Il1b.pdf",width =6,height = 5)
plot_density(afd, features = c("MYH11"),reduction = "umap",
             pal = 'magma', raster = T, size = 0.8) &
  theme_blank()&#
  theme(legend.frame = element_rect(colour = "black"),#
        legend.ticks = element_line(colour = "black", linewidth  = 0),#
        legend.key.width = unit(0.3, "cm"),#
        legend.key.height = unit(0.8, "cm"),#
        legend.title = element_text(color = 'black', face = "bold", size=8))#
dev.off()

table(af.ST$Cell_Anno)


##### 
source("./R/Singlecellratio_plotstat.R")
table(af.ST$Cell_Anno)

af.all<-af
af<-af.ST
af<- subset(af.ST, Cell_Anno %in% c('FMC_C1','FMC_2','FMC_C3','FMC_4','FMC_C5','FMC_C6','FMC_7'))

af$celltype<-af$Celltype
af$orig.ident<-af$Type########


Singlecellratio_plotstat(af, group_by = "Type",###
                         meta.include = c("Type","orig.ident"),#
                         color_by = 'cell.type'###
)+scale_fill_manual(values=colors_meta )#colors_prcc  LEC_color

ggsave("carotid_cell_proportion_annotated.pdf",width = 5,height = 4.5)

######Second method###############################
source(file = "./R/ks_cellration_plot.R")
library(jjAnno)
af.LEC$LEC


ks_scFraction(group = af$predicted.id,
              celltype = af$Type,pal = colors_type,#CN_color4<-c("#c79000", "#93CBAE" ,"#00b1da" ,"#cdb6ff")
              colorLab =  colors_FMC, Total = T)#colors5 group  colors_prcc
ggsave("carotid_scRNA_proportion.pdf",width = 5.5,height = 3.5)

afd<-af.plaque
########
afd$orig.ident<-afd$group###
afd$celltype<-afd$domain  ####
my_comparisons<-list(c("IWD","CWD"))
Singlecellratio_plotstat(afd, group_by = "Type",
                         custom_fill_colors = c("royalblue",  "indianred"),
                         meta.include = c("Type","orig.ident"),
                         comparisons = my_comparisons, color_by = 'Type',
                         group_by.point = "orig.ident",label.x = 1, pt.size = 3,
                         label = 'p.format', ncol =4)
ggsave("single niche cellratio_diff.pdf",width = 12,height =3)

# Mouse single-cell metabolic pathways
library(stringr)
library(Seurat)
library(rsvd)
library(qs)
library(scMetabolism)
library(BiocParallel)
register(MulticoreParam(workers = 4, progressbar = TRUE))
library(Seurat)
library(nichenetr)
library(dplyr)
source('./R/Mouse_sc_metabolism.R')
# Requires V4 data
res <- Mouse.sc.metabolism(af, metabolism.type = 'KEGG')
#########
#
################################
#########################
afd

af[["RNA"]] <- as(af[["RNA"]], "Assay")

sub_scRNA <- af

res <-sc.metabolism.Seurat(obj = af,
                           method = "AUCell", # VISION, AUCell, ssgsea and gsva
                           imputation =F, ncores = 12, 
                           metabolism.type = "KEGG") # 
countexp.Seurat <- res

#
score <- countexp.Seurat@assays$METABOLISM$score
score[1:4,1:4]
#
score_change <- score %>% 
  select_all(~str_replace_all(., "\\.", "-"))  #
#
identical(colnames(score_change) , rownames(countexp.Seurat@meta.data))
#[1] TRUE
countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data,t(score_change) )
colnames(countexp.Seurat@meta.data)
#
p1 <- FeaturePlot(countexp.Seurat,features = "Pyruvate metabolism")
p2 <- VlnPlot(countexp.Seurat,features = "Pyruvate metabolism")
p1 + p2
p1

pathways <- res@assays$METABOLISM$score
head(rownames(pathways))
# [1] "Terpenoid backbone biosynthesis"                
# [2] "Caffeine metabolism"                            
# [3] "Neomycin, kanamycin and gentamicin biosynthesis"
# [4] "Metabolism of xenobiotics by cytochrome P450"   
# [5] "Drug metabolism - cytochrome P450"              
# [6] "Drug metabolism - other enzymes" 

# 
DimPlot.metabolism(obj = countexp.Seurat, 
                   pathway = "Citrate cycle (TCA cycle)", 
                   dimention.reduction.type = "umap",  
                   dimention.reduction.run = F, size = 1)

# DotPlot


rownames(pathways)
########
input.pathway <- c("Glycolysis / Gluconeogenesis",
                   "Citrate cycle (TCA cycle)",
                   "Pyruvate metabolism",
                   "Oxidative phosphorylation",
                   'Nicotinate and nicotinamide metabolism')


#input.pathway<-head(rownames(pathways))

pdf(file = "1.pdf",width =6,height = 4.5)
DotPlot.metabolism(obj = res,
                   pathway = input.pathway,
                   phenotype = "Celltype", #
                   norm = "y")
dev.off()
#BoxPlot

pdf(file = "2.pdf",width =8,height = 6)
BoxPlot.metabolism(obj = countexp.Seurat,
                   pathway = input.pathway, 
                   phenotype = "Celltype", #
                   ncol =2)+scale_fill_manual(values =colors2)
dev.off()


genes <- list(c("PTEN","PIK3CA","KRAS","ARID1A","RCA1","WNT5A"))
names(genes) <- 'gene'
afe<-afd
#afe <- AddModuleScore_UCell(uterus,features=genes,name="_score")

#
library(ggpubr)



df1 <- data.frame(
  Cell_ID = countexp.Seurat$orig.ident,  # 
  Gene_Expression = countexp.Seurat$`Pyruvate metabolism`, # 
  Cell_Annotation = countexp.Seurat$predicted.id  # 
)


#
my_comparisons1 <- list(c("FMC_2", "vSMC")) 
my_comparisons2 <- list(c("FMC_2", "FMC_4"))
my_comparisons3 <- list(c("FMC_2", "FMC_7"))
my_comparisons4 <- list(c("vSMC", "FMC_4"))
colors<-colors_meta
#
ggplot(df1,aes(x=Cell_Annotation,y=Gene_Expression,fill=Cell_Annotation))+ 
  geom_violin(color='black',size=1)+#
  theme_classic() +   theme(text = element_text(size=10, colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15), 
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.y = element_text(color = 'black', size = 12),
        axis.line = element_line(size = 1))+   
  theme(legend.position="none") +  scale_fill_manual(values = colors)+ 
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+ #
stat_compare_means(method="wilcox.test",hide.ns = F,  
                   comparisons =c(my_comparisons1,my_comparisons2,my_comparisons3,my_comparisons4), 
                   label="p.signif",                
                   bracket.size=0.8,   
                   size=6)#
ggsave(filename ="glycolysis_gluconeogenesis_pathway_violin.png",height = 5.2,width = 5.5,dpi = 300,bg="white")

input.pathway <- c("Glycolysis / Gluconeogenesis",
                   "Citrate cycle (TCA cycle)",
                   "Pyruvate metabolism",
                   "Oxidative phosphorylation",
                   'Nicotinate and nicotinamide metabolism')


###################
#

genes <- list(c("nc_ratio"))
names(genes) <- 'gene'

markers=c('Creb3l1','Sox9','Klf4','Fosl1')
#
library(ggpubr)
afe<- subset(af.ST, Cell_Anno %in% c('vSMC','FMC_2','FMC_4','FMC_7'))
af.ST$domain
#afe<-af.ST
cell_anno <- afd$Celltype
#gene_expr <-  afe$nc_ratio
gene_expr <- GetAssayData(afd, assay = "SCT", slot = "data")["FOSL1", ]
df1 <- data.frame(
  Cell_ID = names(gene_expr),  # 
  Gene_Expression = gene_expr, # 
  Cell_Annotation = cell_anno  # 
)

#
my_comparisons1 <- list(c("FMC_2", "vSMC")) 
my_comparisons2 <- list(c("FMC_2", "FMC_4")) 
my_comparisons3 <- list(c("FMC_2", "FMC_7")) 
my_comparisons4 <- list(c("FMC_4", "vSMC")) 
#

ggplot(df1, aes(x = Cell_Annotation, y = Gene_Expression, fill = Cell_Annotation)) + 
  # 
  geom_boxplot(
    color = 'black',      # 
    size = 0.8,           # 
    width = 0.6,          # 
    outlier.shape = 16,   # 
    outlier.size = 1.5,   # 
    outlier.color = "black",  # 
    alpha = 0.8           # 
  ) +
 # coord_cartesian(ylim = c(0, 8)) +
  # 
  stat_compare_means(
    method = "t.test",
    hide.ns = F,
    comparisons = c(my_comparisons1,my_comparisons2,my_comparisons3,my_comparisons4),
    label = "p.signif",
    bracket.size = 0.8,
    size = 6,y.position = 0.5 
  ) +
  
  # 
  coord_cartesian(ylim = c(-0.4,0.5)) +  # 
  
  scale_fill_manual(values = c('#00a7fa','#354c85','#31aaaa','#9d9ac6') ) +
  
  theme_classic() 



markers=c('Myh11','Tagln','Acta2','Lum','Dcn','Bgn','Fn1','Col1a1','Col1a2',
          'Col3a1' ,'Sox9','Mmp2','Mmp14','Vcam1','Icam1','Angptl4',
          'Slc43a3','Thbs1','C3','C4b','Csf1',
          'Cx3cl1','Cxcl12','Mif','Lgals3','Spp1' ,#FMC_2
          'Tcf4','Runx1'
)
markers <- c(
  'MYH11', 'TAGLN', 'ACTA2', 'LUM', 'DCN', 'FN1', 'COL1A1', 'COL1A2',
  'COL3A1', 'VCAM1', 'ICAM1', 'SPP1','LDHA',"TCF4"  , "MIF","LGALS3"
  
)

markers=c('CREB3L1','SOX9','KLF4','FOSL1')
markers=c('Creb3l1','Sox9','Klf4','Fosl1')
markers=c('Ldha','Pkm','Hif1a','Ldhb')
library(scCustomize)
pdf(file = "FMC DOT .pdf",width =10,height = 5)
Clustered_DotPlot(seurat_object = afd, 
                  colors_use_exp = c( '#64A2CC',"white",'#E07B56',"#731929"),
                  colors_use_idents=c(  "#00a7fa" ,"#354c85","#31aaaa","#9d9ac6"  ),
                  flip=T,x_lab_rotate=90,
                  features = markers)
dev.off()


