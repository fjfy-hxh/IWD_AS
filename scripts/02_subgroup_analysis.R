# =============================================================================
# 02 - Subgroup re-clustering and annotation
# -----------------------------------------------------------------------------
# Description : Subset a cell group of interest, re-run normalization, PCA, Harmony integration, UMAP and clustering, then annotate the subclusters and compare groups.
# Inputs      : A Seurat object together with the subgroup defined in script 01.
# Outputs     : Subgroup UMAPs, violin plots of signature genes, cluster tables.
# Depends on  : ./R/vnplot.R, ./R/ks_VlnExp.R
# Note        : All paths are relative to the repository root. Open the .Rproj
#               file (or set the working directory to the repository root)
#               before running this script.
# =============================================================================



###############
#group diff
source(file = "./R/vnplot.R")
gene_sig <- 'Pomc'
#cell diff

comparisons <- list()
comp=combn(names(table(af$Type)),2)
af$Type
#names(table(af$cellType))
for(j in 1:ncol(comp)){comparisons[[j]]<-comp[,j]}
afvp(afe,gene_signature = gene_sig, file_name = "Slc2a1", test_sign = comparisons,pta=0,
     cols=colors6,label="p.signif",group="Type",widplot=4,heiplot=4,ak=1)

af$fat_acid
######Second method
source("./R/ks_VlnExp.R")


my_comparisons3 <- list(c("IWD", "CWD"))

ks_VlnExp(object = afd, group="Type",group_order=c("IWD",'CWD'),
          features='Acly',cols=c(  "indianred","royalblue"),
          comparisons=c(list(c("IWD", "CWD"))
          ))



table(afd$Hip)

afd[["RNA"]] <- as(afd[["RNA"]], "Assay")
table(af1$Celltype)
af<- subset(af, Celltype %in% c(  'FBs like' ))


#DefaultAssay(afd) <- "RNA"

afd <- NormalizeData(afd, normalization.method = "LogNormalize", scale.factor = 10000)
afd <- FindVariableFeatures(afd, selection.method = "vst", nfeatures = 2000)
scale.genes<-VariableFeatures(afd)
afd <-ScaleData(afd,features = scale.genes)
afd <- Seurat::RunPCA(afd, features = VariableFeatures(object = afd))
pct <- afd[["pca"]]@stdev / sum(afd[["pca"]]@stdev) * 100 ; cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
ElbowPlot(afd)$data %>% ggplot() +
  geom_point(aes(x = dims,y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred") +
  theme_bw() + labs(title = "Elbow plot: quantitative approach")
afd <- Seurat::RunPCA(afd,npcs = 13)

afd <- RunHarmony(afd, group.by.vars = "Type")####


afd=FindNeighbors(afd, dims = 1:13, reduction = "harmony")
afd <- RunUMAP(afd, reduction = "harmony", dims = 1:13)

afd <- FindClusters(afd, resolution = 0.7)

table(afd@active.ident)
afd$sub.cluster<-Idents(afd)



afgenes=c('Tnf','Il1b','Il6','Spp1','Ly6C','Fcgr4','Itgam','S100a9','S100a8',
          'Trem2','Apoe','Fabp5','Lyve1','Folr2','Nampt')

afgenes=c('ACTA2','MYH11','LDHA','DCN','LUM','COL1A1','SPP1','TCF4','MIF','LGALS3')

DotPlot(afd, features = afgenes)+theme_bw()+#
  
  theme(panel.grid = element_blank(),
        
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))+#
  
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#64A2CC',"white",'#E07B56',"#731929"))+#
  # c('#64A2CC',"white",'#E07B56',"#731929")
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+coord_flip()


ann.ids <- c("FMC_2",  #0
             "FMC_2",  #1
             "FMC_7",  #2
             "FMC_2", #3
             "FMC_4", #4
             "vSMC",#5
             "FMC_4", 
             "FMC_7"#7
)
Idents(afd)

afidens=mapvalues(Idents(afd), from = levels(Idents(afd)), to = ann.ids)
Idents(afd)=afidens
afd$sub.cluster=Idents(afd)
afd$seurat_clusters
pdf(file = "subgroup_anno_UMAP.pdf",width =7,height = 5.5)
DimPlot(afd, reduction = "umap", label = F, label.size = 2.5,
       group.by = "sub.cluster",
        pt.size = 1,cols= colors2)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
table(Idents(af.all))

Idents(afd)<-afd$sub.cluster
################
Idents(af, cells = colnames(afd)) <- Idents(afd)
af$Celltype<-Idents(af)

af<-af_subset
af$celltype
table(af$celltype)
 ####
#
af$afd_subcluster <- Idents(af)
Idents(af1)<-af1$afd_subcluster
table(Idents(af1))

af1<-af

table(Idents(af))

af.bei<- subset(afd, subset = sub.cluster != "7")
table(afd$ec)
#Idents(af.1)<-"Proximal tubule cell"
# 
af.ST
af<- subset(af1, idents = selected_idents)

af1$afd_subcluster<-Idents(af1)
Idents(af.bei)<-'FBs like'
#######
Idents(af, cells = colnames(af.bei)) <- Idents(af.bei)


###
new_levels <- c(
  "Inf EC", "EndMT EC", "Angio EC", "Res EC", 
  "T cell", "B cell", "Res Mac", "LAMs", 
  "Inf Mac", "vSMC", 
  "FMC_1", "FMC_2", "FMC_3", "FMC_4", "FMC_5", "FMC_6", "FMC_7"
)

# 
af.ST$Celltype <- factor(af.ST$Celltype, levels = new_levels)

