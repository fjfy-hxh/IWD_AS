# =============================================================================
# 03 - Pseudotime trajectory analysis with Monocle3
# -----------------------------------------------------------------------------
# Description : CytoTRACE2 stemness scoring, Monocle3 trajectory learning and ordering, graph_test gene modules, smoothed pseudotime heatmaps and GO term annotation of modules.
# Inputs      : A Seurat object with cell-type annotation.
# Outputs     : CytoTRACE2 / cell-type trajectory plots, Monocle3 cds RDS, pseudotime heatmaps, module GO term text plot, gene trend curves.
# Depends on  : ./R/spatial_distance_heatmap_addflag.R
# Note        : All paths are relative to the repository root. Open the .Rproj
#               file (or set the working directory to the repository root)
#               before running this script.
# =============================================================================

library(Seurat)
library(monocle3)
library(dplyr)
library(ComplexHeatmap) 
library(ClusterGVis)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(Hmisc)
library(dplyr)
library(future.apply)
library(ggplot2)
library(scales)
library(viridis)
library(org.Hs.eg.db) 
library(org.Mm.eg.db) # org.Mm.eg.db
library(CytoTRACE2)

#
##   
ks_run_monocle3 <- function(object,#seurat obj
                            idents,#cellanno celltype
                            redu_Method = "UMAP", #c("UMAP","TSNE")
                            use_partition,
                            learn_graph_control,
                            define_root=F,
                            know_root=F,
                            root_state
){
  
  require(Seurat)
  require(monocle3)
  require(SeuratWrappers)
  
  Idents(object) <- idents
  #creat cds obj
  cds <- as.cell_data_set(object)
  
  #add the cluster and cellannotation information to cds
  partition <- c(rep(1, length(cds@colData@rownames)))
  names(partition) <- cds@colData@rownames
  partition <- as.factor(partition)
  
  if(redu_Method=="UMAP"){
    
    cds@clusters$UMAP$partitions <- partition
    active_cellanno <- object@active.ident
    cds@clusters$UMAP$clusters <- active_cellanno
    
    #
    cds@int_colData@listData$reducedDims$UMAP <-
      object@reductions$umap@cell.embeddings
    
    
    get_earliest_principal_node <- function(cds){
      cell_ids <- dim(cds)[2]
      closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
      closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
      root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
      root_pr_nodes
    }
    
    
    
  }else{
    
    cds@clusters$TSNE$partitions <- partition
    active_cellanno <- object@active.ident
    cds@clusters$TSNE$clusters <- active_cellanno
    
    #
    cds@int_colData@listData$reducedDims$TSNE <-
      object@reductions$tsne@cell.embeddings
    
    
    get_earliest_principal_node <- function(cds){
      cell_ids <- dim(cds)[2]
      closest_vertex <- cds@principal_graph_aux[["TSNE"]]$pr_graph_cell_proj_closest_vertex
      closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
      root_pr_nodes <- igraph::V(principal_graph(cds)[["TSNE"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
      root_pr_nodes
    }
    
  }
  
  #Learn the trajectory graph
  cds <- learn_graph(cds, use_partition, learn_graph_control=learn_graph_control)
  
  #Order the cells in pseudotime
  if(define_root==F){
    
    cds <-  order_cells(cds)
    
  }else{
    
    if(know_root==F){
      
      cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
      
      
    }else{
      
      if(redu_Method=="UMAP"){
        
        cds <-
          order_cells(cds,
                      reduction_method = 'UMAP',
                      root_cells = colnames(cds[, clusters(cds) == root_state]))
        
      }else{
        
        cds <-
          order_cells(cds,
                      reduction_method = 'TSNE',
                      root_cells = colnames(cds[, clusters(cds) == root_state]))
        
      }
      
      
    }
    
    
  }
  
  return(cds)
  
}

af
scRNA<-af
scRNA <- RenameAssays(scRNA, Segmentations = "RNA")
scRNA<-JoinLayers(scRNA)
scRNA<-cytotrace2(scRNA,is_seurat=TRUE,slot_type="counts",
                  species='mouse'#mouse,human
                  )

#scRNA$cytotrace2 <-Idents(scRNA)
table(scRNA$cytotrace2)
scRNA$celltype<-scRNA$Celltype
##run
#cds_data1 <- ks_run_monocle3(object=scRNA,
#                             idents="celltype",
#                             use_partition=F,
#                             learn_graph_control=list(minimal_branch_len=9.5,euclidean_distance_ratio=1),
#                             define_root=T,
#                             know_root=T,
#                             root_state="FBs")


#########
cds_data <- ks_run_monocle3(object=scRNA,
                            idents="Celltype",
                            use_partition=F,
                            learn_graph_control=list(minimal_branch_len=9.5,euclidean_distance_ratio=1),
                            define_root=F)

###
pseudotime$`pseudotime(cds_data)`
pseudotime <- as.data.frame(pseudotime(cds_data))
# 
common_cells <- intersect(rownames(pseudotime), colnames(scRNA))
cat("Number of common cells: ", length(common_cells), "/", ncol(scRNA), "\n")
# 
pseudotime_matched <- pseudotime[common_cells, "pseudotime(cds_data)"]
scRNA_matched <- scRNA[, common_cells]

# 
scRNA_matched$pseudotime <- pseudotime_matched

####################################################
cds<-cds_data
scRNA$CytoTRACE2_Score
#plot
pdf(file = "pseudotime_trajectory.pdf",width =6,height = 5)
plot_cells(cds, label_cell_groups = F, 
           color_cells_by = "CytoTRACE2_Score", #pseudotime
           label_branch_points = F, 
           label_roots =T,
           label_leaves =F,
           graph_label_size = 0, 
           cell_size=1, 
           trajectory_graph_color='black',
           trajectory_graph_segment_size = 1.5)#+
        # scale_color_gradientn(colors = c('#64A2CC',"white",'#E07B56',"#731929"))
dev.off()



pdf(file = "pseudotime_FB_celltype.pdf",width =6,height = 5)
plot_cells(cds, label_cell_groups = F, 
           color_cells_by = "Celltype", 
           label_branch_points = F, 
           label_roots =T,
           label_leaves =F,
           graph_label_size = 0, 
           cell_size=1.5, 
           trajectory_graph_color='black',
           trajectory_graph_segment_size = 1.5)+
  scale_color_manual(values = c( '#cce4b2','#6eb881','#eec88a','#FCCC3A','#f3758a','#f4bbb2','#b9e6fb')) 
dev.off()

#-----------------------------------------------------------------------------
saveRDS(cds, file = "monocle3_cds_object.rds")
#   

#-----------------------------------------------------------------------------
#

colData(cds)
branch_cells <- colData(cds)$Celltype %in% c("FMC_2", "FMC_4",'vSMC')
cds_branch <- cds[, branch_cells]
branch_test <- graph_test(cds_branch, neighbor_graph = "principal_graph", cores = 14)
plot_genes <- row.names(subset(branch_test, q_value < 0.01 & morans_I > 0.2))
plot_genes
aaa<-as.data.frame(pseudotime(cds))
#
selected_cells <- choose_graph_segments(cds)####
selected_cells <- preprocess_cds(selected_cells, method = "PCA")
selected_cells <- reduce_dimension(selected_cells, reduction_method = "UMAP")
selected_cells <- cluster_cells(selected_cells, resolution = 1e-5)  #
selected_cells <- learn_graph(selected_cells)
plot_cells(selected_cells, reduction_method = "UMAP")

modulated_genes <- graph_test(selected_cells, neighbor_graph = "principal_graph", 
                              cores = 14)
plot_genes <- row.names(subset(modulated_genes, q_value < 0.01 & morans_I > 0.2))

plot_genes1 <- c(plot_genes, 
                 c("FAP", "ACTA2",'PDGFRB','CXCL12','IL6',
                   'ELF1','CEBPD','CEBPA','JUND','KLF10',
                   'FOSB','RARG','IRF1','KLF4','EGR1','HES1',
                   'JUN','JUNB'))
graph_gene<-plot_genes


###########
genes <- plot_genes
mycds2<-cds_branch
#
plot_matrix <- exprs(mycds2)[match(genes,
                                   rownames(rowData(mycds2))),
                             order(pseudotime(mycds2))]#

#
plot_matrix <- t(apply(plot_matrix,1,function(x){smooth.spline(x,df=3)$y}))
plot_matrix <- t(apply(plot_matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(plot_matrix) <- genes;
dim(plot_matrix)


cutColumn_Means <- function(data_exp,#
                            cut#
){
  plot_matrix_combin <- list()
  nums <- ncol(data_exp)/cut
  if (nums-round(nums, 0)==0){
    
    for (i in 1:length(seq(1, ncol(data_exp), cut))){
      num <- seq(1, ncol(data_exp), cut)
      A <- as.data.frame(rowMeans(data_exp[,num[i]:(cut+num[i]-1)]))[,1]
      plot_matrix_combin[[i]] <- A
      
    }
    plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
    rownames(plot_matrix_combin) <- rownames(data_exp)
    colnames(plot_matrix_combin) <- seq(1,ncol(plot_matrix_combin),1)
    return(plot_matrix_combin)
    
  }else{
    
    for (i in 1:length(seq(1, ncol(data_exp)-cut, cut))){
      num <- seq(1, ncol(data_exp)-cut, cut)
      A <- as.data.frame(rowMeans(data_exp[,num[i]:(cut+num[i]-1)]))[,1]
      plot_matrix_combin[[i]] <- A
    }
    
    plot_matrix_combin[[length(seq(1, ncol(data_exp)-cut, cut))+1]] <- as.data.frame(rowMeans(data_exp[,(max(num)+cut):ncol(data_exp)]))                       
    plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
    rownames(plot_matrix_combin) <- rownames(data_exp)
    colnames(plot_matrix_combin) <- seq(1,ncol(plot_matrix_combin),1)
    return(plot_matrix_combin)
  }
  
  
}


plot_test <- cutColumn_Means(plot_matrix,cut = 150)
plot_matrix_combin<-plot_test
p1 <- pheatmap::pheatmap(plot_matrix_combin, 
                         useRaster = T,
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=F, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         cutree_rows=4,
                         filename=NA,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                         clustering_callback = callback)


callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}


annotation_row <- data.frame(Cluster=factor(cutree(p1$tree_row, 4)))
row.names(annotation_row) <- rownames(plot_matrix_combin)

rowcolor <- c("#85B22E","#E29827","#922927",'#57C3F3') 
names(rowcolor) <- c("1","2","3","4") # Type colors

annotation_row$Cluster <- as.character(annotation_row$Cluster)
names(rowcolor) <- unique(annotation_row$Cluster)

#
ann_colors <- list(
  Cluster = rowcolor, # 
  Order = colorRampPalette(c('#000075', '#2E00FF', '#9408F7', '#C729D6', '#FA4AB5',
                             '#FF6A95', '#FF8B74', '#FFAC53', '#FFCD32', '#FFFF60'))(100)
)

# 
if (is.null(colnames(plot_matrix_combin))) {
  colnames(plot_matrix_combin) <- paste0("C", 1:ncol(plot_matrix_combin))
}
col_ann <- data.frame(Order = 1:ncol(plot_matrix_combin))
rownames(col_ann) <- colnames(plot_matrix_combin)


p4 <- pheatmap::pheatmap(plot_matrix_combin, 
                         cluster_cols = FALSE, 
                         cluster_rows = TRUE, 
                         show_rownames = TRUE, 
                         show_colnames = FALSE, # 
                         clustering_method = "ward.D2",
                         border_color = NA,
                         fontsize_row = 8,
                         fontsize_main = 22,    # 
                         color = colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                         
                         # 
                         annotation_row = annotation_row,
                         annotation_col = col_ann,
                         annotation_colors = ann_colors,
                         
                         annotation_names_row = FALSE,
                         annotation_names_col = FALSE,
                         clustering_callback = callback,
                         main = "Pseudotime")

source('./R/spatial_distance_heatmap_addflag.R')
add.flag(p4,kept.labels = markers,repel.degree = 0.2)

library(clusterProfiler)
library(org.Mm.eg.db)

module_gene <- as.data.frame(cutree(p4$tree_row, k=4))
colnames(module_gene) <- "Module"
module_gene$gene <- rownames(module_gene)
# 
p_final <- plot_module_text_boxes(
  gene_list = module_gene, 
  species = "mouse", #human  mouse
  rowcolor = rowcolor,
  total_terms = 20
)

print(p_final)

ggsave("pseudotime_GO_term_text.pdf", p_final, 
       width = 6.5, height = 15,
       device = cairo_pdf)

###########

pseudotime <- pseudotime(cds_data) %>% as.data.frame()#
pseudotime$cell <- rownames(pseudotime)
colnames(pseudotime)[1] <- "pseudotime"

celltype<-NULL
# 
umap_coords <- reducedDims(cds_data)$UMAP
celltype <- as.data.frame(umap_coords)
celltype$cell <- rownames(celltype)
###########################
celltype$celltype <- cds_data@colData$Celltype ####celltype

merge<-NULL
merge$celltype
merge <-merge(pseudotime, celltype, by = 'cell')#
#merge <- merge[order(merge$peu), ]#
str(merge)
table(merge$celltype)
merge$celltype<-as.character(merge$celltype)
#View(merge)
colors_meta
##plot1==================================================================================
p=ggplot(merge,aes(pseudotime,fill=celltype, color=celltype)) +  
  geom_density(alpha = 0.5,size=1.2)+#
  scale_fill_manual(values=c("#00a7fa" ,"#354c85", "#31aaaa" ,"#9d9ac6"  ))+  # 
  scale_color_manual(values=c("#00a7fa", "#354c85", "#31aaaa", "#9d9ac6" ))+  
  theme_classic()+ #
  scale_y_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0, 0), limits = c(0, 0.3)) +
  theme(axis.line = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x=element_text(colour='black',size=10))+
  labs(title = "Pseudotime")
p  

mycds_order <- cds
pData(mycds_order)$pseudotime <- pseudotime(mycds_order)#
mycds_order <- mycds_order[,is.finite(pseudotime(mycds_order))]
#
a <- as.numeric(mycds_order$pseudotime)#

table(scRNA$celltype)
mycds_order$pseudotime_bin <-  cut2(as.numeric(mycds_order$pseudotime),
                                    seq(min(a),max(a), by = 0.15)) #
table(mycds_order@colData@listData[["pseudotime_bin"]])#
table(mycds_order@colData@listData[["celltype"]],  #####  
      mycds_order@colData@listData[["pseudotime_bin"]])#

bins <- levels(mycds_order$pseudotime_bin)#
proportion_df <- data.frame()#

for(i in 1:length(bins)){  
  bin = bins[i]  
  meta <- subset(as.data.frame(mycds_order@colData@listData), pseudotime_bin == bin) 
if(nrow(meta) == 0) { #
  df <- data.frame(orig.ident = c("FBs", "my_iCAF"),  ####
                   proportion = c(0,0), 
                   bin_num = c(i, i),  
                   bin= c(bins[i], bins[i])) 
  proportion_df <- rbind(proportion_df, df)  
}      
  
  else{    
  df <- data.frame(table(meta$celltype) / nrow(meta)) #####
  df <- df %>% dplyr::rename(c(orig.ident = Var1, proportion = Freq))  
  df$bin_num <- i  
  df$bin <- bin  
  proportion_df <- rbind(proportion_df, df)  }  
}

#
proportion_cor <- cor.test(  subset(proportion_df[,], 
                                    orig.ident == 'my_iCAF') %>% .$bin_num, ##
                             subset(proportion_df[,],
                                    orig.ident == 'my_iCAF') %>% .$proportion,#
                             method='pearson')

#

p <- ggscatter(  subset(proportion_df, orig.ident == 'my_iCAF'),###
                 color = "#D43D4F",  x = 'bin_num', y = 'proportion',  add ='reg.line',
                                     add.params = list(color = '#5c4b4d', fill = 'lightgray'),
                 conf.int=TRUE,  cor.coef=TRUE,  
                 cor.coeff.args = list(method = 'pearson', label.sep = '\n'), 
                 alpha = 0.8,  
                 size = 2.5,  
                 title = "Propotion of my_iCAF along the pseudotime") + 
                 xlab('pseudotime') + 
                 ylab('proportion')

p
#################

cds@colData$Type
peuGene_exp <- normalized_counts(mycds2, norm_method = "log")
peuGene_exp <- as.data.frame(peuGene_exp)
dim(peuGene_exp)# [1] 11199  6025
pseudotime <- as.data.frame(pseudotime(mycds2))#
colnames(pseudotime) <- "pseudotime"#f
pseudotime <- pseudotime[rowSums(abs(pseudotime)==Inf)==0,,drop=F]
peuGene_exp <- peuGene_exp[,colnames(peuGene_exp) %in% rownames(pseudotime)]
dim(peuGene_exp)# [1] 11199  5951

#
A_cells <- rownames(cds@colData)[cds@colData$Celltype == "FMC_2"]
peuGene_exp_A <- peuGene_exp[, colnames(peuGene_exp) %in% A_cells]
pseudotime_A <- pseudotime[colnames(peuGene_exp_A),]

B_cells <- rownames(cds@colData)[cds@colData$Celltype == "FMC_4"]
peuGene_exp_B <- peuGene_exp[, colnames(peuGene_exp) %in% B_cells]
pseudotime_B <- pseudotime[colnames(peuGene_exp_B),]

C_cells <- rownames(cds@colData)[cds@colData$Celltype == "vSMC"]
peuGene_exp_C <- peuGene_exp[, colnames(peuGene_exp) %in% C_cells]
pseudotime_C <- pseudotime[colnames(peuGene_exp_C),]
#
A_anxal1 <- peuGene_exp_A['Zbtb20',]
A_anxal1 <- as.data.frame(t(A_anxal1))
A_anxal1$pseudotime <- pseudotime_A
A_anxal1$group <- 'FMC_2'

B_anxal1 <- peuGene_exp_B['Zbtb20',]
B_anxal1 <- as.data.frame(t(B_anxal1))
B_anxal1$pseudotime <- pseudotime_B
B_anxal1$group <- 'FMC_4'

C_anxal1 <- peuGene_exp_C['Zbtb20',]
C_anxal1 <- as.data.frame(t(C_anxal1))
C_anxal1$pseudotime <- pseudotime_C
C_anxal1$group <- 'vSMC'

peu_trend <- rbind( A_anxal1,B_anxal1,C_anxal1)

#
ggplot(peu_trend, aes(x=pseudotime, y=Zbtb20, color=group))+ 
  geom_smooth(aes(fill=group),method = "loess")+ #
  xlab('pseudotime') +   ylab('Relative expression') +  ggtitle(' ')+ 
  theme_classic(base_size = 12)+  theme(axis.text = element_text(color = 'black',size = 12),
                                        axis.title = element_text(color = 'black',size = 14))+
  scale_color_manual(name=NULL, values = c("#00a7fa" ,"#354c85" , "#9d9ac6"))+#
  scale_fill_manual(name=NULL, values = c("#00a7fa", "#354c85", "#9d9ac6"))#


mouse_sc <- scRNA[, rownames(mycds_order@colData)]#
mouse_sc$pseudotime <- mycds_order@colData$pseudotime
mouse_sc$pseudotime_bin <- mycds_order@colData$pseudotime_bin


PMN_state_gene <- c('PCDH7','IFITM1')########
mouse_sc <- AddModuleScore(mouse_sc, 
                           features=list('ISC' = PMN_state_gene), 
                           pool = rownames(mouse_sc), k=F, nbin=24,
                           name = c('ISC'))



modules <- select(mouse_sc@meta.data, c(pseudotime_bin, ISC1))#
modules$pseudotime_bin_num <- as.numeric(modules$pseudotime_bin)
features <- c('ISC1')
tmp <- lapply(1:max(modules$pseudotime_bin_num), 
              function(i){  df <- modules %>% subset(pseudotime_bin_num == i) 
              data.frame(    value = as.numeric(sum(df[,'ISC1']) / nrow(df)),  
                             bin_num = i,    feature = features)}) #
plot_df <- Reduce(rbind, tmp)#
plot_df <- na.omit(plot_df)#
#

ggplot(plot_df, aes(bin_num, value)) +
                geom_hline(yintercept = 0, 
                           linetype = 'dashed', 
                           color = 'gray', alpha = 0.75) + 
                geom_point(size = 2, color = '#7DBFA7') +  # 
                geom_smooth(aes(color = I("#3C77AF"), fill = I("#AECDE1")), se = TRUE,size = 2) + # Set the smooth line color to yellow
                xlab('pseudotime') +
                ylab('Module score') + 
                theme_bw() +  
                theme(panel.grid.major = element_blank(),  
                      panel.grid.minor = element_blank(),  
                      axis.ticks.x = element_blank(),     
                      axis.text.x = element_blank()) +
                ggtitle('Gene module along the pseudotime')


