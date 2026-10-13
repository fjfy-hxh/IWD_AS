# =============================================================================
# 04 - Differential gene volcano plot and Venn diagram
# -----------------------------------------------------------------------------
# Description : FindMarkers-based differential expression, custom volcano plot with labelled genes, and Venn diagrams of overlapping gene sets.
# Inputs      : A Seurat object and gene sets of interest.
# Outputs     : Volcano plot PDF, Venn diagram.
# Depends on  : No repository function; packages: Seurat, ggplot2, ggrepel, ggvenn.
# Note        : All paths are relative to the repository root. Open the .Rproj
#               file (or set the working directory to the repository root)
#               before running this script.
# =============================================================================


focus_cell='FB C4'
afd$cell_cluster = ifelse(afd$mye.anno==focus_cell, focus_cell,"other_cells")###


Idents(afd)<-afd$Hip

table(af.hip$Cell_Anno)
table(afe$group)
afe<-subset(af.ARH, Cell_Anno %in% c('Neurons_IN'))

af<- subset(af.hip, DG %in% c('niche_2','niche_3'))
afe<- subset(af, Cell_Anno %in% c("Astro"))

af<- subset(af.hip, DG %in% c("DG mo",'DG po','niche_2','niche_3'))
afe<- subset(af, Cell_Anno %in% c("Astro"))
Idents(afe)<-afe$Type

table(af$Cell_Anno)

af<-JoinLayers(af)
##
#markers_EX<-markers
markers <- FindMarkers(afe,ident.1 = "IWD", only.pos = F, assay = "SCT",
                       min.pct = 0.25, logfc.threshold = 0)#

#
markers$sign <- ifelse(markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) > 0.3,
                       rownames(markers),NA)
#
k <- c("Hspa5",'Sfpq','Clk1')
markers <- markers %>%
  mutate(sign = ifelse(rownames(markers) %in% k, rownames(markers), NA))


#
log2FC = 0.1
padj = 0.05
markers$change <- cut(markers$avg_log2FC, 
                      breaks = c(-Inf, -log2FC, log2FC, Inf), 
                      labels = c("down", "ns", "up"))
rownames(markers)[markers$change == "up"]
rownames(markers)[markers$change == "down"]
library(ggplot2)
library(ggrepel) #
markers <- markers %>%
  mutate(Difference = pct.1 - pct.2) %>% 
  rownames_to_column("gene")


#
ggplot(markers, aes(x=Difference, y=avg_log2FC, color = change)) + 
  geom_point(size=1.2) + 
  scale_color_manual('change',labels=c(paste0("down(",table(markers$change)[[1]],')'),
                                       'ns',
                                       paste0("up(",table(markers$change)[[3]],')' )),
                     values=c("#0d1b46", "grey","tomato" ))+
  geom_label_repel(data=subset(markers,  p_val <= 0.05), 
                   aes(label=sign),  #
                   color="black", #
                   segment.colour = "black",#
                   label.padding = 0.2, 
                   segment.size = 0.3,  #
                   size=4,
                   max.overlaps = 20) +  # 
  geom_vline(xintercept = 0,linetype = 2) +
  geom_hline(yintercept = 0,linetype = 2) +
  labs(x="△ Percentage difference",y="Log-Fold change") + 
  theme_bw()
ggsave("volcanoplot_DG-mo Astro.pdf",width = 6,height = 4.5)



af.stromaLEC
focus_cell='CN_3'
af.stromaLEC$cell_cluster = ifelse(af.stromaLEC$CN==focus_cell, focus_cell,"other_cells")

# 
venn_list <- list(
  spital_region = fat,
  #Disturbed_gene =Augur_gene#,
 s_vSMC_C1=scDist_genes
)
#
jacc_colors <- c( "#9d9ac6","#009E73")#,'#BFA6C9'
# 


venn(
  venn_list,
  zcolor = jacc_colors,
  opacity = 0.3,
  box = FALSE,
  ilcs = 1,
  sncs = 1,
  counts=T # 
)
dev.off()
nichenet_genes
intersection <- Reduce(intersect, list(marker_genes_diff, Tumor_marker_gene))#,
intersection



