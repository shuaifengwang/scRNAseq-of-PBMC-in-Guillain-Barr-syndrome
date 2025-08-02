# figure

## figure2a
feats <- c("nFeature_RNA", "nCount_RNA", "percent.MT")
  p1_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) +
    NoLegend()

## figure2b-f
DimPlot(pbmc, reduction = "umap", group.by = "doublet_info")

## figure2g
DimPlot(seuratObj,reduction = "harmony",label=T ,group.by = "Sample")

## figure2h
DimPlot(sce, reduction = "tsne", label=T, group.by = "seurat_clusters")
###################
## figure3a
DimPlot(sce, reduction = "umap", group.by ="cell_type", label = TRUE, pt.size = 1.2)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

## figure3b
DotPlot(sce2, features = feature) + RotatedAxis() + 
        scale_colour_gradientn(
        colours = c("blue", "red"),  # 定义颜色渐变（低→高：蓝→红）
        name = "Expression"           # 色标名称（可选）
        )
## figure3c
library(patchwork)

bar.df=sce2@meta.data
active.ident = data.frame(sce2@active.ident)
colnames(active.ident) <- c("maintype")
bar.df = cbind(bar.df, active.ident$maintype)
text.df=as.data.frame(table(bar.df$Sample))
p1=bar.df%>%ggplot(aes(x=Sample))+geom_bar(aes(fill=active.ident$maintype))+
  scale_x_discrete("")+
  scale_y_continuous("cell number",expand = c(0.02,0))+
  geom_text(data = text.df,aes(x=Var1,y=Freq,label=Freq),size=4)+
  theme_classic()+
  theme(
    legend.position = "none"
  )

p2=bar.df%>%ggplot(aes(x=Sample))+geom_bar(aes(fill=active.ident$maintype),position = "fill")+
  scale_x_discrete("")+
  scale_y_continuous("cell ratio",expand = c(0.02,0))+
  theme_classic()

p1+p2

## figure3d
DimPlot(sce, reduction = "umap", group.by ="cell_type", split.by = 'group', label = TRUE, pt.size = 1.2)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

figure3e

#################
## figure4a,4c
DimPlot(T, reduction = "umap", group.by ="cell_type", label = TRUE, pt.size = 1.2)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

## figure4b,4d
colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
  "#A65628", "#F781BF", "#999999", "#6BAED6", "#3182BD", "#08519C",
  "#31A354", "#74C476", "#A1D99B", "#FFB579", "#FF851B", "#E6550D"
)
T$seurat_annotation = T$cell_type
ggplot(T@meta.data,aes(Group, fill = seurat_annotation))+ # Group
  geom_bar(position = "fill", alpha = 0.9)+
  scale_fill_brewer(palette = "coul")+ # scale_fill_manual(values = colors)
  theme_classic()
ggsave("cellraito.pdf", width = 6,height = 6)

## figure4e-h
detailed in monocle3
######################
## figure5a
DimPlot(Mac, reduction = "umap", group.by ="cell_type", label = TRUE, pt.size = 1.2)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

## figure5b
DoHeatmap(subset(Mac, downsample = 100), features = features, size = 3)

## figure5c
colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
  "#A65628", "#F781BF", "#999999", "#6BAED6", "#3182BD", "#08519C",
  "#31A354", "#74C476", "#A1D99B", "#FFB579", "#FF851B", "#E6550D"
Mac$seurat_annotation = Mac$cell_type
ggplot(Mac@meta.data,aes(Group, fill = seurat_annotation))+ # Group
  geom_bar(position = "fill", alpha = 0.9)+
  scale_fill_brewer(palette = "coul")+ # scale_fill_manual(values = colors)
  theme_classic()
ggsave("cellraito.pdf", width = 6,height = 6)

## figure5d,5f
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)

diff<-read.csv("degenes_unfiled.csv", header= T)
gene.df <- bitr(diff$"Gene.Symbol",fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
gene <- gene.df$ENTREZID
head(gene)

ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

dotplot(ego_BP, showCategory=20)

## figure5e,5g
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
dotplot(kk , showCategory=20)
########################
## figure6a-6j
detailed in cellchat

## figure6k
stem: all of the analysis were used OE Online analysis system
url: https://cloud.oebiotech.com/#/home

## figure6l
gene <- read.csv("kegg.id.csv")
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 1,
                 qvalueCutoff = 1)
library(enrichplot)
kkx <- pairwise_termsim(kk)
p2 <- treeplot(kkx, hclust_method = "average")
ggsave(plot = p1, "KEGG/KEGG.tree.pdf", width = 12, height = 6)

