# Data integration, cluster analysis, cell annotation

## Data integration and cluster analysis
```samples = c('HPBMC1','HPBMC2','PBMC1','PBMC2','PBMC3')
HPBMC1 <- readRDS("./Singlet/HPBMC1.rds")
HPBMC2 <- readRDS("./Singlet/HPBMC2.rds")
PBMC1 <- readRDS("./Singlet/PBMC1.rds")
PBMC2 <- readRDS("./Singlet/PBMC2.rds")
PBMC3 <- readRDS("./Singlet/PBMC3.rds")


sce.all <- merge(x = HPBMC1, 
                 y = c(HPBMC2, PBMC1, PBMC2, PBMC3),
                 add.cell.ids = samples)
names(sce.all@assays$RNA@layers)

sce.all[["RNA"]]$counts 
# Alternate accessor function with the same result
LayerData(sce.all, assay = "RNA", layer = "counts")
sce.all 
sce.all <- JoinLayers(sce.all)

head(sce.all@meta.data, 10)
table(sce.all$orig.ident)

#### add group in meta.data
library(stringr)
phe = sce.all@meta.data
table(phe$orig.ident)
View(phe)

phe$group = str_split(phe$orig.ident, "(?<=[A-Z])(?=\\d)", simplify = TRUE)[[,1]]
phe$group <- gsub("\\d+$", "", phe$orig.ident)
sce.all@meta.data$group <- gsub("\\d+$", "", phe$orig.ident)
saveRDS(sce.all, "sce.all.rds")
write.csv(phe, 'phe.csv')

### QC
if(F){
  dir.create("./1-QC")
  setwd("./1-QC")
  sce.all <- readRDS("../sce.all.rds")
  mito_genes = rownames(sce.all)[grep("^MT-", rownames(sce.all),ignore.case = T)]
  print(mito_genes) #mit_genes
  sce.all[["percent.MT"]] <- PercentageFeatureSet(sce.all, features = mito_genes)
  # # ribo_ratio
  ribo_genes=rownames(sce.all)[grep("^RP[SL]", rownames(sce.all),ignore.case = T)]
  print(ribo_genes)
  sce.all=PercentageFeatureSet(sce.all, features = ribo_genes, col.name = "percent.ribo")
  fivenum(sce.all@meta.data$percent_ribo)
  # calculate Hb_ratio
  Hb_genes=rownames(sce.all)[grep("^HB[^(p)]", rownames(sce.all),ignore.case = T)]
  print(Hb_genes)
  sce.all=PercentageFeatureSet(sce.all, features = Hb_genes,col.name = "percent.hb")
  fivenum(sce.all@meta.data$percent.hb)
  head(sce.all@meta.data)
  feats <- c("nFeature_RNA", "nCount_RNA", "percent.MT")
  p1=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) +
    NoLegend()
  p1
  
  ggsave(filename="Vlnplot1.pdf",plot=p1, width = 8, height = 4)
  
  p2=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0)
  p2
  ggsave(filename="Scatterplot.pdf",plot=p2,  width = 4, height = 4)
  
  
  sce.all.filt <- subset(sce.all, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.MT < 5)   
  dim(sce.all)
  
  dim(sce.all.filt)

  
  # visual
  feats <- c("nFeature_RNA", "nCount_RNA", "percent.MT")
  p1_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) +
    NoLegend()
  p1_filtered
  ggsave(filename="Vlnplot1_filtered.pdf",plot=p1_filtered, width = 8, height = 4)
  
  # cell-cycle score
  sce.all.filt = NormalizeData(sce.all.filt)
  s.genes=Seurat::cc.genes.updated.2019$s.genes
  g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
  sce.all.filt=CellCycleScoring(object = sce.all.filt,
                                s.features = s.genes,
                                g2m.features = g2m.genes,
                                set.ident = TRUE)
  p3=VlnPlot(sce.all.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
             ncol = 2, pt.size = 0)
  p3
  ggsave(filename="Vlnplot3_cycle.pdf",plot=p3, width = 6, height = 4)
  
  sce.all.filt@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
    theme_minimal()
  ggsave(filename="cycle_details.pdf" , width = 4, height = 4)
  # S.Score较高的为S期，G2M.Score较高的为G2M期，都比较低的为G1期
  
  dim(sce.all)
  dim(sce.all.filt)
  saveRDS(sce.all.filt, "sce.all_qc.rds")
}
#################
## harmony
if(T){
  dir.create("../2-harmony")
  setwd("../2-harmony/")
  sce <- readRDS("../1-QC/sce.all_qc.rds")
  sce <- NormalizeData(sce)
  sce<- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 4000)
  sce <- ScaleData(sce, features = rownames(sce))
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce@meta.data$Sample = sce@meta.data$orig.ident
  sce <- SetIdent(sce, value = "Sample")
  # sce@meta.data$orig.ident=sce@meta.data$Sample   
  table(sce@meta.data$orig.ident)
  
  seuratObj <- RunHarmony(sce, "orig.ident")
  names(seuratObj@reductions)
  seuratObj <- RunUMAP(seuratObj,  dims = 1:20,
                       reduction = "harmony")
  DimPlot(seuratObj,reduction = "umap",label=F ,group.by = "Sample")
  ggsave(file="harmony_samples.umap.pdf",width = 6,height = 4)
  DimPlot(seuratObj,reduction = "harmony",label=T ,group.by = "Sample")
  sce=seuratObj
  sce <- FindNeighbors(sce, reduction = "harmony",
                       dims = 1:20)
  #for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce=FindClusters(sce, #graph.name = "CCA_snn", 
                   resolution = res, algorithm = 1)
}
sce <- FindClusters(object= sce, resolution=0.5) 
sce <- RunTSNE(sce, reduction = "harmony", dims = 1:20)
DimPlot(sce, reduction = "tsne", label=T, group.by = "seurat_clusters")
ggsave(file="seurat_cluster_tens.pdf",width = 5,height = 4)
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:20)
DimPlot(sce, reduction = "umap",label=T,  group.by = "seurat_clusters")
ggsave(file="seurat_cluster_umap.pdf",width = 5,height = 4)
saveRDS(sce,file = "sce.all_int.rds")
}
```

## cell annotation

```if(F){
  # find markergens
  dir.create("../3-cell");setwd("../3-cell")
  sce <- readRDS("../2-harmony/sce.all_int.rds")
  allmarkers <- FindAllMarkers(sce,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)
  top5_markers <- allmarkers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)
  top10_markers <- allmarkers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  
  # top genes
  top5_markers
  
  # Dot Plot
  DotPlot(sce, features = unique(top5_markers$gene)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) + NoLegend()
  ggsave(file="Dotplot_top5markers.pdf",width = 12, height = 6)
  
  # Heatmap
  DoHeatmap(sce, features = unique(top5_markers$gene)) +
    theme(axis.text.y = ggplot2::element_text(size = 8)) + NoLegend()
  ggsave(file="DoHeatmap_top5markers.pdf",width = 15,height = 10)
  
  # save
  save(allmarkers, top5_markers, file = "markers_integrated.Rdata")
  write.csv(top5_markers, "top5_markers.csv")
}
###############
library(easybio)
head(pbmc.markers)

head(allmarkers)
(marker<-matchCellMarker2(marker=allmarkers,n=50,spc='Human')[,head(.SD,2),by=cluster])
plotPossibleCell(marker)

library(mLLMCelltype)
cache_dir <- "./mllmcelltype_cache"
dir.create(cache_dir, showWarnings =FALSE, recursive =TRUE)

pbmc_markers = allmarkers

single_model_results <- annotate_cell_types(
  input = allmarkers,
  tissue_name ="Single cell  sequencing of PBMC", # 
  model = "deepseek-reasoner",  # 指定单个模型
  api_key = "sk-cb2154221c7b4d7aae6776d02ea28559",  # 直接提供API密钥
  top_gene_count = 10
)

print(single_model_results)
###########
## Artificial correction
## Manual 
anno = as.data.frame(single_model_results)
write.csv(anno, 'anno.csv', row.names = F)

anno = read.csv('anno_top10.csv', header = T)

new.cluster.ids <- anno$celltype
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
sce$cell_type = sce@active.ident
##########
# visual
##########
DimPlot(sce, reduction = "umap", group.by ="cell_type",label = TRUE, pt.size = 1.2) +
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()

ggsave(file="seurat_anno_umap.pdf",width = 6,height = 6)
```