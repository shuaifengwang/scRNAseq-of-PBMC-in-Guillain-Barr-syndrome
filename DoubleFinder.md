# Remove double cells

## DoubletFinder

```
## remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

library(dplyr)
library(cowplot)
library(Seurat)
library(ggplot2)
library(psych)
library(qgraph)
library(igraph)
library(Matrix)
library(SeuratWrappers)
library(pryr)

## create folder
dir.create('../Singlet')


Find_doublet <- function(data){
  # find best pk value, if paramSweep_v3 error, please set paramSweep
  sweep.res.list <- paramSweep(data, PCs = 1:20, sct = FALSE) 
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  
  homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)   Doubletrate <- 0.05
  nExp_poi <- round(Doubletrate*ncol(data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  data <- doubletFinder(data, PCs = 1:20, pN = 0.25, pK = p, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  data
}

## Dimensionality reduction clustering
PBMC3 <- Read10X("./soupX/PBMC3/",gene.column=1)
PBMC3 = CreateSeuratObject(PBMC3,
                                project = 'PBMC3',
                                min.cells = 5,
                                min.features = 300)
PBMC3 <- NormalizeData(PBMC3, normalization.method = "LogNormalize", scale.factor = 10000)
PBMC3 <- FindVariableFeatures(PBMC3, selection.method = "vst", nfeatures = 3000)
PBMC3.genes <- rownames(PBMC3)
PBMC3 <- ScaleData(PBMC3, features = PBMC3.genes)
PBMC3 <- RunPCA(PBMC3, features = VariableFeatures(PBMC3), npcs = 30, verbose = F)
PBMC3 <- FindNeighbors(PBMC3, dims = 1:20)
PBMC3 <- FindClusters(PBMC3, resolution = 0.5)

PBMC3 <- RunUMAP(PBMC3, dims = 1:20)
DimPlot(PBMC3, reduction = "umap")+NoLegend()

PBMC3 <- RunTSNE(PBMC3, dims = 1:20, n.neighbors = 30L, min.dist = 0.3)
DimPlot(PBMC3, reduction = "tsne")+NoLegend()


gc()


PBMC3 <- Find_doublet(PBMC3)

# View the position of the determined double cell in the UMAP map
DimPlot(PBMC3, reduction = "umap", group.by = "doublet_info")
ggsave('./Singlet/PBMC3.pdf', width = 6,height = 6)

# save rds
saveRDS(PBMC3, './Singlet/PBMC3.rds')
```