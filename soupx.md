# code in this article

## soupx
```
library(SoupX)
library(Seurat)
library(HDF5Array)
library(DropletUtils)
# Perform background contamination removal for each sample
##Parameter 
#toc is an analysis matrix, which is a filtered matrix
#tod is full matrix, i.e. a matrix without any filtering
#rho is the pollution ratio coefficient, which can be set by oneself. If not set, it will be automatically calculated

  toc <- Read10X("outs/filtered_feature_bc_matrix",gene.column=1)
  tod <- Read10X("outs/raw_feature_bc_matrix",gene.column=1) 
  
  toc <- Read10X_h5("outs/filtered_feature_bc_matrix.h5")
  tod <- Read10X_h5("outs/raw_feature_bc_matrix.h5")
  
  tod <- tod[rownames(toc) ,]
  
  all <- toc
  all <- CreateSeuratObject(all)
  all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
  all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
  all.genes <- rownames(all)
  all <- ScaleData(all, features = all.genes)
  all <- RunPCA(all, features = VariableFeatures(all), npcs = 30, verbose = F)
  all <- FindNeighbors(all, dims = 1:20)
  all <- FindClusters(all, resolution = 0.5)
  all <- RunUMAP(all, dims = 1:20)
  
  matx <- all@meta.data
  
  sc = SoupChannel(tod, toc)
  sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
  
  if (is.null(rho)) {
    tryCatch(
      {sc = autoEstCont(sc)}, 
      error=function(e) {
        # If error,please rho为0.2
        print("autoEstCont Error !")
        sc = setContaminationFraction(sc, 0.2)} 
    )
  }else{
    
    sc = setContaminationFraction(sc, rho)
  }
  # adjust
  out = adjustCounts(sc)
  PBMC3 = out
  # save
  saveRDS(sc,"./PBMC3.rds")
  # save as，10X format
  DropletUtils:::write10xCounts("../soupX/PBMC3", PBMC3,version="3")
```