# monocle3

library(monocle3)
library(SeuratWrappers) 
library(dplyr)

cds <- as.cell_data_set(sce)
cds <- estimate_size_factors(cds)
cds <- cluster_cells(cds = cds, k = 5, reduction_method = "UMAP",resolution = 0.0000001)
## cds <- cluster_cells(cds,resolution = 0.0000001)

cds <- learn_graph(cds, use_partition = TRUE)
plot_cells(
  cds, 
  color_cells_by = "seurat_clusters", 
  label_groups_by_cluster=TRUE,
  label_leaves=TRUE, 
  label_branch_points=TRUE,
  graph_label_size=3,group_label_size = 5
)

dir.create("./monocle")
setwd("./monocle")
ggsave("EC_umap_monocle.pdf",width = 7,height = 7)


# root cells
rootcell <- rownames(sce@meta.data[which(sce$seurat_clusters == "0"),]) ## 
cds <- order_cells(cds, reduction_method = "UMAP", root_cells=rootcell)

plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_branch_points=FALSE,
  label_leaves = FALSE,
  # graph_label_size = 3
)

ggsave("umap_monocle.root.pdf",width = 7,height = 7)

# diff analysis
subset_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)  # neighbor_graph default knn
head(subset_pr_test_res)

# expression
library(dplyr)
topgene <- subset_pr_test_res %>% top_n(n=12, morans_I) %>% row.names() %>% as.character()
topgene

plot_genes_in_pseudotime(cds[topgene,], color_cells_by="celltype.group", ncol=2,
                         label_by_short_name = F)
ggsave("top12genes_morans_I.pdf", height = 30, width = 15)
plot_genes_in_pseudotime(cds[topgene,], color_cells_by="celltype", ncol=2,
                         label_by_short_name = F)
ggsave("top4genes_morans_I_cellgroup.pdf", height = 30, width = 15)