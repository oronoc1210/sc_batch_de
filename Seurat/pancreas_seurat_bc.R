# Libraries
library(Seurat)
library(ggplot2)
library(cowplot)

# Load data
pancreas_data <- readRDS(file="../Data/pancreas_v3_files/pancreas_expression_matrix.rds")
metadata <- readRDS(file="../Data/pancreas_v3_files/pancreas_metadata.rds")
pancreas <- CreateSeuratObject(pancreas_data, meta.data=metadata)

# Split, normalize, fvf
pancreas_list <- SplitObject(pancreas,split.by="tech")
pancreas_list <- lapply(X=pancreas_list, FUN=function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method="vst", nfeatures=2000)
})

# features, anchors
features <- SelectIntegrationFeatures(object.list = pancreas_list)
pancreas_anchors <- FindIntegrationAnchors(object.list=pancreas_list, anchor.features=features)

# CCA, scale
pancreas_combined <- IntegrateData(anchorset=pancreas_anchors)
DefaultAssay(pancreas_combined) <- "integrated"
pancreas_combined <- ScaleData(pancreas_combined, verbose=FALSE)

# PCA, elbow plot
pancreas_combined <- RunPCA(pancreas_combined, verbose=FALSE)
png("pancreas_seurat_bc_elbow.png")
ElbowPlot(pancreas_combined)
dev.off()

# Find neighbors, clusters
pancreas_combined <- FindNeighbors(pancreas_combined, reduction="pca", dims=1:20)
pancreas_combined <- FindClusters(pancreas_combined, resolution=0.5)

# UMAP
pancreas_combined <- RunUMAP(pancreas_combined, reduction="pca", dims=1:20)
png("pancreas_seurat_bc_umap_cluster.png")
DimPlot(pancreas_combined, reduction="umap")
dev.off()

# Save rds
saveRDS(pancreas_combined, file="pancreas_seurat_bc.rds")

# UMAP by tech/celltype
png("pancreas_seurat_bc_umap_cell.png", width=1000)
p1 <- DimPlot(pancreas_combined, reduction='umap', group.by="tech")
p2 <- DimPlot(pancreas_combined, reduction='umap', group.by="celltype", label=TRUE)
plot_grid(p1,p2)
dev.off()

