# Libraries
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(cowplot)

# Load data
bm <- readRDS(file="nuig_bone_marrow.Rds")
bm_s <- as.Seurat(bm, counts="counts", data=NULL)

# Split, normalize, fvf
bm_s_list <- SplitObject(bm_s,split.by="sample")
bm_s_list <- lapply(X=bm_s_list, FUN=function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method="vst", nfeatures=2000)
})

# features, anchors
features <- SelectIntegrationFeatures(object.list = bm_s_list)
bm_anchors <- FindIntegrationAnchors(object.list=bm_s_list, anchor.features=features)

# CCA, scale
bm_combined <- IntegrateData(anchorset=bm_anchors)
DefaultAssay(bm_combined) <- "integrated"
bm_combined <- ScaleData(bm_combined, verbose=FALSE)

# PCA, elbow plot
bm_combined <- RunPCA(bm_combined, verbose=FALSE)
png("bm_seurat_bc_elbow.png")
ElbowPlot(bm_combined)
dev.off()

# Find neighbors, clusters
bm_combined <- FindNeighbors(bm_combined, reduction="pca", dims=1:20)
bm_combined <- FindClusters(bm_combined, resolution=0.5)

# UMAP
bm_combined <- RunUMAP(bm_combined, reduction="pca", dims=1:20)
png("bm_seurat_bc_umap_cluster.png")
DimPlot(bm_combined, reduction="umap")
dev.off()

# Save rds
saveRDS(bm_combined, file="bm_seurat_bc.rds")

# UMAP by tech/celltype
png("bm_seurat_bc_umap_cell.png", width=1000, height=1000)
p1 <- DimPlot(bm_combined, reduction='umap', group.by="sample")
p2 <- DimPlot(bm_combined, reduction='umap', group.by="patient")
p3 <- DimPlot(bm_combined, reduction='umap', group.by="timepoint")
p4 <- DimPlot(bm_combined, reduction='umap', group.by="celltype", label=TRUE)
plot_grid(p1,p2,p3,p4, ncol=2)
dev.off()

