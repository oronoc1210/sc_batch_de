# Libraries
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(cowplot)

# Load data
spl <- readRDS(file="../Data/splatter_6g_3b_complex.Rds")
spl_s <- as.Seurat(spl, counts="counts", data=NULL)

# Split, normalize, fvf
spl_s_list <- SplitObject(spl_s,split.by="Batch")
spl_s_list <- lapply(X=spl_s_list, FUN=function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method="vst", nfeatures=2000)
})

# features, anchors
features <- SelectIntegrationFeatures(object.list = spl_s_list)
spl_anchors <- FindIntegrationAnchors(object.list=spl_s_list, anchor.features=features)

# CCA, scale
spl_combined <- IntegrateData(anchorset=spl_anchors)
DefaultAssay(spl_combined) <- "integrated"
spl_combined <- ScaleData(spl_combined, verbose=FALSE)

# PCA, elbow plot
spl_combined <- RunPCA(spl_combined, verbose=FALSE)
png("spl_x_seurat_bc_elbow.png")
ElbowPlot(spl_combined)
dev.off()

# Find neighbors, clusters
spl_combined <- FindNeighbors(spl_combined, reduction="pca", dims=1:20)
spl_combined <- FindClusters(spl_combined, resolution=0.5)

# UMAP
spl_combined <- RunUMAP(spl_combined, reduction="pca", dims=1:20)
png("spl_x_seurat_bc_umap_cluster.png")
DimPlot(spl_combined, reduction="umap")
dev.off()

# Save rds
saveRDS(spl_combined, file="spl_x_seurat_bc.rds")

# UMAP by tech/celltype
png("spl_x_seurat_bc_umap_cell.png", width=1000)
p1 <- DimPlot(spl_combined, reduction='umap', group.by="Batch")
p2 <- DimPlot(spl_combined, reduction='umap', group.by="Group", label=TRUE)
plot_grid(p1,p2)
dev.off()

