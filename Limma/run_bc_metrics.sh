library(SingleCellExperiment)
library(Seurat)
library(scPOP)

# Pancreas
# Seurat UMAP
panc <- readRDS("panc_seurat_limma_bc.Rds")
panc <- NormalizeData(panc)
panc <- FindVariableFeatures(panc, selection.method="vst", nfeatures=2000)
panc <- ScaleData(panc, features=rownames(panc))
panc <- RunPCA(panc, features=VariableFeatures(panc))
panc <- RunUMAP(panc, dims=1:13)
saveRDS(panc, file="panc_seurat_limma_bc_umap.Rds")
png("panc_limma_bc_umap.png")
DimPlot(panc, reduction="umap")
dev.off()

# scPOP
panc_sce <- as.SingleCellExperiment(panc)
panc_metrics <- run_all_metrics(reduction=reducedDim(panc_sce), metadata=colData(panc_sce),
                                batch_key="tech", label1_key="celltype", label2_key="tech",
                                run_name="panc_limma")
write.csv(panc_metrics, "panc_limma_scpop.csv")

# bone marrow
# Seurat UMAP
bm <- readRDS("bm_seurat_limma_bc.Rds")
bm <- NormalizeData(bm)
bm <- FindVariableFeatures(bm, selection.method="vst", nfeatures=2000)
bm <- ScaleData(bm, features=rownames(bm))
bm <- RunPCA(bm, features=VariableFeatures(bm))
bm <- RunUMAP(bm, dims=1:15)
saveRDS(bm, file="bm_seruat_limma_bc_umap.Rds")
png("bm_limma_bc_umap.png")
DimPlot(bm, reduction="umap")
dev.off()

#scPOP
bm_sce <- as.SingleCellExperiment(bm)
bm_metrics <- run_all_metrics(reduction=reducedDim(bm_sce), metadata=colData(bm_sce),
                              batch_key="sample", label1_key="celltype", label2_key="sample",
                              run_name="bm_limma")
write.csv(bm_metrics, "bm_limma_scpop.csv")

# Merge metrics
metrics <- rbind(panc_metrics, bm_metrics)
write.csv(metrics, "combined_limma_scpop.csv")
