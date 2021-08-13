library(SingleCellExperiment)
library(Seurat)
library(scPOP)

panc <- readRDS("panc_seurat_limma_bc_umap.Rds")
panc_sce <- as.SingleCellExperiment(panc)
metrics <- run_all_metrics(reducedDim(panc_sce), metadata=colData(panc_sce),
                           batch_key='tech', label1_key='celltype', label2_key='tech',
                           run_name="panc_limma")
metrics
ari_score <- ari(panc$tech, panc$celltype)
ari_score
