# Libraries
library(Seurat)
library(limma)

# Load data
bm <- readRDS(file="../Data/nuig_bone_marrow_data/nuig_bone_marrow_seurat_norm.Rds")
bm_data <- GetAssayData(object = bm, slot = "data")
bm_limma_bc <- removeBatchEffect(bm_data, bm$sample)
bm <- SetAssayData(object = bm, slot = "data", new.data = bm_limma_bc)
saveRDS(spl, file="../Data/nuig_bone_marrow_data/bm_seurat_limma_bc.Rds")
