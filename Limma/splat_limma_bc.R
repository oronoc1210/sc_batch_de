# Libraries
library(Seurat)
library(limma)
library(ggplot2)
library(cowplot)

# Load data
spl <- readRDS(file="../Data/splatter_6g_3b_complex_seurat_norm.Rds")
spl_data <- GetAssayData(object = spl, slot = "data")
spl_limma_bc <- removeBatchEffect(spl_data, spl$Batch)
spl <- SetAssayData(object = spl, slot = "data", new.data = spl_limma_bc)
saveRDS(spl, file="splatter_6g_3b_complex_seurat_limma_bc_2.Rds")
