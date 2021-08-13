# Libraries
library(Seurat)
library(dplyr)

# Load data
bm <- readRDS(file="../Data/lugh_data/nuig_bone_marrow_seurat_norm.Rds")
dim(bm)

# Seurat LR
Idents(bm) <- bm$celltype
lr_genes <- FindAllMarkers(bm, test.use='LR', latent.vars='sample')
head(lr_genes)
top5_genes <- lr_genes %>% group_by(cluster) %>% top_n(5, avg_log2FC)
write.csv(top5_genes, "bm_seurat_lr_top5.csv")
write.csv(lr_genes, "bm_seurat_lr_all.csv")
