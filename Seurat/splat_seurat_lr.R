# Libraries
library(Seurat)
library(dplyr)

# Load data
spl <- readRDS(file="../Data/splatter_seurat_norm.Rds")
dim(spl)

# Seurat LR
Idents(spl) <- spl$Group
lr_genes <- FindAllMarkers(spl, test.use='LR', latent.vars='Batch')
head(lr_genes)
top5_genes <- lr_genes %>% group_by(cluster) %>% top_n(5, avg_log2FC)
write.csv(top5_genes, "spl_seurat_lr_top5.csv")
write.csv(lr_genes, "spl_seurat_lr_all.csv")
