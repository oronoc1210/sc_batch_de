# Libraries
library(Seurat)
library(dplyr)

# Load data
panc <- readRDS(file="../Data/pancreas_seurat_norm.Rds")
dim(panc)

# Seurat LR
Idents(panc) <- panc$celltype
lr_genes <- FindAllMarkers(panc, test.use='LR', latent.vars='tech')
head(lr_genes)
top5_genes <- lr_genes %>% group_by(cluster) %>% top_n(5, avg_log2FC)
write.csv(top5_genes, "panc_seurat_lr_top5.csv")
write.csv(lr_genes, "panc_seurat_lr_all.csv")
