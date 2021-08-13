library(Seurat)
library(MAST)
library(dplyr)

bm <- readRDS("../Data/lugh_data/nuig_bone_marrow_seurat_norm.Rds")
Idents(bm) <- bm$celltype
bm_mast_markers <- FindAllMarkers(object=bm, latent.vars=c("patient", "timepoint"),
				   test.use="MAST")
head(bm_mast_markers)
top5_markers <- bm_mast_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
write.csv(top5_markers, "bm_mast_de_top5.csv")
write.csv(bm_mast_markers, "bm_mast_de_all.csv")
