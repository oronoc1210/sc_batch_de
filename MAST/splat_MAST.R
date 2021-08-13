library(Seurat)
library(MAST)
library(dplyr)

spl <- readRDS("../Data/splatter_6g_3b_complex_seurat_norm.Rds")
Idents(spl) <- spl$Group
spl_mast_markers <- FindAllMarkers(object=spl, latent.vars="Batch",
				   test.use="MAST")
top5_markers <- spl_mast_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
write.csv(top5_markers, "splat_mast_de_top5.csv")
write.csv(spl_mast_markers, "splat_mast_de_all.csv")
