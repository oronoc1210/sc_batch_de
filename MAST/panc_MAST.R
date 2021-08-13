library(Seurat)
library(MAST)
library(dplyr)

panc <- readRDS("../Data/pancreas_seurat_norm.Rds")
Idents(panc) <- panc$celltype
panc_mast_markers <- FindAllMarkers(object=panc, latent.vars="tech",
				   test.use="MAST")
top5_markers <- panc_mast_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
write.csv(top5_markers, "panc_mast_de_top5.csv")
write.csv(panc_mast_markers, "panc_mast_de_all.csv")
