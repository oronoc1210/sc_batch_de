library(SingleCellExperiment)
library(Seurat)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
argsLen <- length(args)
ident <- args[1]
rds <- args[2]
outbase <- args[3]

if (argsLen<3){
	stop('please provide ident and filename arguments')
}
if (argsLen>4){
	stop('too many arguments')
}
if (argsLen==4){
	sce <- args[4]
	data <- readRDS(file=rds)
	data <- as.Seurat(data)
} else{
	data <- readRDS(file=rds)
}

Idents(data) <- ident
data_markers <- FindAllMarkers(object=data)
top5_markers <- data_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
saveRDS(data_markers, paste(outbase, "_seurat_de.Rds"))
write.csv(top5_markers, paste(outbase, "_seurat_de_top5.csv", sep=""))
write.csv(data_markers, paste(outbase, "_seurat_de_all.csv", sep=""))
