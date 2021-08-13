#!/bin/sh

singularity exec --bind /data/codonoghue/SingleCell /data/codonoghue/SingleCell/sc_batch_de.img Rscript bm_seurat_lr.R
