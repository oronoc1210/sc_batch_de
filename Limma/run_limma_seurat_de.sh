#!/bin/sh

# Bone Marrow
singularity exec --bind /data/codonoghue/SingleCell \
  /data/codonoghue/SingleCell/sc_batch_de.img Rscript \
  /data/codonoghue/SingleCell/Seurat_de.R celltype \
  bm_seurat_limma_bc.Rds  bm_limma_bc

# Pancreas
singularity exec --bind /data/codonoghue/SingleCell \
  /data/codonoghue/SingleCell/sc_batch_de.img Rscript \
  /data/codonoghue/SingleCell/Seurat_de.R celltype \
  panc_seurat_limma_bc.Rds  panc_limma_bc

# Splatter
singularity exec --bind /data/codonoghue/SingleCell \
  /data/codonoghue/SingleCell/sc_batch_de.img Rscript \
  /data/codonoghue/SingleCell/Seurat_de.R Group \
  splatter_6g_3b_complex_seurat_limma_bc.Rds splat_limma_bc
