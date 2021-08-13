#!/bin/sh

singularity exec --bind /data/codonoghue/SingleCell /data/codonoghue/SingleCell/sc_batch_de.img Rscript bm_limma_bc.R
singularity exec --bind /data/codonoghue/SingleCell /data/codonoghue/SingleCell/sc_batch_de.img Rscript panc_limma_bc.R
