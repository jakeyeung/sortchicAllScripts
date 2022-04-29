# Jake Yeung
# Date of Creation: 2022-04-13
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/8b-clean_H3K27me3_eryths_after_keeping.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_k27me3_remove_bad_eryths"
dir.create(outdir)


# Load k27me3 LDAs --------------------------------------------------------



# Load k27me3 UMAPs -------------------------------------------------------


inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata/repressive_cleaned_keep_eryths/metadata_plate_experi_batch.merged_with_old_dynbins.dynamicbins.k27me3.2022-04-13.txt"
dat.meta <- fread(inf.meta)


# Filter out bad cluster in eryths selectively ----------------------------------------



# Save outputs ------------------------------------------------------------





