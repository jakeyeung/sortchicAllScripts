# Jake Yeung
# Date of Creation: 2019-04-29
# File: ~/projects/scchic/scripts/scripts_analysis/make_databank_objects/make_LDA_objects.R
# Make LDA objects for uploading to GEO

rm(list=ls())


# BINS
inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
load(inf, v=T)

# TSS
traj.mix.lst <- GetTrajMixed()


# Save specific objects LDA ---------------------------------------------------

save(tm.result.lst, file = paste0("~/data/scchic/databank_objects/LDA_matrix_outputs.", Sys.Date(), ".RData"))

