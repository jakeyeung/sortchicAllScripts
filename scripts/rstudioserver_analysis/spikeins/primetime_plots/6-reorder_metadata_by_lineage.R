# Jake Yeung
# Date of Creation: 2021-01-02
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/6-reorder_metadata_by_lineage.R
# Reorder meta data by lineage

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


ctypes.arranged <-  c("Eryths", "Bcells", "NKs", "pDCs", "Granulocytes", "DCs", "Basophils", "HSPCs")
ctypes.arranged.k9 <-  c("Eryths", "Bcells", "Granulocytes", "HSPCs")

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# Load old metas ----------------------------------------------------------

indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28")
outdir.meta <- file.path(indir.meta, "rearranged_by_lineage")
dir.create(outdir.meta)

dat.merged.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.meta <- file.path(indir.meta, paste0("metadata_batch_corrected.", jmark, ".2020-12-28.txt"))
  print(inf.meta)
  dat.meta <- fread(inf.meta)
})



# Rerrange ----------------------------------------------------------------

dat.merged.rearranged.lst <- lapply(jmarks, function(jmark){
  dat.meta <- dat.merged.lst[[jmark]]
  if (jmark != "H3K9me3"){
    dat.meta$cluster <- factor(dat.meta$cluster, levels = ctypes.arranged)
  } else {
    dat.meta$cluster <- factor(dat.meta$cluster, levels = ctypes.arranged.k9)
  }
  dat.meta <- dat.meta %>%
    arrange(cluster, cuts_in_peak / spikein_cuts)
  return(dat.meta)
})


# Write new metas ---------------------------------------------------------

for (jmark in jmarks){
  print(jmark)
  jdat <- dat.merged.rearranged.lst[[jmark]]
  outftmp <- file.path(outdir.meta, paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".", Sys.Date(), ".txt"))
  fwrite(jdat, file = outftmp, quote = FALSE, sep = "\t", na = "NA")
}
