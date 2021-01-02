# Jake Yeung
# Date of Creation: 2020-12-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/2b-write_K562_metadata_to_output.R
# Write metadata K562 to output

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
outdir.k562 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/K562_metadata")
dir.create(outdir.k562)

# K562 G1 filt  -----------------------------------------------------------

indir.g1filt <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2")

# keep good cells
dat.metas.g1 <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(indir.g1filt, paste0("K562_QC_plots.2020-11-30.add_fracnonzeros.", jmark, ".good_cells.txt"))
  dat.meta <- fread(inf.tmp)
  return(dat.meta)
}) %>%
  bind_rows()

# Ncells data  ------------------------------------------------------------



jdate <- "2020-07-24"
infrds <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins/GenomeWideFits.", jdate, ".logncells.rds")
jsub.sum <- readRDS(infrds)

# jsub.sum <- readRDS(inf.meta)
dat.metas.ncells <- jsub.sum %>%
  rowwise() %>%
  mutate(ncells.lin = ncells,
         ncells.log = log(ncells))
 


# Write outputs -----------------------------------------------------------


# K562 g1 filts

for (jmark in jmarks){
  
  fname.g1filt <- paste0("K562_G1analysis_metadata.", jmark, ".", Sys.Date(), ".txt")
  outf.g1filt <- file.path(outdir.k562, fname.g1filt)
  fwrite(subset(dat.metas.g1, mark == jmark), file = outf.g1filt, sep = "\t", quote = FALSE, na = "NA")
  
  if (jmark %in% c("H3K4me3", "H3K27me3")){
    fname.ncells <- paste0("K562_ncells_metadata.", jmark, ".", Sys.Date(), ".txt")
    outf.ncells <- file.path(outdir.k562, fname.ncells)
    fwrite(subset(dat.metas.ncells, mark == jmark), file = outf.ncells, sep = "\t", quote = FALSE, na = "NA")
  }
}






