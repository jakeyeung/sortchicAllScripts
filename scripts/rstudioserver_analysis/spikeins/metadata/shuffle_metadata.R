# Jake Yeung
# Date of Creation: 2021-02-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/metadata/shuffle_metadata.R
# 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Load metadata -----------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/shuffled_cells"
dir.create(outdir)

ctypes <- c("Eryths", "Bcells", "NKs", "DCs", "Granulocytes", "Basophils", "pDCs" , "HSPCs")

dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  fread(file.path(indir, fname))
})


# Shuffle order of cells --------------------------------------------------

dat.metas.shuffled <- lapply(dat.metas, function(jdat){
  jdatnew <- jdat %>%
    group_by(cluster) %>%
    sample_frac(., size = 1, replace = FALSE) %>%
    mutate(idx = seq(length(cell)))
  jdatnew$cluster <- factor(jdatnew$cluster, levels = ctypes)
  jdatnew <- jdatnew %>%
    arrange(cluster, idx)
  return(jdatnew)
})

jtest <- dat.metas.shuffled$H3K4me1
plot(x = seq(nrow(jtest)), y = as.factor(jtest$jrep),  type = "o")


# write 
for (jmark in jmarks){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.shuffled.", jmark, ".2021-02-19.txt")
  outf <- file.path(outdir, fname)
  fwrite(dat.metas.shuffled[[jmark]], file = outf, sep = "\t")
}





