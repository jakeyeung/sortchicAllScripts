# Jake Yeung
# Date of Creation: 2021-01-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/DE_analysis/6-check_DE_and_high_bins_for_GC_content_and_distance.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

binsize <- 50000

# ExtendBinSize <- function(coord){
#   dat.coord <- GetBedFromCoords(coord, add.chr = FALSE, strip.chr = FALSE)
# }

# Load metas  -------------------------------------------------------------

indir.metas <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage"
dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  inf <- file.path(indir.metas, fname)
  fread(inf)
})



# Load high bins  ---------------------------------------------------------

indir.bins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned"

dat.bins <- lapply(jmarks, function(jmark){
  fname <- paste0("High_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.", jmark, ".2021-01-30.txt")
  inf.bins <- file.path(indir.bins, fname)
  fread(inf.bins)
})

dat.bins <- lapply(dat.bins, function(jdat){
  jdat <- jdat %>%
    rowwise() %>%
    mutate(startExtend = start + 1 - binsize / 2,
           endExtend = end - 1 + binsize / 2)
  jdat$region_coordExtend <- paste(jdat$seqnames, paste(jdat$startExtend, jdat$endExtend, sep = "-"), sep = ":")
  return(jdat)
})

dat.de.bins <- lapply(jmarks, function(jmark){
  fname <- paste0("DE_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.", jmark, ".2021-01-30.txt")
  inf.bins <- file.path(indir.bins, fname)
  fread(inf.bins)
})

dat.de.bins <- lapply(dat.de.bins,function(jdat){
  jdat <- jdat %>%
    rowwise() %>%
    mutate(startExtend = start + 1 - binsize / 2,
           endExtend = end - 1 + binsize / 2)
  jdat$region_coordExtend <- paste(jdat$seqnames, paste(jdat$startExtend, jdat$endExtend, sep = "-"), sep = ":")
  return(jdat)
})



# Load K27me3 HSPC bins  --------------------------------------------------

inf.k27me3 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/tables_H3K27me3_HSPC_lost_bins/H3K27me3_HSPCs_lost_bins_DE_estimates.2021-01-29.txt"
dat.k27me3 <- fread(inf.k27me3)


# Count features within each bins -----------------------------------------

# GCs
load("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_gc_analysis/gcs_genomewide.RData", v=T)

gr.gc.dat.dedup <- gr.gc.dat[!duplicated(gr.gc.dat), ]


# Calculate GCs  ----------------------------------------------------------


dat.bins.annot.long <- lapply(jmarks, function(jmark){
  jdat <- left_join(dat.bins[[jmark]], gr.gc.dat.dedup, by = c("region_coordExtend" = "bname")) 
}) %>%
  bind_rows() %>%
  mutate(type = "HighBins") %>%
  dplyr::rename(bin = region_coordExtend)

dat.de.bins.annot.long <- lapply(jmarks, function(jmark){
  jdat <- left_join(dat.de.bins[[jmark]], gr.gc.dat.dedup, by = c("region_coordExtend" = "bname")) 
}) %>%
  bind_rows() %>%
  mutate(type = "DiffBins") %>%
  dplyr::rename(bin = region_coordExtend)

dat.de.k27me3 <- left_join(dat.k27me3, gr.gc.dat.dedup, by = c("bin" = "bname")) %>%
  bind_rows() %>%
  mutate(mark = "H3K27me3",
         type = "StemCellLostBins") 

dat.bins.merge <- rbind(dat.bins.annot.long, dat.de.bins.annot.long)

jmark.keep <- "H3K27me3"
dat.k27me3.merge <- rbind(subset(dat.bins.annot.long, mark == jmark.keep, select = c(mark, type, bin, gc)), 
                          subset(dat.de.bins.annot.long, mark == jmark.keep, select = c(mark, type, bin, gc)), 
                          subset(dat.de.k27me3, select = c(mark, type, bin, gc)))

ggplot(dat.bins.annot.long, aes(x = mark, y = gc)) + 
  geom_point() + 
  geom_boxplot() + 
  ggtitle("High bins") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.de.bins.annot.long, aes(x = mark, y = gc)) + 
  geom_point() + 
  geom_boxplot() + 
  ggtitle("DE bins") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.bins.merge, aes(x = interaction(type, mark), y = gc)) + 
  geom_point() + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.k27me3.merge, aes(x = interaction(type, mark), y = gc)) + 
  geom_point() + 
  geom_boxplot() + 
  ggtitle("DE bins") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())






