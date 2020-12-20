# Jake Yeung
# Date of Creation: 2020-11-18
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/11-create_metadata_for_bonemarrow.R
# Create meta data containing spikeincuts, totalcuts, cluster, and plate for each cell

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# Load cluster data -------------------------------------------------------

dat.meta <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM/clusterfilt.2020-11-04/cell_cluster_table.old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.txt"))
  dat.meta <- fread(inf)
  dat.meta$mark <- jmark
  return(dat.meta)
}) %>%
  bind_rows()


# Load spikein and totalcuts data -----------------------------------------



jchromos <- paste(c(seq(19), "X", "Y"), sep = "")
spikeinchromo <- "J02459.1"
indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds/counts_in_peaks_vs_nonpeaks_vs_blacklist.BM-AllMerged3")

dat.spikeins <- lapply(jmarks, function(jmark){
  print(jmark)
  ctypes <- unique(sapply(list.files(indir, pattern = paste0(".*", jmark, ".*genome.csv")), function(x){
    # BM_round1_round2_merged_H3K4me1_pDCs.HiddenDomains.Eryths_2500_H3K4me1.cuts_in_genome.csv
    xsplit <- strsplit(strsplit(x, "\\.")[[1]][[1]], "_")[[1]][[6]]
    return(xsplit)
  }))
  names(ctypes) <- ctypes
  print(ctypes)
  
  dat.merge.annot.lst.bymark <- lapply(ctypes, function(jctype){
    print(jctype)
    
    # paste0("BM_round1_round2_merged_H3K9me3_HSPCs.HiddenDomains.Lymphoid_2500_H3K9me3.cuts_in_peaks.csv")
    inf1 <- list.files(indir, pattern = paste0("BM_round1_round2_merged_", jmark, "_", jctype, ".HiddenDomains.", jctype, ".*.cuts_in_peaks.csv"), full.names = TRUE)
    inf2 <- list.files(indir, pattern = paste0("BM_round1_round2_merged_", jmark, "_", jctype, ".HiddenDomains.", jctype, ".*.cuts_in_genome.csv"), full.names = TRUE)
    
    # jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
    dat1 <- fread(inf1) %>%
      group_by(samp) %>%
      summarise(cuts_in_peak = sum(cuts_in_peak))
    
    dat2 <- fread(inf2) %>%
      filter(chromosome %in% jchromos) %>%
      group_by(samp) %>%
      summarise(cuts_total = sum(cuts_in_chromo))
    
    spikeins <- fread(inf2) %>%
      filter(chromosome == spikeinchromo) %>%
      group_by(samp) %>%
      summarise(spikein_cuts = sum(cuts_in_blacklist))
    
    # NA spikeincuts if round1, nice! 
    dat.merge <- dat1 %>%
      left_join(., dat2, by = "samp")  %>%
      left_join(., spikeins, by = "samp") %>%
      ungroup() %>%
      mutate(ctype = jctype)
    return(dat.merge)
  }) %>%
    bind_rows()
})  %>%
  bind_rows()



# Get plate information ---------------------------------------------------


dat.annot.merge <- left_join(dat.meta, dat.spikeins, by = c("cell" = "samp")) %>%
  rowwise() %>%
  mutate(plate.orig = ClipLast(cell, jsep = "_"))

dat.annots.merge2 <- dat.annot.merge %>%
  mutate(Cluster = ifelse(cluster == "HSPCs", "aHSPCs", cluster)) %>%
  rowwise() %>%
  mutate(batch = IsRound1(cell, mark = mark)) %>%
  mutate(plate = ifelse(batch == "Round2", plate.orig, "Round1"))



# Remove duplicate cells maybe a grep bug ---------------------------------

ncells.unique <- length(unique(dat.annots.merge2$cell))
ncells.total <- length(dat.annots.merge2$cell)

cells.unique <- sort(unique(dat.annots.merge2$cell))

dat.annots.merge2.dupfilt <- subset(dat.annots.merge2, ctype == cluster)

cells.unique.check  <- sort(dat.annots.merge2.dupfilt$cell)

assertthat::assert_that(identical(cells.unique, cells.unique.check))


# Resplit by mark  --------------------------------------------------------

dat.annots.merge2.lst <- split(x = dat.annots.merge2.dupfilt, f = dat.annots.merge2.dupfilt$mark)


# Save to output ----------------------------------------------------------

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins"

for (jmark in jmarks){
  print(jmark)
  dat.tmp <- subset(dat.annots.merge2.lst[[jmark]], select = -c(experi, rowcoord, colcoord))
  fname <- paste0("cell_cluster_table_with_spikeins.", jmark, ".", Sys.Date(), ".dupfilt.txt")
  outf <- file.path(outdir, fname)
  fwrite(dat.tmp, file = outf, quote = FALSE, sep = "\t", na = "NA")
}
