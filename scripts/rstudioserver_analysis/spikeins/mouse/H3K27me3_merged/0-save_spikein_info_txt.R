# Jake Yeung
# Date of Creation: 2020-10-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse/H3K27me3_merged/0-save_spikein_info_txt.R
# Needed for glmpca and downstream analysis

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)





AnnotateSortFromLayout <- function(plate, rowcoord, colcoord){
  assertthat::assert_that(is.numeric(plate))
  assertthat::assert_that(is.numeric(rowcoord))
  assertthat::assert_that(is.numeric(colcoord))
  if (rowcoord >= 1 & rowcoord <= 8 & colcoord == 1){
    print("Empty cell")
    warning("Empty cell")
    ctype <- "Empty"
    return(ctype)
  }
  if (plate >= 1 & plate <= 7){
    if (colcoord >= 1 & colcoord <= 11){
      ctype <- "Unenriched"
    } else if (colcoord >= 12 & colcoord <= 18){
      ctype <- "LinNeg"
    } else if (colcoord >= 19 & colcoord <= 24){
      ctype <- "LSK"
    }
    
  } else if (plate >= 8 & plate <= 13){
    if (colcoord >= 1 & colcoord <= 12){
      ctype <- "Unenriched"
    } else if (colcoord >= 13 & colcoord <= 24){
      ctype <- "LinNeg"
    }
  } else {
    print(paste("Unknown plate:", plate))
    ctype <- "Unknown"
  }
  return(ctype)
}




jmarks <- c("H3K27me3" = "H3K27me3")
jmark <- jmarks[[1]]

inf.chromo <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merged_across_runs/countTablesAndRZr1only_ByChromo.NewFilters/PZ-ChIC_H3K27me3_merged.VAN5046_VAN5230.sorted.tagged.countTable.ByChromo.WithSpikeIns.NoChromo.csv"
inf.bin <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merged_across_runs/countTablesAndRZr1only_ByChromo.NewFilters/PZ-ChIC_H3K27me3_merged.VAN5046_VAN5230.sorted.tagged.countTable.binsize_50000.csv"
inf.rz <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merged_across_runs/countTablesAndRZr1only_ByChromo.NewFilters/PZ-ChIC_H3K27me3_merged.VAN5046_VAN5230.sorted.tagged.countTable.RZ.csv"



# Get good cells ----------------------------------------------------------

jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")

dat.rz <- ReadLH.SummarizeTA(inf.rz)
# dat.chromo <- ReadChrReads(inf.chromo)
dat.chromo <- GetChromoCounts(inf.chromo, spikeinchromo = jspikeinchromo, chromos.keep = jchromos) %>%
  filter(chromo == "1")
dat.bin <- ReadMatSlideWinFormat(inf.bin)

dat.rzchromo <- left_join(dat.rz, dat.chromo)

dat.rzchromo <- dat.rzchromo %>%
  rowwise() %>%
  mutate(rowcoord = AddPlateCoordinates(samp)$rowcoord,
         colcoord = AddPlateCoordinates(samp)$colcoord,
         is.empty = rowcoord <= 8 & colcoord == 1,
         plate = as.numeric(strsplit(experi, "-")[[1]][[6]]),
         stype = AnnotateSortFromLayout(plate, rowcoord, colcoord),
         mark = jmark)

dat.rz.merged <- dat.rzchromo
# dat.rz.merged$mark <- jmark

dat.out.lst <- list()
dat.out.lst[[jmark]] <- list(dat.rz = dat.rzchromo, dat.mat = dat.bin)


ggplot(dat.rzchromo, aes(x = log2(chromocounts / spikeincounts), y = TA.frac, color = is.empty)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  coord_cartesian(ylim = c(0, 1))

ggplot(dat.rzchromo, aes(x = log10(chromocounts), y = TA.frac, color = is.empty)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  coord_cartesian(ylim = c(0, 1))


dat.rz.merged <- lapply(dat.out.lst, function(jdat){
  jdat$dat.rz
}) %>%
  bind_rows()


outf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.chromo2spikeinfilt.merged_across_runs/spikeins_dat_H3K27me3_merged.txt"
fwrite(x = dat.rz.merged, file = outf, sep = ",")


