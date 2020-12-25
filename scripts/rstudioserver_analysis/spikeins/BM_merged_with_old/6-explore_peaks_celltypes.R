# Jake Yeung
# Date of Creation: 2020-11-10
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/6-explore_peaks_celltypes.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

GetCtypeFromStr <- function(x){
  # x <- "maxcount_60x2500.BM_round1_round2_merged_H3K4me1_Bcells.2500.cutoff"
  xsplit <- strsplit(strsplit(x, "_")[[1]][[7]], "\\.")[[1]][[1]]
  return(xsplit)
}

GetMarkFromStr <- function(x){
  # x <- "maxcount_60x2500.BM_round1_round2_merged_H3K4me1_Bcells.2500.cutoff"
  xsplit <- strsplit(x, "_")[[1]][[6]]
  return(xsplit)
}

GetMinlengthFromStr <- function(x){
  # x <- "maxcount_60x2500.BM_round1_round2_merged_H3K4me1_Bcells.2500.cutoff"
  xsplit <- strsplit(x, "\\.")[[1]][[3]]
  return(xsplit)
}

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"

inmain <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds/hiddendomains_outputs_minlength_2500.FromR.maxcount_40_60_R"))
outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/peaks_analysis/size_of_peaks.", Sys.Date(), ".pdf")

pdf(outpdf, file = outpdf)

# jmark <- "H3K27me3"

peaksbed.mark.long <- lapply(jmarks, function(jmark){
  indir.vec <- list.files(inmain, pattern = paste0("*_", jmark, "_*"))
  names(indir.vec) <- indir.vec
  peaksbed.mark.long.mark <- lapply(indir.vec, function(dname){
    ctype <- GetCtypeFromStr(dname)
    jmark <- GetMarkFromStr(dname)
    jminlength <- GetMinlengthFromStr(dname)
    # fname2 <- paste0("BM_round1_round2_merged_", jmark, "_", ctype, ".", jminlength, ".cutoff_analysis.bed")
    fname2 <- paste0("BM_round1_round2_merged_", jmark, "_", ctype, ".", jminlength, ".cutoff_vis.bed")
    
    print(dname)
    print(jmark)
    
    inf.tmp <- file.path(inmain, dname, fname2)
    assertthat::assert_that(file.exists(file.path(inf.tmp)))
    if (endsWith(x = fname2, suffix = "_analysis.bed")){
      dat.bed <- fread(inf.tmp, col.names = c("chromo", "jstart", "jend", "peak"), header = FALSE)
    } else if (endsWith(x = fname2, suffix = "_vis.bed")){
      dat.bed <- fread(inf.tmp, col.names = c("chromo", "jstart", "jend", "peak", "score", "strand", "thickStart", "thickEnd", "itemRgb"), header = TRUE)
    } else {
      print(paste("Warning: fname2 should end with _analysis.bed or _vis.bed:", fname2))
    }
    dat.bed$chromo <- paste("chr", dat.bed$chromo, sep = "")
    dat.bed$mark <- jmark
    dat.bed$cluster <- ctype
    dat.bed <- dat.bed %>%
      rowwise() %>%
      mutate(jlength = jend - jstart)
    return(dat.bed)
  }) %>%
    bind_rows()
}) 


# Plot length of peaks across clusters ------------------------------------

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
mlst <- lapply(jmarks, function(jmark){
  m <- ggplot(peaksbed.mark.long[[jmark]], aes(x = jlength, fill = cluster)) + 
    geom_density(alpha = 0.25) + 
    theme_bw() + 
    scale_x_log10() + 
    facet_wrap(~cluster) + 
    ggtitle(jmark) + 
    scale_fill_manual(values = cbPalette) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})

print(mlst)

# Calculate fraction of genome covered ------------------------------------

peaksbed.mark.sum <- peaksbed.mark.long %>%
  bind_rows() %>%
  group_by(cluster, mark) %>%
  summarise(jlength = sum(jlength))
 
ggplot(peaksbed.mark.sum, aes(x = cluster, y = log2(jlength))) + 
 geom_point() + facet_wrap(~mark) + 
 theme_bw() + 
 theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

dev.off()