# Jake Yeung
# Date of Creation: 2020-10-29
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_round2/5-K562_cellcycle_from_peaks.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

hubprefix <- "/home/jyeung/hub_oudenaarden"

indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/bams_CellCycleSorted_split_by_cellcycle.NoChrPrefix/countTablesAndRZr1only_ByBed.NewFilters")

# K562_CellCycleSorted_H3K4me3.merged.sorted.tagged.X0_G1.sorted.bam.countTable.bedfile.csv
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
cellcycles <- c("X0_G1", "X1_S", "X2_G2_M"); names(cellcycles) <- cellcycles

mats.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  outs.lst <- lapply(cellcycles, function(cellcycle){
    print(cellcycle)
    bname <- paste0("K562_CellCycleSorted_", jmark, ".merged.sorted.tagged.", cellcycle, ".sorted.bam.countTable.bedfile.csv")
    inf.tmp <- file.path(indir, bname)
    jmat <- ReadMatTSSFormat(inf.tmp, add.coord = TRUE, sort.rnames = TRUE)
    jdat.annot <- data.frame(cell = colnames(jmat), cellcycle = cellcycle, stringsAsFactors = FALSE)
    return(list(jmat = jmat, jdat.annot = jdat.annot))
  })
  mats.merge.lst <- lapply(outs.lst, function(x) x$jmat)
  rows.common <- Reduce(f = intersect, lapply(mats.merge.lst, function(jmat) rownames(jmat)))
  mats.merge.filt.lst <- lapply(mats.merge.lst, function(jmat) jmat[rows.common, ])
  mats.merge <- do.call(cbind, mats.merge.filt.lst)
  dat.annot.lst <- lapply(outs.lst, function(x) x$jdat.annot) %>%
    bind_rows()
  return(list(mats.merge = mats.merge, dat.annot.lst = dat.annot.lst))
})



# Plot slopes -------------------------------------------------------------

# merge across cells

dat.peak.cuts <- lapply(jmarks, function(jmark){
  print(jmark)
  jdat <- data.frame(peakcounts = colSums(mats.lst[[jmark]]$mats.merge), cell = colnames(mats.lst[[jmark]]$mats.merge), stringsAsFactors = FALSE)
}) 


# Load spikeins -----------------------------------------------------------

inf.spikeins.cc.hoescht <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/K562_RData_objs/K562_cellcycle_with_hoescht_4_marks.rds"
# inf.spikeins.cc <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2/dat_spikeins_all.round2.WithCellCycleLabel.rds"
dat.spikeins.cc.hoescht <- readRDS(inf.spikeins.cc.hoescht)

dat.merge.lst <- lapply(jmarks, function(jmark){
  dat.merge <- left_join(dat.spikeins.cc.hoescht[[jmark]], dat.peak.cuts[[jmark]])
})

m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.merge.lst[[jmark]], aes(x = cellcycle.str, y = log2(peakcounts / spikeincounts), fill = cellcycle.str)) + 
    geom_boxplot(alpha = 0.5) + 
    geom_point()  + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})
JFuncs::multiplot(m.lst$H3K4me1, m.lst$H3K4me3, m.lst$H3K27me3, m.lst$H3K9me3, cols = 4)


m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.merge.lst[[jmark]], aes(x = hoesch.scale, fill = cellcycle.str)) + 
    geom_density(alpha = 0.25) + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})
print(m.lst)

JFuncs::multiplot(m.lst$H3K4me1, m.lst$H3K4me3, m.lst$H3K27me3, m.lst$H3K9me3, cols = 4)

m <- ggplot(dat.merge.lst$H3K4me1, aes(x = hoesch.scale, y = log2(peakcounts / spikeincounts), color = cellcycle.str)) + 
  geom_point()  + 
  theme_bw() + 
  ggtitle("H3K4me1") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)



m <- ggplot(dat.merge.lst$H3K4me1, aes(x = hoesch, y = log2(peakcounts / spikeincounts), color = cellcycle.str)) + 
  geom_point()  + 
  theme_bw() + 
  ggtitle("H3K4me1") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)


# Do fits -----------------------------------------------------------------

dat.glmpca.peak.lst <- lapply(dat.merge.lst, function(jsub){
  jsub$ltr <- log2(jsub$peakcounts / jsub$spikeincounts)
  jsub$ltr.wins <- DescTools::Winsorize(jsub$ltr, probs = c(0.01, 0.99))
  return(jsub)
})

jslopes.peak <- lapply(jmarks, function(jmark){
  jsub <- dat.glmpca.peak.lst[[jmark]]
  jfit <- lm(formula = log2(peakcounts / spikeincounts) ~ hoesch.scale, data = jsub)
  jslope <- summary(jfit)$coefficients["hoesch.scale", "Estimate"]
  jslope.dat <- data.frame(slope = jslope, mark = jmark, stringsAsFactors = FALSE)
  return(jslope.dat)
}) 

jslopes.peak.discrete <- lapply(jmarks, function(jmark){
  jsub <- dat.glmpca.peak.lst[[jmark]]
  jfit <- lm(formula = log2(peakcounts / spikeincounts) ~ cellcycle.str, data = jsub)
  # jslope <- summary(jfit)$coefficients[c("cellcycle.str1_S", "cellcycle.str2_G2/M"), "Estimate"]
  jslope <- summary(jfit)$coefficients[c("cellcycle.str2_G2/M"), "Estimate"]
  jslope.dat <- data.frame(mark = jmark, slope = jslope, mark = jmark, stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(foldchange = 2 ^ jslope)
  return(jslope.dat)
}) %>%
  bind_rows()

jslopes.peak.winsorized <- lapply(jmarks, function(jmark){
  jsub <- dat.glmpca.peak.lst[[jmark]]
  jfit <- lm(formula = ltr.wins ~ hoesch.scale, data = jsub)
  jslope <- summary(jfit)$coefficients["hoesch.scale", "Estimate"]
  jslope.dat <- data.frame(slope = jslope, mark = jmark, stringsAsFactors = FALSE)
  return(jslope.dat)
}) 





