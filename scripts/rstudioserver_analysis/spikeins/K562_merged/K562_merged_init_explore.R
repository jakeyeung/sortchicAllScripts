# Jake Yeung
# Date of Creation: 2020-09-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_merged/K562_merged_init_explore.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

outpdf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_K562.round2.AllMerged.WithDoubleCells/K562_merged_G1_vs_G2_across_chromosomes.pdf"


pdf(outpdf, useDingbats = FALSE)

# Load data  --------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

jchromos <- c(seq(19), "X", "Y")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


dat.lst <- lapply(jmarks, function(jmark){
  
  inf.bins <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/countTablesAndRZr1only_ByBins.NewFilters/K562_CellCycleSorted_", jmark, ".merged.sorted.tagged.countTable.binsize_50000.csv"))
  inf.chromo <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/countTablesAndRZr1only_ByChromo.NewFilters/K562_CellCycleSorted_", jmark, ".merged.sorted.tagged.countTable.ByChromo.csv"))
  inf.lh <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/countTablesAndRZr1only_TAfrac.NewFilters/K562_CellCycleSorted_", jmark, ".merged.sorted.tagged.LH_counts.demuxbugfixed_mergeplates.csv"))

  dat.bins <- ReadMatSlideWinFormat(inf = inf.bins, as.sparse = TRUE, sort.rnames = TRUE, add.chromo = TRUE)
  dat.chromos <- GetChromoCounts(inf.chromo, spikeinchromo = "J02459.1", chromos.keep = jchromos) %>%
    filter(chromo == "1") %>%
    mutate(mark = jmark)
  dat.lh <- ReadLH.SummarizeTA(inf.lh) %>%
    left_join(., dat.chromos, by = "samp") %>%
    rowwise() %>%
    mutate(cellid = paste("cell", strsplit(samp, "_")[[1]][[2]], sep = ""),
           indx = strsplit(samp, "_")[[1]][[2]],
           rowcoord = GetPlateCoord(cell = cellid, platecols = 24, is.zero.base = FALSE)[[1]],
           colcoord = GetPlateCoord(cell = cellid, platecols = 24, is.zero.base = FALSE)[[2]],
           is.empty = rowcoord <= 8 & colcoord == 1) %>%
    mutate(mark = jmark)
  return(list(dat.bins = dat.bins, dat.chromos = dat.chromos, dat.lh = dat.lh))
})



dat.lh <- lapply(dat.lst, function(jdat) jdat$dat.lh) %>%
  bind_rows()
jtitle <- paste("K562 merged")

ggplot(dat.lh, aes(x = log2(chromocounts / spikeincounts), y = TA.frac, color = is.empty)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + 
  ggtitle(jtitle) + 
  facet_wrap(~mark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(dat.lh, aes(x = log10(chromocounts), y = TA.frac, color = is.empty)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + 
  ggtitle(jtitle) + 
  facet_wrap(~mark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(dat.lh, aes(x = log10(spikeincounts), y = TA.frac, color = is.empty)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + 
  ggtitle(jtitle) + 
  facet_wrap(~mark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(dat.lh, aes(x = log10(spikeincounts), y = TA.frac, color = is.empty)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + 
  ggtitle(jtitle) + 
  facet_wrap(~mark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(dat.lh, aes(x = log10(spikeincounts), y = log10(chromocounts)), color = is.empty) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + 
  geom_abline(slope = 1, intercept = 0) + 
  ggtitle(jtitle) + 
  facet_wrap(~mark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


# Load good cells ---------------------------------------------------------


fracmin <- 0.5
# chromocountmin <- list(H3K4me1 = 2000, H3K4me3 = 2000, H3K27me3 = 2000, H3K9me3 = 2000)
chromocountmin <- list(H3K4me1 = 1000, H3K4me3 = 1000, H3K27me3 = 1000, H3K9me3 = 1000)
log2fcmin <- 2

dat.lh.cutoffs <- lapply(jmarks, function(jmark){
  dat.lst[[jmark]]$dat.lh %>%
    mutate(chromocutoff = chromocountmin[[jmark]])
}) %>%
  bind_rows() %>%
  mutate(is.good = !is.empty & chromocounts >= chromocutoff & log2(chromocounts / spikeincounts) > log2fcmin)

ggplot(dat.lh.cutoffs, aes(x = log10(chromocounts), y = TA.frac, color = is.good)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + 
  ggtitle(jtitle) + 
  facet_wrap(~mark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


# Load hoescht annots -----------------------------------------------------

inf.hoescht <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/spikein_fits_cellcycle.hoesch_staining/hoesch_on_K562_plates.rds"
dat.hoescht.all <- readRDS(inf.hoescht)


unique(dat.hoescht.all$experi)
# jmark <- "H3K4me1"
# jmark.hoescht <- "4me1.*.4.*_index$"



# jmark <- "H3K27me3"

jprefix <- "H3K"
# jsuffix <- "4me3"
# jsuffix <- "9me3"
jsuffix <- "27me3"
jsuffix <- "4me1"
jsuffix <- "9me3"
jsuffix <- "4me3"

jmark <- paste0(jprefix, jsuffix)
# jmark.hoescht <- "27me3.*.4.*_index$"
if (jsuffix == "4me1"){
  jmark.hoescht <- paste0(jsuffix, ".*.4.*.index$")
} else {
  jmark.hoescht <- paste0(jsuffix, ".*.4_index$")
}

dat.lh.sub <- subset(dat.lh.cutoffs, mark == jmark) %>%
  mutate(experi2 = strsplit(samp, "_")[[1]][[1]])

dat.hoescht <- subset(dat.hoescht.all, grepl(jmark.hoescht, x = experi))
assertthat::assert_that(nrow(dat.hoescht) > 0)

experi.str <- unique(dat.lh.sub$experi2)
assertthat::assert_that(length(experi.str) == 1)

# rename experi to match glmpca
dat.hoescht <- dat.hoescht %>%
  rowwise() %>%
  mutate(experi2 = experi.str)

jdat <- left_join(dat.lh.sub, dat.hoescht, by = c("rowcoord" = "row.indx", "colcoord" = "col.indx", "experi2" = "experi2")) %>%
  filter(is.good)

ggplot(jdat %>% filter(is.good), aes(x = hoesch, y = log2(chromocounts / spikeincounts), color = cellcycle)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jmark, jmark.hoescht) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jdat %>% filter(is.good), aes(x = cellcycle, y = log2(chromocounts / spikeincounts), color = cellcycle)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1) + 
  theme_bw() + 
  ggtitle(jmark, jmark.hoescht) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jdat %>% filter(is.good), aes(x = cellcycle, y = log10(chromocounts), color = cellcycle)) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jmark, jmark.hoescht) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jdat %>% filter(is.good), aes(x = cellcycle, y = log10(spikeincounts), color = cellcycle)) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jmark, jmark.hoescht) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# fit 

jdat.forfit <- subset(jdat, cellcycle != "1_S")
jfit <- lm(formula = log2(chromocounts / spikeincounts) ~ cellcycle, data = jdat.forfit)

print(jfit)

ggplot(jdat %>% filter(is.good & cellcycle != "1_S"), aes(x = cellcycle, y = log2(chromocounts / spikeincounts), color = cellcycle)) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jmark, jmark.hoescht) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Check differences by chromosomes  ---------------------------------------

dat.bychromo <- lapply(jmarks, function(jmark){
  inf.chromo <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/countTablesAndRZr1only_ByChromo.NewFilters/K562_CellCycleSorted_", jmark, ".merged.sorted.tagged.countTable.ByChromo.csv"))
  dat.chromos <- GetChromoCounts(inf.chromo, spikeinchromo = "J02459.1", chromos.keep = jchromos)
  return(dat.chromos)
}) %>%
  bind_rows()

dat.bychromo.merge <- left_join(dat.bychromo, 
                                subset(dat.lh.cutoffs, select = c(samp, mark, cellid, indx, rowcoord, colcoord, is.empty, is.good)), by = "samp") %>%
  filter(is.good) %>%
  rowwise() %>%
  mutate(cellcycle = AddCellCycleLabel(colcoord))
    
ggplot(dat.bychromo.merge, aes(x = chromo, y = log2(counts / spikeincounts), fill = cellcycle)) + 
  geom_boxplot() + 
  facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.bychromo.merge %>% filter(chromo == "9"), aes(x = mark, y = log2(counts / spikeincounts), fill = cellcycle)) + 
  geom_boxplot() + 
  facet_wrap(~chromo) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")



# Fit linear model (robust) for every chromosome for every mark  ----------



jchromo <- "1"
jmark <- "H3K4me1"

chromos.keep <- c(seq(19), "X")
# jchromos <- unique((dat.bychromo.merge %>% filter(chromo %in% chromos.keep))$chromo)

names(chromos.keep) <- chromos.keep

jfits.all <- lapply(chromos.keep, function(jchromo){
  print(jchromo)
  jout <- lapply(jmarks, function(jmark){
    print(jmark)
    jsub.tmp <- subset(dat.bychromo.merge, chromo == jchromo & mark == jmark)
    # winsorize?
    jsub.tmp <- jsub.tmp %>%
      group_by(cellcycle) %>%
      mutate(counts.winsorized = DescTools::Winsorize(x = counts, probs = c(0.05, 0.95)))
    
    jvar.tmp <- jsub.tmp %>%
      group_by(cellcycle) %>%
      summarise(jvar = var(counts.winsorized))
    
    # print(jvar.tmp)
    
    if (any(jvar.tmp$jvar == 0)){
      print("Skipping...")
      m <- ggplot(jsub.tmp, aes(x = cellcycle, y = log2(counts.winsorized / spikeincounts), fill = cellcycle)) + 
        geom_boxplot() + 
        geom_point()  + 
        theme_bw() + 
        ggtitle(paste("Problematic:", jmark, jchromo)) + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      print(m)
      return(NULL)
    }
    
    
    jfit.tmp <- lm(formula = log2(counts.winsorized / spikeincounts) ~ 1 + cellcycle, data = jsub.tmp)
    jfit.null.tmp <- lm(formula = log2(counts.winsorized / spikeincounts) ~ 1, data = jsub.tmp)
    jfit.anova <- anova(jfit.null.tmp, jfit.tmp)
    
    # summarize fits
    jfit.sum <- data.frame(param = rownames(summary(jfit.tmp)$coefficients), as.data.frame(summary(jfit.tmp)$coefficients), stringsAsFactors = FALSE)
    jfit.sum$mark <- jmark
    jfit.sum$chromo <- jchromo
    jfit.sum$pval <- jfit.anova$`Pr(>F)`[[2]]
    return(jfit.sum)
  })
  return(jout)
})

jfits.all.merged <- lapply(jfits.all, function(jdat){
  return(jdat %>% bind_rows())
}) %>%
  bind_rows()

jfits.all.merged$mark <- factor(jfits.all.merged$mark, jmarks)

# plot estimates
ggplot(jfits.all.merged %>% filter(param == "cellcycle2_G2/M"), aes(x = chromo, y = Estimate)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  facet_wrap(~mark)

ggplot(jfits.all.merged %>% filter(param == "cellcycle2_G2/M"), aes(x = mark, y = Estimate)) + 
  geom_boxplot() + 
  geom_point(mapping = aes(color = chromo)) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(jfits.all.merged %>% filter(param != "(Intercept)"), 
       aes(x = chromo, y = Estimate, color = param, ymin = Estimate - 1.96 * Std..Error, ymax = Estimate + 1.96 * Std..Error)) +  
  geom_point() + 
  geom_errorbar() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
  facet_wrap(~mark, nrow = 1)

ggplot(jfits.all.merged %>% filter(param == "cellcycle2_G2/M"), 
       aes(x = chromo, y = Estimate, ymin = Estimate - 1.96 * Std..Error, ymax = Estimate + 1.96 * Std..Error)) +  
  geom_point() + 
  geom_errorbar() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
  facet_wrap(~mark, nrow = 1)

dev.off()
# # get cellcycle
# dat.cellcycle <- dat.lh.cutoffs %>%
#   rowwise() %>%
#   mutate(cellcycle = AddCellCycleLabel(colcoord)) %>%
#   dplyr::select(-total.count) %>%
#   dplyr::rename(cell = samp)

# # for buys
# for (jmark in jmarks){
#   jtmp <- subset(dat.cellcycle, mark == jmark)
#   fwrite(dat.cellcycle, file = paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/tables_K562_for_Buys/cell_to_cellcycle.", jmark, ".txt"))
# }
# 
# # 


  