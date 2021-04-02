# Jake Yeung
# Date of Creation: 2020-11-07
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/1-QC_plots_K562.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(ggrastr)

library(topicmodels)

library(JFuncs)
library(scchicFuncs)

# Load K562 for QC plots --------------------------------------------------

jchromos <- paste("chr", c(seq(22), "X", "Y"), sep = "")
jchromos.nochr <- paste(c(seq(22), "X", "Y"), sep = "")
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

outdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/K562_again")
dir.create(outdir)
outpdf <- file.path(outdir, paste0("K562_QC_plots.", Sys.Date(), ".add_fracnonzeros.pdf"))
outtxt <- file.path(outdir, paste0("K562_QC_plots.", Sys.Date(), ".add_fracnonzeros.txt"))

pdf(outpdf, useDingbats = FALSE)

jsuffix <- "G1filt"
hubprefix <- "/home/jyeung/hub_oudenaarden"

indir.chromo <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/tagged_bams/G1_sorted/merged_bams/countTablesAndRZr1only_ByChromo.NewFilters")
indir.lh <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/tagged_bams/G1_sorted/merged_bams/countTablesAndRZr1only_TAfrac.NewFilters")
indir.counts <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/tagged_bams/G1_sorted/merged_bams/countTablesAndRZr1only_CountTables.NewFilters")









jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")

fracmin <- 0.5
chromocountmin <-list(H3K4me1 = 500, H3K4me3 = 500, H3K27me3 = 1000, H3K9me3 = 1000)
log2fcmin <- 2.5
dat.cutoffs <- data.frame(mark = names(chromocountmin), chromocountsmin = unlist(chromocountmin), stringsAsFactors = FALSE)

infs.chromo <- list.files(indir.chromo, pattern = "K562-EtOH-.*.csv", full.names = TRUE)
assertthat::assert_that(length(infs.chromo) > 0)


dat.chromos <- lapply(infs.chromo, function(inf){
  dat.filt.long <- GetChromoCounts(inf, spikeinchromo = jspikeinchromo, chromos.keep = NA)
}) %>%
  bind_rows()



# Load mats ---------------------------------------------------------------

mats <- lapply(jmarks, function(jmark){
  fname <- paste0("K562-EtOH-", jmark, ".G1sorted.merged.sorted.tagged.countTable.csv")
  inf.tmp <- file.path(indir.counts, fname)
  ReadMatSlideWinFormat(inf.tmp)
})

# get frac nonzeros
fracnonzeros.global <- lapply(mats, function(jmat){
  apply(jmat, 2, function(jcol) nnzero(jcol) / length(jcol))
})

# make into dat
dat.fracnonzeros <- lapply(jmarks, function(jmark){
  data.frame(cell = names(fracnonzeros.global[[jmark]]), mark = jmark, fracnonzeros = fracnonzeros.global[[jmark]], stringsAsFactors = FALSE)
}) %>%
  bind_rows()

# calculate thresholds
jthres <- 2
bad.cells.lst <- lapply(jmarks, function(jmark){
  nonzeros.frac.tmp <- fracnonzeros.global[[jmark]]
  jsd.global <- mad(nonzeros.frac.tmp)
  jmean.global <- median(nonzeros.frac.tmp)
  jthres.frac <- jmean.global + jsd.global * jthres
  plot(density(nonzeros.frac.tmp), main = jmark, xlab = "Fraction of nonzero cuts globally")
  abline(v = jmean.global)
  abline(v = jthres.frac, lty = 2)
  bad.cells <- names(nonzeros.frac.tmp)[which(nonzeros.frac.tmp > jthres.frac)]
  return(list(bad.cells = bad.cells, jthres.frac = jthres.frac))
})
bad.cells.merged <- lapply(bad.cells.lst, function(x) x$bad.cells) %>%
  unlist()
jthres.frac.lst <- lapply(bad.cells.lst, function(x) x$jthres.frac)


# Load LH counts ----------------------------------------------------------

# indir.lh <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN6969/K562/tagged_bams/RZcounts.NewFilters")
infs.lh <- list.files(indir.lh, pattern = "K562-EtOH-.*.csv", full.names = TRUE)

dat.lh <- lapply(infs.lh, function(inf){
  dat.filt.long.lh <- ReadLH.SummarizeTA(inf)
}) %>%
  bind_rows() %>%
  mutate(experi = ClipLast(samp, jsep = "_"))


chromocounts <- subset(dat.chromos, chromo == "1", select = c(samp, chromocounts, spikeincounts))

dat.lh <- left_join(dat.lh, chromocounts)
print(dim(dat.lh))

dat.lh$mark <- sapply(dat.lh$experi, function(x) strsplit(x, "-")[[1]][[3]])

dat.lh$mark <- factor(dat.lh$mark, levels = c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"))


dat.lh <- dat.lh %>%
  rowwise() %>%
  mutate(cellid = paste("cell", strsplit(samp, "_")[[1]][[2]], sep = ""),
         indx = strsplit(samp, "_")[[1]][[2]], 
         rowcoord = GetPlateCoord(cell = cellid, platecols = 24, is.zero.base = FALSE)[[1]],
         colcoord = GetPlateCoord(cell = cellid, platecols = 24, is.zero.base = FALSE)[[2]],
         is.empty = rowcoord <= 8 & colcoord == 1)

dat.lh <- dat.lh %>%
  mutate(is.good = chromocounts > chromocountmin[[mark]] & TA.frac > fracmin & log2(chromocounts / spikeincounts) > log2fcmin & !is.empty)


# good.cells <- subset(dat.lh, is.good & !is.empty)$samp

ggplot(dat.lh, aes(x = log10(chromocounts), y = TA.frac)) + 
  geom_point_rast(alpha = 0.25) + 
  facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.lh, aes(x = log2(chromocounts / spikeincounts), y = TA.frac)) + 
  geom_point_rast(alpha = 0.25) + 
  facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.lh, aes(x = log2(chromocounts / spikeincounts), y = TA.frac, color = is.good)) + 
  geom_point_rast(alpha = 0.25) + 
  facet_wrap(~mark, nrow = 1) + 
  geom_vline(xintercept = log2fcmin, linetype = "dotted") + 
  geom_hline(yintercept = fracmin, linetype = "dotted") + 
  xlab("log2(cuts in genome/spikein cuts)") + ylab("Fraction of cuts starting with TA") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

dat.lh <- dat.lh %>%
  mutate(is.good2 = is.good & !samp %in% bad.cells.merged)
  # mutate(is.good3 = !cell %in% bad.cells.merged)


dat.lh <- left_join(dat.lh, subset(dat.fracnonzeros, select = -mark), by = c("samp" = "cell"))
print(dim(dat.lh))



for (jmark in jmarks){
  m.points <- ggplot(dat.lh %>% filter(mark == jmark), aes(x = log2(chromocounts / spikeincounts), y = TA.frac, color = is.good2)) + 
    geom_point_rast(alpha = 0.25) + 
    geom_vline(xintercept = log2fcmin, linetype = "dotted") + 
    geom_hline(yintercept = fracmin, linetype = "dotted") + 
    xlab("log2(cuts in genome/spikein cuts)") + ylab("Fraction of cuts starting with TA") + 
    ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.points)
  
  m.points.test <- ggplot(dat.lh %>% filter(mark == jmark), aes(x = log2(chromocounts / spikeincounts), y = fracnonzeros, color = is.good2)) + 
    geom_point(alpha = 0.25) +
    geom_vline(xintercept = log2fcmin, linetype = "dotted") + 
    geom_hline(yintercept = jthres.frac.lst[[jmark]], linetype = "dotted") + 
    xlab("log2(cuts in genome/spikein cuts)") + ylab("Fraction of bins with nonzero number of cuts") + 
    ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.points.test)
  
  m.points.test <- ggplot(dat.lh %>% filter(mark == jmark), aes(x = log2(chromocounts / spikeincounts), y = TA.frac, color = fracnonzeros)) + 
    geom_point(alpha = 0.25) +
    geom_vline(xintercept = log2fcmin, linetype = "dotted") + 
    geom_hline(yintercept = fracmin, linetype = "dotted") + 
    xlab("log2(cuts in genome/spikein cuts)") + ylab("Fraction of cuts starting with TA") + 
    scale_color_viridis_c() + 
    ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.points.test)
  
  m.density <- ggplot(dat.lh %>% filter(mark == jmark), aes(x = log2(chromocounts / spikeincounts), fill = is.good2)) + 
    geom_density(alpha = 0.5) + 
    facet_wrap(~mark, nrow = 1) + 
    xlab("log2(cuts in genome/spikein cuts)") + ylab("Fraction of cuts starting with TA") + 
    ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.density)
  
  m.points2 <- ggplot(dat.lh %>% filter(mark == jmark), aes(x = log10(chromocounts), y = TA.frac, color = is.good2)) + 
    geom_point_rast(alpha = 0.25) + 
    geom_vline(xintercept = log2fcmin, linetype = "dotted") + 
    geom_hline(yintercept = fracmin, linetype = "dotted") + 
    xlab("log10(cuts in genome)") + ylab("Fraction of cuts starting with TA") + 
    ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.points2)
  
  m.density2 <- ggplot(dat.lh %>% filter(mark == jmark), aes(x = log10(chromocounts), fill = is.good2)) + 
    geom_density(alpha = 0.5) + 
    facet_wrap(~mark, nrow = 1) + 
    xlab("log10(cuts in genome)") + ylab("Fraction of cuts starting with TA") + 
    ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.density)
  
  
}
# 
# 
# ggplot(dat.lh, aes(x = log2(chromocounts / spikeincounts), y = TA.frac, color = is.good)) + 
#   geom_point_rast(alpha = 0.25) + 
#   facet_wrap(~mark, nrow = 1) + 
#   xlab("log2(cuts in genome/spikein cuts)") + ylab("Fraction of cuts starting with TA") + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# 
# ggplot(dat.lh, aes(x = log2(chromocounts / spikeincounts), y = log2(chromocounts), color = is.good)) + 
#   geom_point_rast(alpha = 0.25) + 
#   facet_wrap(~mark, nrow = 1) + 
#   geom_vline(xintercept = log2fcmin) + 
#   geom_hline(yintercept = fracmin) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(dat.lh, aes(x = log10(chromocounts), y = TA.frac, color = mark)) + 
#   geom_point_rast(alpha = 0.25) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~experi)
# 
# ggplot(dat.lh, aes(x = log10(chromocounts), y = TA.frac, color = mark)) + 
#   geom_point_rast(alpha = 0.25) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~mark)



# 
# # load K562 count mat with spikeins
# 
# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# 
# # "/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/countTablesAndRZr1only_ByBed.NewFilters.peaks.unfiltbed/K562_AllMerged_H3K27me3.merged.sorted.tagged.bam.countTable.bedfile.BL_with_chr.csv"
# 
# hubprefix <- "/home/jyeung/hub_oudenaarden"
# 
# infs.peaks <- lapply(jmarks, function(jmark){
#   indir.peaks <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/countTablesAndRZr1only_ByBed.NewFilters.peaks.unfiltbed")
#   # fname.peaks <- paste0("K562_AllMerged_", jmark, ".merged.sorted.tagged.bam.countTable.bedfile.BL_with_chr.csv")
#   fname.peaks <- paste0("K562_AllMerged_", jmark, ".merged.sorted.tagged.bam.countTable.bedfile.BL_with_chr.csv")
#   inf.peaks <- file.path(indir.peaks, fname.peaks)
#   assertthat::assert_that(file.exists(inf.peaks))
#   return(inf.peaks)
# })
# 
# mats.peaks <- lapply(infs.peaks, function(inf){
#   ReadMatTSSFormat(inf, as.sparse = TRUE, add.coord = TRUE, sort.rnames = TRUE)
# })
# 
# dat.peaks <- lapply(jmarks, function(jmark){
#   jmat <- mats.peaks[[jmark]]
#   data.frame(peakcounts = colSums(jmat), cell = colnames(jmat), stringsAsFactors = FALSE)
# }) %>%
#   bind_rows()
# 
# 
# ggplot(dat.lh.merge, aes(y = peakcounts / chromocounts, x = log2(chromocounts), color = is.good)) + 
#   geom_point_rast(alpha = 0.25)  + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   coord_cartesian(ylim = c(0, 1)) + 
#   facet_wrap(~mark) 
# 
# ggplot(dat.lh.merge %>% filter(mark == "H3K27me3"), aes(y = log2(peakcounts / chromocounts), x = log2(chromocounts), color = is.good)) + 
#   geom_point_rast(alpha = 0.25)  + 
#   theme_bw() + 
#   facet_wrap(~experi) + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
# 
# 
# 


# Recalculate fraction of cuts in peaks  ----------------------------------


# dat.lh.merge <- left_join(dat.lh, dat.peaks, by = c("samp" = "cell")) %>%
#   left_join(., chromocounts %>% dplyr::select(samp, chromocounts))

# indir.peaks <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/counts_in_peaks_vs_nonpeaks_vs_blacklist.faster.clean")
indir.peaks <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/counts_in_peaks_vs_nonpeaks_vs_blacklist")  # by chromo
assertthat::assert_that(dir.exists(indir.peaks))


dat.peaks <- lapply(jmarks, function(jmark){
  infs <- list.files(path = indir.peaks, pattern = paste0("K562_AllMerged_", jmark, ".merged.sorted.tagged.bam.cuts_in_peaks_vs_nonpeaks.*.csv"), full.names = TRUE)
  dats.csv <- lapply(infs, function(inf){
    dat.csv <- fread(inf, sep = "\t", col.names = c("cell", "cuts_in_peaks", "cuts_notin_peaks", "cuts_in_blacklist", "chromo", "badreads"))
    dat.csv$chromo <- as.character(dat.csv$chromo)
    return(dat.csv)
  }) %>%
    bind_rows()%>%
    filter(chromo %in% jchromos.nochr) %>%
    group_by(cell) %>%
    summarise(cuts_in_peaks = sum(cuts_in_peaks), cuts_notin_peaks = sum(cuts_notin_peaks), cuts_in_blacklist = sum(cuts_in_blacklist))
  dats.csv$mark <- jmark
  return(dats.csv)
}) %>%
  bind_rows()

# add spikeins
dat.peaks.merge <- left_join(dat.peaks, subset(dat.lh, select = c(samp, spikeincounts, is.good, is.good2)), by = c("cell" = "samp"))
 
ggplot(dat.peaks.merge, aes(x = log10(cuts_in_peaks + cuts_notin_peaks), y = cuts_in_peaks / (cuts_in_peaks + cuts_notin_peaks)))  +
  geom_point_rast(alpha = 0.25) + 
  theme_bw() + facet_wrap(~mark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.peaks.merge, aes(x = log2(cuts_in_peaks / spikeincounts), y = cuts_in_peaks / (cuts_in_peaks + cuts_notin_peaks)))  +
  geom_point_rast(alpha = 0.25) + 
  facet_wrap(~mark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.peaks.merge, aes(x = log2(cuts_in_peaks / spikeincounts), y = log2(cuts_in_peaks / cuts_notin_peaks)))  +
  geom_point_rast(alpha = 0.25) + 
  facet_wrap(~mark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.peaks.merge, aes(x = log10(cuts_in_peaks), y = log2(cuts_in_peaks / (cuts_in_peaks + cuts_notin_peaks))))  +
  geom_point_rast(alpha = 0.25) + 
  facet_wrap(~mark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.peaks.merge %>% filter(!is.na(is.good)), aes(x = log10(cuts_in_peaks), y = cuts_in_peaks / (cuts_in_peaks + cuts_notin_peaks), color = is.good2))  +
  geom_point_rast(alpha = 0.25) + 
  facet_wrap(~mark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.peaks.merge %>% filter(!is.na(is.good)), aes(x = cuts_in_peaks / (cuts_in_peaks + cuts_notin_peaks), fill = is.good2))  +
  geom_density(alpha = 0.5) + 
  facet_wrap(~mark, nrow = 1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

# ggplot(dat.peaks.merge, aes(x = cuts_in_peaks / (cuts_in_peaks + cuts_notin_peaks)))  +
#   geom_histogram(alpha = 0.5, fill = "red") + 
#   facet_wrap(~mark, nrow = 1) + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Load LDA  ---------------------------------------------------------------

jmark <- "H3K9me3"
inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_K562_spikein_round2.merged/lda_outputs.K562_count_tables_50000.", jmark, ".AllMerged.K-30.binarize.FALSE/ldaOut.K562_count_tables_50000.", jmark, ".AllMerged.K-30.Robj"))
load(inf.lda, v=T)

tm.result <- posterior(out.lda)
dat.impute.log <- t(log2(tm.result$topics %*% tm.result$terms))
dat.var <- CalculateVarAll(dat.impute.log, jchromos = jchromos)


dat.peaks.merge.var <- left_join(dat.var, dat.peaks.merge)

ggplot(dat.peaks.merge.var, aes(y = cell.var.within.sum.norm, x = (cuts_in_peaks / (cuts_in_peaks + cuts_notin_peaks)))) + 
  geom_point_rast()  + 
  ggtitle(jmark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(dat.peaks.merge.var, aes(y = cell.var.within.sum.norm, x = log2(cuts_in_peaks / cuts_notin_peaks))) + 
  geom_point_rast()  + 
  ggtitle(jmark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

dev.off()


# Write to output ---------------------------------------------------------

fwrite(dat.lh, file = outtxt)

# keep good cells
for (jmark in jmarks){
  print(jmark)
  outtxt.tmp <- file.path(outdir, paste0("K562_QC_plots.", Sys.Date(), ".add_fracnonzeros.", jmark, ".good_cells.txt"))
  dat.lh.tmp <- dat.lh %>% filter(mark == jmark & is.good2)
  print(dim(dat.lh.tmp))
  fwrite(dat.lh.tmp, file = outtxt.tmp)
}



# cutoffs used for analysis



# cells used for analysis








# Plot fraction of reads in peaks -----------------------------------------








# Plot ncells per well analysis -------------------------------------------





