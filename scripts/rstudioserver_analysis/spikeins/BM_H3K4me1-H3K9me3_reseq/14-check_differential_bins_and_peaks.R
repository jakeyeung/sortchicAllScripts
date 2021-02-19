# Jake Yeung
# Date of Creation: 2021-01-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/14-check_differential_bins_and_peaks.R
#


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(forcats)

library(topicmodels)
library(JFuncs)
library(scchicFuncs)

options(scipen=0)

make.plots <- FALSE

pvalcutoff <- 1e-10
# padjcutoff <- 1e-9
padjcutoff <- 1e-50

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned"
fname <- paste0("DE_bins_all_marks_padjcutoff_dists_to_TSS.tweak.", Sys.Date(), ".padj_", padjcutoff, ".pdf")
outpdf <- file.path(outdir, fname)


if (make.plots){
  pdf(outpdf, useDingbats = FALSE)
}


# Constants ---------------------------------------------------------------


jkeeptop <- 150
jlow.in.k9 <- TRUE

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load MATs  --------------------------------------------------------------

# ctypes <- c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs")
ctypes <- c("Eryths", "Bcells", "NKs", "pDCs", "Granulocytes", "DCs", "Basophils", "HSPCs")
ctypes.k9me3 <- c("Eryths", "Bcells", "Granulocytes", "HSPCs")

dat.metas <- lapply(jmarks, function(jmark){
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final/umaps_final_get_marks.", jmark, ".metadata.2020-12-17.txt")
  fread(inf)
})

dat.metas <- lapply(jmarks, function(jmark){
  dat.metas.tmp <- dat.metas[[jmark]]
  if (jmark != "H3K9me3"){
    dat.metas.tmp$jrep2 <- sapply(dat.metas.tmp$jrep, function(x) ifelse(x == "rep1old", "zold", "anew"))
  } else {
    dat.metas.tmp$jrep2 <- sapply(dat.metas.tmp$jrep, function(x) ifelse(x != "rep1old", "zold", "anew"))
  }
  return(dat.metas.tmp)
})

jmetas.pretty.lst <- lapply(jmarks, function(jmark){
  jmeta <- dat.metas[[jmark]]
  if (jmark != "H3K9me3"){
    jmeta$cluster <- factor(jmeta$cluster, levels = ctypes)
  } else { 
    jmeta$cluster <- factor(jmeta$cluster, levels = ctypes.k9me3)
  }
  jmeta <- jmeta %>% arrange(cluster, jrep)
})

cells.keep.lst <- lapply(jmetas.pretty.lst, function(jdat){
  jdat$cell
})


# Get DE outputs ----------------------------------------------------------

jfits.lst.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jinf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.RData"))
  load(jinf, v=T)
  return(jfits.lst)
})


params.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfits.lst <- jfits.lst.lst[[jmark]]
  jnames <- names(jfits.lst); names(jnames) <- jnames
  params.dat.all <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xkeep <- grepl("^Cluster.*.Estimate$", x = names(x))
    jparams <- x[xkeep]
    data.frame(bin = jname, param = names(jparams), estimate = unlist(jparams), mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  if (jmark == "H3K9me3"){
    params.dat.all <- params.dat.all %>% 
      mutate(param = gsub("Eryth", "Eryths", param),
             param = gsub("Lymphoid", "Bcells", param))
  }
  # make params more readable
  params.dat.all$ctype <- params.dat.all$param
  params.dat.all$ctype <- gsub("Cluster", "", params.dat.all$ctype)
  params.dat.all$ctype <- gsub(".Estimate", "", params.dat.all$ctype)
  return(params.dat.all)
})


pvals.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfits.lst <- jfits.lst.lst[[jmark]]
  jnames <- names(jfits.lst); names(jnames) <- jnames
  lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, padj = p.adjust(xvec), stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
})

pvals.lst2 <- lapply(jfits.lst.lst$H3K9me3, function(x) x$pval)
k9.bins <- which(pvals.lst2 < pvalcutoff)

pvals.adj.lst <- lapply(pvals.lst, function(jdat){
  jdat$padj <- p.adjust(jdat$pval, method = "BH")
  return(jdat)
})

pfilt <- lapply(pvals.adj.lst, function(jdat){
  subset(jdat, bin %in% names(k9.bins))
})

# plot(-log10(pfilt$H3K9me3$pval), -log10(pfilt$H3K9me3$padj), log = "xy")



# Print statemetns --------------------------------------------------------


print("Range of padj values for significant K9 bins")
print(lapply(pfilt, function(jdat) range(jdat$padj)))

print("Number of significant k9 bins")
print(lapply(pfilt, function(jdat) dim(jdat)))


# How many total heterochromatin regions do we have?  ---------------------

# from peaks

# for K27me3 
inf.lda.lst <- lapply(jmarks, function(jmark){
  if (jmark != "H3K9me3"){
    # inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/tm_results_from_projs/tm_result_old_to_new.", jmark, ".2020-12-28.RData"))
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins.from_sitecount_mat.from_same_annot_file/lda_outputs.count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.K-30.binarize.FALSE/ldaOut.count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.K-30.Robj"))
  } else {
    # inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/tm_results_from_projs/tm_result_new_to_old.", jmark, ".2020-12-28.RData"))
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins/lda_outputs.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.binarize.FALSE/ldaOut.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.Robj"))
  }
  assertthat::assert_that(file.exists(inf.lda.tmp))
  
  return(inf.lda.tmp)
})

out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.lda <- inf.lda.lst[[jmark]]
  load(inf.lda, v=T)  # out.lda, count.mat
  tm.result <- posterior(out.lda)
  return(list(tm.result = tm.result, count.mat = count.mat))
})

count.mat.lst <- lapply(out.lst, function(jout){
  jout$count.mat
})

# From bins ---------------------------------------------------------------



# load glmpca outputs? 

inf.lda.bins.lst <- lapply(jmarks, function(jmark){
  if (jmark != "H3K27me3"){
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
  } else {
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_", jmark, ".cleaned.varfilt_2.K-30.Robj"))
  }
  return(inf.lda.tmp)
})


count.mat.bins.lst <- lapply(inf.lda.bins.lst, function(inf.bins){
  print(inf.bins)
  inf.bins <- inf.bins
  load(inf.bins, v=T)  # out.lda, count.mat
  return(count.mat)
})


# Define heterochromatin regions  -----------------------------------------

# jthres.lst <- c(3, 3, 3, 4)
jthres.lst <- c(3, 4, 2.3, 4)
names(jthres.lst) <- jmarks

cvec.norm.log.lst <- lapply(jmarks, function(jmark){
  jthres <- jthres.lst[[jmark]]
  print(jmark)
  cvec <- rowSums(count.mat.bins.lst[[jmark]])
  cvec.norm <- 10^6 * cvec / sum(cvec)
  cvec.norm.log <- log2(cvec.norm + 1)
  # cvec <- DescTools::Winsorize(cvec.norm, probs = c(0.01, 0.99))
  plot(density(cvec.norm.log), main = jmark)
  abline(v = jthres)
  names(cvec.norm.log) <- sapply(names(cvec.norm.log), function(x) strsplit(x, ";")[[1]][[1]])
  # remake them centers
  
  coord <- names(cvec.norm.log)
  jstart <- sapply(coord, JFuncs::GetStart, returnAsInt = TRUE)
  jend <- sapply(coord, JFuncs::GetEnd, returnAsInt = TRUE)
  jchromo <- sapply(coord, JFuncs::GetChromo, add.chr = FALSE)
  jmid <- jstart + round((jstart - jend) / 2)
  # new coord
  coordnew <- paste(jchromo, paste(jmid - 1, jmid + 1, sep = "-"), sep = ":")
  names(cvec.norm.log) <- coordnew
  
  return(cvec.norm.log)
})

bins.high.lst <- lapply(jmarks, function(jmark){
  cvec.filt.i <- cvec.norm.log.lst[[jmark]] > jthres.lst[[jmark]]
  bins.high.tmp <- names(cvec.norm.log.lst[[jmark]])[cvec.filt.i]
})

print(lapply(bins.high.lst, length))
print(lapply(bins.high.lst, head))



# Calculate total number of base pairs by peaks ---------------------------

rnames.dat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jcount <- count.mat.lst[[jmark]]
  coord <- sapply(rownames(jcount), function(x) strsplit(x, ";")[[1]][[1]])
  jstart <- sapply(coord, JFuncs::GetStart, returnAsInt = TRUE)
  jend <- sapply(coord, JFuncs::GetEnd, returnAsInt = TRUE)
  jchromo <- sapply(coord, JFuncs::GetChromo, add.chr = TRUE)
  jdat <- data.frame(Coord = paste("chr", coord, sep = ""), Chromo = jchromo, Start = jstart, End = jend, stringsAsFactors = FALSE)  %>%
    ungroup() %>%
    mutate(Dist = End - Start,
           mark = jmark)
  return(jdat)
})

# get distributoin of size of peaks
mdist <- ggplot(rnames.dat.lst %>% bind_rows(), aes(x = Dist, fill = mark)) + 
  geom_density(alpha = 0.25) +
  scale_x_log10() + 
  facet_wrap(~mark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(mdist)

rnames.dat.sum <- rnames.dat.lst %>%
  bind_rows() %>%
  group_by(mark) %>%
  summarise(Dist = sum(Dist)) %>%
  mutate(gbtotal = Dist / 10^9)

mbar <- ggplot(rnames.dat.sum, aes(x = mark, y = Dist)) + 
  geom_col() 

print(mbar)


# Are regions near genes? -------------------------------------------------

inf.tsspretty <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.50000.again.cut.pretty.bed"
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
# bins.annot <- AnnotateCoordsFromList(coords.vec = rnames.dat.lst$H3K9me3$Coord, inf.tss = inf.tsspretty, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)
bins.annot <- AnnotateCoordsFromList(coords.vec = pvals.adj.lst$H3K9me3$bin, inf.tss = inf.tsspretty, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

head(bins.annot$regions.annotated)

plot(density(bins.annot$regions.annotated$distanceToTSS))

bins.annot.filt <- subset(bins.annot$regions.annotated, region_coord %in% names(k9.bins))

plot(density(abs(bins.annot.filt$distanceToTSS + 1)), log = "x")
# plot(density(abs(bins.annot$regions.annotated$distanceToTSS + 1)), log = "x")
plot(density(abs(bins.annot$regions.annotated$distanceToTSS + 1)), log = "x")
plot(density(abs(bins.annot.filt$distanceToTSS + 1)), log = "x")

# Are there regions that aren't inside "peaks" ?  -------------------------


dat.regions <- data.frame(region_coord = pvals.lst$H3K9me3$bin, stringsAsFactors = FALSE) %>%
  ungroup() %>%
  mutate(seqnames = sapply(region_coord, JFuncs::GetChromo),
         start = sapply(region_coord, JFuncs::GetStart, returnAsInt = TRUE),
         end = sapply(region_coord, JFuncs::GetEnd, returnAsInt = TRUE))
regions.gr <- makeGRangesFromDataFrame(dat.regions, keep.extra.columns = TRUE)

dat.peaks <- data.frame(rnames.dat.lst$H3K9me3) %>%
  dplyr::rename(seqnames = Chromo, 
                start = Start,
                end = End,
                region_coord = Coord)

peaks.gr <- makeGRangesFromDataFrame(dat.peaks, keep.extra.columns = TRUE)

# heterochromatin.regions <- findOverlaps(regions.gr, peaks.gr, type = "within")
gr.out <- findOverlaps(regions.gr, peaks.gr, type = "any")

dat.overlap <- data.frame(regions.gr[queryHits(gr.out), ], peaks.gr[subjectHits(gr.out), ], stringsAsFactors = FALSE)

bins.overlap <- unique(dat.overlap$region_coord)

# heterochromatin.regions2 <- findOverlaps(peaks.gr, regions.gr, type = "any")
subset(dat.regions, start > 10117500 & end < 1022500)


# Distance of all regions to TSS  -----------------------------------------

jsub <- bins.annot$regions.annotated %>% filter(region_coord %in% bins.overlap)
jsub.hits <- bins.annot$regions.annotated %>% filter(region_coord %in% names(k9.bins))
plot(density(abs(jsub$distanceToTSS)))
plot(density(abs(jsub.hits$distanceToTSS)))

plot(density(log10(abs(jsub$distanceToTSS) + 1)))
plot(density(log10(abs(jsub.hits$distanceToTSS) + 1)))

plot(density(log10(abs(jsub$distanceToTSS) + 1)))
plot(density(log10(abs(jsub.hits$distanceToTSS) + 1)))


# Check significant bins in K4me3, K4me1, and K27me3  ---------------------

signif.bins.lst <- lapply(jmarks, function(jmark){
  jdat <- pvals.lst[[jmark]]
  if (jmark != "H3K9me3"){
    jdat.filt <- subset(jdat, padj < padjcutoff)
  } else {
    jdat.filt <- subset(jdat, pval < pvalcutoff)
  }
  coord <- jdat.filt$bin
  jstart <- sapply(coord, JFuncs::GetStart, returnAsInt = TRUE)
  jend <- sapply(coord, JFuncs::GetEnd, returnAsInt = TRUE)
  jchromo <- sapply(coord, JFuncs::GetChromo, add.chr = FALSE)
  jmid <- jstart + round((jstart - jend) / 2)
  # new coord
  coordnew <- paste(jchromo, paste(jmid - 1, jmid + 1, sep = "-"), sep = ":")
  
  jdat <- data.frame(Coord = coordnew, 
                     Chromo = jchromo, 
                     Start = jstart, End = jend, 
                     Mid = jmid, 
                     stringsAsFactors = FALSE)  %>%
    ungroup() %>%
    mutate(Dist = End - Start,
           mark = jmark)
  
  
})

bins.annot.filt.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  bins.annot.tmp <- AnnotateCoordsFromList(coords.vec = signif.bins.lst[[jmark]]$Coord, inf.tss = inf.tsspretty, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)
  bins.annot.tmp <- bins.annot.tmp$regions.annotated
  bins.annot.tmp$mark <- jmark
  return(bins.annot.tmp)
})

bins.annot.high.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  bins.annot.tmp <- AnnotateCoordsFromList(coords.vec = bins.high.lst[[jmark]], inf.tss = inf.tsspretty, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)
  bins.annot.tmp <- bins.annot.tmp$regions.annotated
  bins.annot.tmp$mark <- jmark
  return(bins.annot.tmp)
})


# get diostributio nof TSS
print("Median distances to TSS")
lapply(bins.annot.filt.lst, function(jdat) median(jdat$distanceToTSS))

# get nbins as str
nbins.str <- paste(lapply(jmarks, function(jmark){
  nbins <- nrow(bins.annot.filt.lst[[jmark]])
  jstr <- paste(jmark, nbins, sep = ":N=")
}) %>%
  unlist(), collapse = ", ")

nbins.str.high <- paste(lapply(jmarks, function(jmark){
  nbins <- nrow(bins.annot.high.lst[[jmark]])
  jstr <- paste("high", jmark, nbins, sep = ":N=")
}) %>%
  unlist(), collapse = ", ")

dat.bins.annot.filt.long <- bins.annot.filt.lst %>%
  bind_rows()

ggplot(dat.bins.annot.filt.long, aes(x = abs(distanceToTSS) + 1, fill = mark)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + 
  scale_x_log10() + 
  ggtitle(paste("Padjcutoff:", padjcutoff), nbins.str) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(bins.annot.high.lst %>% bind_rows(), aes(x = abs(distanceToTSS) + 1, fill = mark)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + 
  scale_x_log10() + 
  ggtitle(nbins.str.high) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.bins.annot.filt.long, aes(x = abs(distanceToTSS) + 1)) + 
  geom_density(alpha = 0.25, fill = "grey50") + 
  theme_bw() + 
  scale_x_log10() + 
  ggtitle(paste("Padjcutoff:", padjcutoff), nbins.str) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(bins.annot.high.lst %>% bind_rows(), aes(x = abs(distanceToTSS) + 1)) + 
  geom_density(alpha = 0.25, fill = "grey50") + 
  theme_bw() + 
  scale_x_log10() + 
  ggtitle(nbins.str.high) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())





ggplot(dat.bins.annot.filt.long, aes(x = forcats::fct_reorder(.f = mark, .x = distanceToTSS, .fun = median, .desc = TRUE), 
                                     y = abs(distanceToTSS) + 1, 
                                     fill = mark)) + 
  geom_violin(alpha = 0.5) +  
  theme_bw() + 
  xlab("") + 
  scale_y_log10() + 
  ggtitle(paste("Padjcutoff:", padjcutoff), nbins.str) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(bins.annot.high.lst %>% bind_rows(), aes(x = forcats::fct_reorder(.f = mark, .x = distanceToTSS, .fun = median, .desc = TRUE), 
                                     y = abs(distanceToTSS) + 1, 
                                     fill = mark)) + 
  geom_violin(alpha = 0.5) +  
  theme_bw() + 
  xlab("") + 
  scale_y_log10() + 
  ggtitle(nbins.str.high) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



ggplot(dat.bins.annot.filt.long, aes(x = forcats::fct_reorder(.f = mark, .x = distanceToTSS, .fun = median, .desc = TRUE), 
                                     y = abs(distanceToTSS) + 1)) + 
  geom_violin(alpha = 0.5, fill = "grey50") +  
  theme_bw() + 
  xlab("") + 
  scale_y_log10() + 
  ggtitle(paste("Padjcutoff:", padjcutoff), nbins.str) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(bins.annot.high.lst %>% bind_rows(), aes(x = forcats::fct_reorder(.f = mark, .x = distanceToTSS, .fun = median, .desc = TRUE), 
                                     y = abs(distanceToTSS) + 1)) + 
  geom_violin(alpha = 0.5, fill = "grey50") +  
  theme_bw() + 
  xlab("") + 
  scale_y_log10() + 
  ggtitle(nbins.str.high) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())





ggplot(dat.bins.annot.filt.long, aes(x = abs(distanceToTSS) + 1, 
                                     fill = mark)) + 
  geom_histogram(alpha = 0.5) +  
  facet_wrap(~mark, ncol = 1) + 
  theme_bw(10) + 
  geom_vline(xintercept = c(25000), linetype = "dotted") + 
  scale_x_log10() + 
  ggtitle(paste("Padjcutoff:", padjcutoff), nbins.str) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(bins.annot.high.lst %>% bind_rows(), aes(x = abs(distanceToTSS) + 1, 
                                     fill = mark)) + 
  geom_histogram(alpha = 0.5) +  
  facet_wrap(~mark, ncol = 1) + 
  theme_bw(10) + 
  geom_vline(xintercept = c(25000), linetype = "dotted") + 
  scale_x_log10() + 
  ggtitle(nbins.str.high) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



ggplot(dat.bins.annot.filt.long, aes(x = abs(distanceToTSS) + 1)) + 
  geom_histogram(alpha = 0.5, fill = "grey50") +  
  facet_wrap(~mark, ncol = 1) + 
  theme_bw(10) + 
  geom_vline(xintercept = c(25000), linetype = "dotted") + 
  scale_x_log10() + 
  ggtitle(paste("Padjcutoff:", padjcutoff), nbins.str) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(bins.annot.high.lst %>% bind_rows(), aes(x = abs(distanceToTSS) + 1)) + 
  geom_histogram(alpha = 0.5, fill = "grey50") +  
  facet_wrap(~mark, ncol = 1) + 
  theme_bw(10) + 
  geom_vline(xintercept = c(25000), linetype = "dotted") + 
  scale_x_log10() + 
  ggtitle(nbins.str.high) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



ggplot(bins.annot.high.lst %>% bind_rows(), aes(x = abs(distanceToTSS) + 1, 
                                     fill = mark)) + 
  geom_histogram(alpha = 0.5, bins = 100) +  
  facet_wrap(~mark, nrow = 1) + 
  theme_bw(10) + 
  geom_vline(xintercept = c(25000), linetype = "dotted") + 
  scale_x_log10() + 
  ggtitle(nbins.str.high) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(dat.bins.annot.filt.long, aes(x = abs(distanceToTSS) + 1)) + 
  geom_histogram(alpha = 0.5, bins = 100, fill = "grey50") +  
  facet_wrap(~mark, nrow = 1) + 
  theme_bw(10) + 
  geom_vline(xintercept = c(25000), linetype = "dotted") + 
  scale_x_log10() + 
  ggtitle(paste("Padjcutoff:", padjcutoff)) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(bins.annot.high.lst %>% bind_rows(), aes(x = abs(distanceToTSS) + 1)) + 
  geom_histogram(alpha = 0.5, bins = 100, fill = "grey50") +  
  facet_wrap(~mark, nrow = 1) + 
  theme_bw(10) + 
  geom_vline(xintercept = c(25000), linetype = "dotted") + 
  scale_x_log10() + 
  ggtitle(nbins.str.high) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


# Make boxplots to compare?  ----------------------------------------------

bins.merge <- rbind(dat.bins.annot.filt.long %>% mutate(is.signif = TRUE), 
                    bins.annot.high.lst %>% bind_rows() %>% mutate(is.signif = FALSE))

bins.avg <- bins.merge %>%
  group_by(mark, is.signif) %>%
  summarise(distanceToTSS.median = median(abs(distanceToTSS)))

medians.str <- paste(apply(bins.avg, MARGIN = 1, function(x) paste(x, collapse = ",")), collapse = "\n")

ggplot(bins.merge, aes(x = abs(distanceToTSS) + 1, fill = is.signif)) + 
  geom_histogram(alpha = 0.5, bins = 100, position = "stack") + 
  facet_wrap(~mark, nrow = 1) + 
  theme_bw(10) + 
  ggtitle(medians.str) + 
  geom_vline(xintercept = c(25000), linetype = "dotted") + 
  scale_x_log10() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom", plot.title = element_text(size = 5))

ggplot(bins.merge, aes(x = abs(distanceToTSS) + 1, fill = is.signif)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~mark, nrow = 1) + 
  theme_bw(10) + 
  ggtitle(medians.str) + 
  geom_vline(xintercept = c(25000), linetype = "dotted") + 
  scale_x_log10() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom", plot.title = element_text(size = 5))

ggplot(bins.merge, aes(x = abs(distanceToTSS) + 1, fill = is.signif)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~mark, nrow = 1) + 
  theme_bw(10) + 
  ggtitle(medians.str) + 
  geom_vline(xintercept = c(25000), linetype = "dotted") + 
  scale_x_log10() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom", plot.title = element_text(size = 5))


ggplot(bins.merge, aes(x = is.signif, y = abs(distanceToTSS) + 1)) + 
  geom_boxplot() + 
  facet_wrap(~mark, nrow = 1) + 
  ggtitle(medians.str) + 
  theme_bw(10) + 
  geom_hline(yintercept = c(25000), linetype = "dotted") + 
  scale_y_log10() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom", plot.title = element_text(size = 5))

if (make.plots){
  dev.off()
}



# write to output
for (jmark in jmarks){
  print(jmark)
  jtmp <- bins.annot.filt.lst[[jmark]]
  print(jtmp)
  fnametmp <- paste0("DE_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.", jmark, ".", Sys.Date(), ".txt")
  outftmp <- file.path(outdir, fnametmp)
  fwrite(jtmp, file = outftmp, quote = FALSE, sep = "\t")
}

# write high bins
for (jmark in jmarks){
  print(jmark)
  jtmp <- bins.annot.high.lst[[jmark]]
  print(jtmp)
  fnametmp <- paste0("High_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.", jmark, ".", Sys.Date(), ".txt")
  outftmp <- file.path(outdir, fnametmp)
  fwrite(jtmp, file = outftmp, quote = FALSE, sep = "\t")
}
