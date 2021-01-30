# Jake Yeung
# Date of Creation: 2021-01-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/3-get_signif_bins_check_features.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

hubprefix <- "/home/jyeung/hub_oudenaarden"




# Load significant bins  --------------------------------------------------

keeptop <- 5000

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

fits.out <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.fits <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData")
  # inf.fits <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.allmarks_spikeins/poisson_fit_bins.", jmark, ".2020-12-25.newannot2.witherrors.MoreBins.newestmeta.RData")
  load(inf.fits, v=T)
  params.long <- SummarizeParamsPvalues(jfits.lst, jmark = jmark, paramname = "Cluster") %>%
    mutate(log2fc = estimate / log(2))
  params.long$padj <- p.adjust(params.long$pval.param)
  means.long <- SummarizeMeanValue(jfits.lst, jmark = jmark)
  jnames <- names(jfits.lst); names(jnames) <- jnames
  pvals.long <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  return(list(params.long = params.long, means.long = means.long, pvals.long = pvals.long))
})


params.long.lst <- lapply(fits.out, function(jfit) jfit$params.long)
means.long.lst <- lapply(fits.out, function(jfit) jfit$means.long)
pvals.long.lst <- lapply(fits.out, function(jfit) jfit$pvals.long)

# Load .metasmetas  -------------------------------------------------------------

indir.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned"
dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  fname.tmp <- paste0("cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
  inf <- file.path(indir.meta, fname.tmp)
  jdat <- fread(inf)
  if (jmark == "H3K9me3"){
    jdat$mark <- jmark
  }
  return(jdat)
})



# Get significant bins ----------------------------------------------------

pval.cutoff.k9me3 <- 10^-10
pval.cutoff.others <- 10^-50

bins.signif.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  if (jmark == "H3K9me3"){
    jbins <- subset(pvals.long.lst[[jmark]], pval < pval.cutoff.k9me3)$bin
  } else {
    jbins <- subset(pvals.long.lst[[jmark]], pval < pval.cutoff.others)$bin
  }
  return(jbins)
})

# try taking top n
bins.signif.lst <- lapply(jmarks, function(jmark){
  (pvals.long.lst[[jmark]] %>% arrange(pval))$bin[1:keeptop]
})



# Count features within each bins -----------------------------------------

# GCs
load("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_gc_analysis/gcs_genomewide.RData", v=T)

gr.gc.dat.dedup <- gr.gc.dat[!duplicated(gr.gc.dat), ]

# plot chr1

library(gtools)
chr1.bins <- subset(gr.gc.dat, grepl("^chr", bname))
chr1.bins <- chr1.bins[gtools::mixedorder(chr1.bins$bname), ]


# Start plots -------------------------------------------------------------




outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K27me3_H3K9me3_analysis"
outpdf <- file.path(outdir, paste0("features_by_bins.", Sys.Date(), ".again.keeptop_", keeptop, ".pdf"))
pdf(file = outpdf, useDingbats = FALSE)

ggplot(chr1.bins, aes(x = bname, y = gc)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# load GC outputs



gc.dat.filt.long <- lapply(jmarks, function(jmark){
  subset(gr.gc.dat.dedup, bname %in% bins.signif.lst[[jmark]]) %>%
    mutate(mark = jmark)
}) %>%
  bind_rows()

# add random
bins.random <- sample(gr.gc.dat.dedup$bname, size = keeptop)
gc.dat.random <- subset(gr.gc.dat.dedup, bname %in% bins.random) %>%
  mutate(mark = "random")

gc.dat.filt.long.withrand <- rbind(gc.dat.filt.long, gc.dat.random)

# ggplot(gc.dat.filt.long, aes(x = mark, y = gc)) + 
ggplot(gc.dat.filt.long.withrand, aes(x = mark, y = gc)) + 
  geom_point() + 
  geom_boxplot() + 
  theme_bw(24) + 
  ylab("GC fraction") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(gc.dat.filt.long.withrand %>% filter(mark != "random"), aes(x = mark, y = gc)) + 
  geom_point() + 
  geom_boxplot() + 
  theme_bw(24) + 
  ylab("GC fraction") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(gc.dat.filt.long.withrand, aes(x = mark, y = log2(gc / (1 - gc)))) + 
  geom_point() + 
  geom_boxplot() + 
  theme_bw(24) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(gc.dat.filt.long.withrand %>% filter(mark != "random"), aes(x = mark, y = log2(gc / (1 - gc)))) + 
  geom_point() + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(gc.dat.filt.long.withrand, aes(x = gc, fill = mark)) + 
  geom_density(alpha = 0.25) + 
  facet_wrap(~mark, ncol = 1) + 
  theme_bw(24) + 
  ylab("GC fraction") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


bins.lst.withrand <- bins.signif.lst
bins.lst.withrand$random <- bins.random




# Check statistical significance?  ----------------------------------------

jfits.gc <- lapply(jmarks, function(jmark){
  jsub <- subset(gc.dat.filt.long.withrand, mark %in% c(jmark, "random")) %>%
    mutate(mark = ifelse(mark == "random", "arandom", mark))
  jsub.fit <- lm(formula = gc ~ mark, data = jsub)
})





# overlaps as function of pcutoff -----------------------------------------

# jmark.compare <- "H3K27me3"
jmark.compare <- "H3K27me3"

jmark.ref <- "H3K9me3"

# jmark.ref <- "H3K27me3"
jmark.ref <- "H3K4me1"
jmark.ref <- "H3K4me3"

pcutoffs <- 10^(-1 * seq(50))


frac.overlaps.all <- lapply(jmarks, function(jmark.ref){
  signif.bins <- bins.signif.lst[[jmark.ref]]
  frac.overlaps <- lapply(jmarks[!jmarks %in% jmark.ref], function(jmark.compare){
    print(paste(jmark.ref, jmark.compare))
    n.overlaps <- lapply(pcutoffs, function(p){
      jbins <- subset(pvals.long.lst[[jmark.compare]], pval < p)$bin
      N <- length(intersect(jbins, signif.bins))
      data.frame(frac.overlaps = N / length(signif.bins), pcutoff = p, mark = jmark.compare, mark.ref = jmark.ref, stringsAsFactors = FALSE)
    }) %>%
      bind_rows()
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

ggplot(frac.overlaps.all, aes(x = -log10(pcutoff), y = frac.overlaps, color = mark)) + 
  geom_point() + 
  facet_wrap(~mark.ref) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# show repressive only
ggplot(frac.overlaps.all %>% filter(mark.ref %in% c("H3K27me3", "H3K9me3")), aes(x = -log10(pcutoff), y = frac.overlaps, color = mark)) + 
  geom_point() + 
  facet_wrap(~mark.ref) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# show active only 
ggplot(frac.overlaps.all %>% filter(mark.ref %in% c("H3K4me1", "H3K4me3")), aes(x = -log10(pcutoff), y = frac.overlaps, color = mark)) + 
  geom_point() + 
  facet_wrap(~mark.ref) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Check pseudbulks of dynamic bins?  -------------------------------------------------------


# load pseudobulks



# load LDA 
outs.all.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  if (jmark != "H3K27me3"){
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
  } else {
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.Robj"))
  }
  load(inf.lda, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})

# jmark.check <- "H3K9me3"
dat.pbulk.lst <- lapply(jmarks, function(jmark.check){
  print(jmark.check)
  cnames.keep.lst <- split(x = dat.metas[[jmark.check]]$cell, f = dat.metas[[jmark.check]]$cluster)
  pbulks.lst <- SumAcrossClusters(outs.all.lst[[jmark.check]]$count.mat, cnames.keep.lst = cnames.keep.lst)
  
  
  mat.pbulk <- bind_rows(pbulks.lst) %>%
    as.data.frame()
  rownames(mat.pbulk) <- rownames(outs.all.lst[[jmark.check]]$count.mat)
  mat.pbulk <- sweep(mat.pbulk, MARGIN = 2, STATS = colSums(mat.pbulk), FUN = "/")
  
  dat.pbulk <- as.matrix(mat.pbulk) %>%
    melt()
  colnames(dat.pbulk) <- c("bin", "ctype", "cuts")
  dat.pbulk$mark <- jmark.check
  return(dat.pbulk)
})


ggplot(dat.pbulk.lst %>% bind_rows(), aes(x = log2(cuts), fill = ctype)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + 
  # facet_wrap(~ctype) + 
  facet_grid(mark~ctype) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pbulk.lst$H3K27me3, aes(x = log2(cuts))) + 
  geom_density(alpha = 0.25, fill = "lightblue") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Check K27me3 dynamic bins and where the K9me3 signal is -----------------

jmark.compare <- "H3K9me3"
jmark.ref <- "H3K27me3"
# jmark.ref <- "H3K4me3"

dat.pbulk.merge.long.lst <- lapply(jmarks, function(jmark.compare){
  print(jmark.compare)
  dat.pbulk.init <- dat.pbulk.lst[[jmark.compare]] %>% mutate(type = paste("AllBins"))
  dat.pbulk.merge.long <- lapply(jmarks, function(jmark.ref){
    dat.pbulk.merge.to.add <- subset(dat.pbulk.lst[[jmark.compare]], bin %in% bins.signif.lst[[jmark.ref]]) %>% mutate(type = paste0(jmark.ref, " DynamicBins"))
  }) %>%
    bind_rows()
  dat.pbulk.merge.long <- rbind(dat.pbulk.init, dat.pbulk.merge.long) %>%
    mutate(mark = jmark.compare)
}) 


m <- ggplot(dat.pbulk.merge.long.lst %>% bind_rows(), aes(x = cuts, fill = type)) + 
  geom_density(alpha = 0.25) + 
  scale_x_log10() + 
  theme_bw() + 
  ggtitle("AllMarks") + 
  facet_grid(type ~ mark) + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)

m <- ggplot(dat.pbulk.merge.long.lst %>% bind_rows(), aes(y = log2(cuts), x = type)) + 
  geom_boxplot() + 
  # geom_violin() + 
  theme_bw() + 
  facet_wrap(~mark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(m)

jsub.tmp <- dat.pbulk.merge.long.lst %>% bind_rows() %>% group_by(bin, mark, type) %>% summarise(cuts = max(cuts))
jsub.tmp.mean <- dat.pbulk.merge.long.lst %>% bind_rows() %>% group_by(bin, mark, type) %>% summarise(cuts = mean(cuts))

m <- ggplot(jsub.tmp, aes(y = log2(cuts), x = type)) + 
  geom_boxplot() + 
  # geom_violin() + 
  theme_bw() + 
  facet_wrap(~mark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(m)


m <- ggplot(jsub.tmp %>% filter(type != "AllBins"), aes(y = log2(cuts), x = mark)) + 
  geom_boxplot() + 
  # geom_violin() + 
  ggtitle("Max across ctypes") + 
  theme_bw() + 
  facet_wrap(~type) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(m)

m <- ggplot(jsub.tmp.mean %>% filter(type != "AllBins"), aes(y = log2(cuts), x = mark)) + 
  geom_boxplot() + 
  # geom_violin() + 
  theme_bw() +  
  ggtitle("Mean across ctypes") + 
  facet_wrap(~type) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(m)

jsub <- dat.pbulk.merge.long.lst %>% bind_rows() %>% group_by(bin, mark, type) %>% summarise(log2.sd.cuts = sd(log2(cuts + 0.001)))

m <- ggplot(jsub, aes(y = log2.sd.cuts, x = type)) + 
  geom_boxplot() + 
  # geom_violin() + 
  theme_bw() + 
  facet_wrap(~mark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(m)





for (jmark in jmarks){
  m <- ggplot(dat.pbulk.merge.long.lst[[jmark]], aes(x = cuts, fill = type)) + 
    geom_density(alpha = 0.25) + 
    scale_x_log10() + 
    theme_bw() + 
    ggtitle(jmark) + 
    facet_wrap(~type, ncol = 1) + 
    theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}

# Show fold changes for each bin 

# jmark <- "H3K4me1"
# jmark <- "H3K9me3"
jmark <- "H3K4me3"

jmark <- "H3K27me3"
# jbin <- bins.signif.lst[[jmark]][[1]]
jbins <- bins.signif.lst[[jmark]]
print(length(jbins))

jsub.long <- lapply(jmarks, function(jmark){
  jsub <- params.long.lst[[jmark]] %>% filter(bin %in% jbins)
}) %>%
  bind_rows()

ggplot(jsub.long, aes(x = mark, y = abs(log2fc))) + 
  geom_boxplot() + 
  coord_cartesian(ylim = c(0, 5)) + 
  ggtitle(jmark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Check TEs ---------------------------------------------------------------

dat.te <- readRDS("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_gc_analysis/te_counts_50kb_bins.rds")
dat.te$Class2 <- as.factor(dat.te$Class2)

jnames <- names(bins.lst.withrand)
names(jnames) <- jnames

dat.te.filt.long <- lapply(jnames, function(jname){
  jsub <- subset(dat.te, bin %in% bins.lst.withrand[[jname]]) %>%
    group_by(bin, Class2, .drop = FALSE) %>%
    tally()  %>%
    mutate(mark = jname)
}) %>%
  bind_rows() %>%
  group_by(mark, Class2) %>%
  mutate(frac = n / sum(n))

# ggplot(dat.te.filt.long, aes(x = mark, y = frac)) + 
ggplot(dat.te.filt.long, aes(x = mark, y = n)) + 
  geom_point() + 
  # geom_boxplot(outlier.shape = NA) + 
  geom_violin() + 
  facet_wrap(~Class2, ncol = 1, scales = "free_y") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# Check a k9me3 bin  ------------------------------------------------------

# bins.lst.withrand$H3K9me3[1:5]
# bins.lst.withrand$H3K9me3[1]
# 
# (jbin <- bins.lst.withrand$H3K27me3[[1]])
# (jbins <- bins.lst.withrand$H3K9me3[[1:100]])
# (jbins <- bins.lst.withrand$H3K27me3[[1:100]])
# 
# subset(dat.te, bin == jbin)
  
# dim(subset(dat.te, bin == jbin))
jsub <- subset(dat.te, bin %in% jbins)

jsub.sum <- jsub %>%
  group_by(bin, Class2, .drop = FALSE) %>%
  tally()

print(jsub.sum)

ggplot(jsub.sum, aes(x = Class2, y = n)) + 
  geom_point() + 
  geom_violin() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Check diff bins make sense  ---------------------------------------------


dev.off()

