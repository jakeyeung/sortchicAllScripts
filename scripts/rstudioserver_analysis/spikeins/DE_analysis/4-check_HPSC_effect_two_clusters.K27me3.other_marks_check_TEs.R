# Jake Yeung
# Date of Creation: 2021-01-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/DE_analysis/4-check_HPSC_effect_two_clusters.K27me3.other_marks_check_TEs.R
# 



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(ggrastr)
library(JFuncs)
library(scchicFuncs)

library(moments)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)



set.seed(0)
nbootstraps <- 10000
nsampsperbootstrap <- 3000
prob.high <- 0.95
prob.low <- 0.05

# pvalcutoff <- 1e-100
padjcutoff <- 1e-50

# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# jsuffix <- "total"
jsuffix <- "HSPCs_vs_nonHSPCs.spikeins"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned"
# outpdf <- file.path(outdir, paste0("DE_downstream_global_log2_HSPC_vs_nonHSPCs.", jsuffix, ".downstream_genes.", Sys.Date(), ".pdf"))
outpdf <- file.path(outdir, paste0("DE_downstream_global_log2_HSPC_vs_nonHSPCs.", jsuffix, ".downstream_genes.padjcutoff_", padjcutoff, Sys.Date(), ".pdf"))


make.plots <- FALSE
write.tables <- FALSE

if (make.plots){
  pdf(outpdf, useDingbats = FALSE)
}

# Load  -------------------------------------------------------------------

indir.fits <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.HSPCs_vs_nonHSPCs"

out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  # fname <- paste0("poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData")
  fname <- paste0("poisson_fit_bins.", jmark, ".2020-12-29.newannot2.witherrors.MoreBins.", jsuffix, ".RData")
  # fname <- paste0("poisson_fit_bins.", jmark, ".2020-12-29.newannot2.witherrors.MoreBins.HSPCs_vs_nonHSPCs.spikeins.RData")
  inf.fits <- file.path(indir.fits, fname)
  load(inf.fits, v=T)
  params.long <- SummarizeParamsPvalues(jfits.lst, jmark = jmark, paramname = "Cluster") %>%
    mutate(log2fc = estimate / log(2))
  params.long$padj <- p.adjust(params.long$pval.param)
  jnames <- names(jfits.lst); names(jnames) <- jnames
  pvals.long <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  return(list(params.long = params.long, pvals.long = pvals.long))
})

params.lst <- lapply(out.lst, function(out){
  out$params.long
})

pvals.lst <- lapply(out.lst, function(out){
  out$pvals.long
})

# Check skewness across marks  --------------------------------------------

params.sum.long <- lapply(jmarks, function(jmark){
  params.sum <- params.lst[[jmark]] %>%
    mutate(mark = jmark)
}) %>%
  bind_rows()

m <- ggplot(params.sum.long %>% filter(abs(log2fc) < 5), aes(x = log2fc, fill = mark)) + 
  geom_density(alpha = 0.25) + 
  theme_bw(8) + 
  ggtitle(jsuffix) + 
  facet_wrap(~mark, ncol = 1) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)

m <- ggplot(params.sum.long %>% filter(abs(log2fc) < 5), aes(x = log2fc, fill = mark)) + 
  geom_density(alpha = 0.65) + 
  theme_bw(8) + 
  ggtitle(jsuffix) + 
  facet_wrap(~mark, nrow = 1) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m)


# plot mark by mark
for (jmark in jmarks){
  m <- ggplot(params.sum.long %>% filter(abs(log2fc) < 5 & mark == jmark), aes(x = log2fc, fill = mark)) + 
    geom_density(alpha = 0.25) + 
    theme_bw(12) + 
    ggtitle(jmark, jsuffix) + 
    facet_wrap(~mark, ncol = 1) + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}

m <- ggplot(params.sum.long %>% filter(abs(log2fc) < 5), aes(x = log2fc, fill = mark)) + 
  geom_density(alpha = 0.25) + 
  theme_bw(8) + 
  ggtitle(jsuffix) + 
  # facet_wrap(~mark, ncol = 1) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)


# m <- ggplot(params.sum.long %>% filter(abs(log2fc) < 5 & pval.param < 10^-10), aes(x = log2fc, fill = mark)) + 
#   geom_density(alpha = 0.25) + 
#   theme_bw() + 
#   facet_wrap(~mark, ncol = 1) + 
#   geom_vline(xintercept = 0, linetype = "dotted") + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# print(m)

m <- ggplot(params.sum.long %>% filter(abs(log2fc) < 5), aes(x = log2fc, y = -log10(pval.param), color = mark)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  ggtitle(jsuffix) + 
  facet_wrap(~mark, ncol = 1) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)


jsum <- params.lst[["H3K27me3"]] %>%
  group_by(bin) %>%
  mutate(estimate.sign = sign(estimate)) %>%
  group_by(estimate.sign) %>%
  summarise(nbins = length(bin))


# Compare K27me3 vs other marks  ------------------------------------------


refmark <- "H3K27me3"
# comparemark <- "H3K4me1"
comparemark <- "H3K9me3"
comparemarks <- jmarks[jmarks != refmark]

for (comparemark in comparemarks){
  
  params1 <- params.lst[[refmark]] %>%
    filter(abs(log2fc) < 5) %>%
    dplyr::rename(log2fc1 = log2fc)
  params2 <- params.lst[[comparemark]] %>%
    filter(abs(log2fc) < 5) %>%
    dplyr::rename(log2fc2 = log2fc)
  
  params.merge <- left_join(params1, params2, by = c("bin"))
  
  
  m <- ggplot(params.merge, aes(x = log2fc1, y = log2fc2)) + 
    geom_point_rast(alpha = 0.25) + 
    geom_density_2d() + 
    ggtitle(jsuffix) + 
    theme_bw() + 
    xlab(refmark) + 
    ylab(comparemark) + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  print(m)
  
  
}



# Calculate skews  --------------------------------------------------------

params.filt.lst <- lapply(params.lst, function(x) subset(x, abs(estimate) < 5))

ivec <- seq(nbootstraps)
names(ivec) <- ivec

skews.lst.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  skews.lst <- lapply(ivec, function(i){
    xsamp <- sample(x = params.filt.lst[[jmark]]$estimate, size = nsampsperbootstrap, replace = TRUE)
    skew.tmp <- moments::skewness(xsamp)
  })
})

dat.skews.lst <- lapply(jmarks, function(jmark){
  skews.vec <- unlist(skews.lst.lst[[jmark]])
  dat.skews.tmp <- data.frame(mark = jmark, skew = skews.vec, i = seq(length(skews.vec)), stringsAsFactors = FALSE)
}) 

# cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
ggplot(dat.skews.lst %>% bind_rows(), aes(x = skew, fill = mark)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.skews.lst %>% bind_rows(), aes(x = skew, fill = mark)) + 
  geom_density(alpha = 0.25) + 
  facet_wrap(~mark, nrow = 1) + 
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# get skewness as mean and 95 % esstimate
# get upper values
dat.mean.ci <- lapply(dat.skews.lst, function(jdat){
  jmean <- mean(jdat$skew)
  jupper <- quantile(jdat$skew, probs = prob.high)
  jlower <- quantile(jdat$skew, probs = prob.low)
  jmark <- jdat$mark[[1]]
  jdat.out <- data.frame(skew.mean = jmean, skew.upper = jupper, skew.lower = jlower, mark = jmark, stringsAsFactors = FALSE)
}) %>%
  bind_rows()

# rearrange
dat.mean.ci$mark <- factor(dat.mean.ci$mark, levels = c("H3K4me3", "H3K27me3", "H3K4me1", "H3K9me3"))

ggplot(dat.mean.ci, aes(x = mark, y = skew.mean, ymin = skew.lower, ymax = skew.upper)) + 
  geom_errorbar() + 
  geom_point() + 
  theme_bw(18) + 
  xlab("") + ylab("Skew") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# Check the bins  ---------------------------------------------------------

head(params.lst$H3K4me1)




jmarkref <- "H3K27me3"
# jmarkcompare <- "H3K4me1"
# jmarkcompare <- "H3K4me3"
jmarkcompare <- "H3K9me3"

# jfilt <- params.lst[[jmarkref]] %>%
#   arrange(estimate) %>%
#   filter(abs(estimate) < 5 & pval.param < pvalcutoff)

jfilt <- params.lst[[jmarkref]] %>%
  arrange(estimate) %>%
  filter(abs(estimate) < 5 & pval.param < padjcutoff)

ggplot(params.lst$H3K27me3, aes(x = estimate, y = -log10(pval.param))) +
# ggplot(params.lst$H3K27me3, aes(x = estimate, y = -log10(padj))) + 
  geom_point(alpha = 0.05) + 
  theme_bw() + 
  coord_cartesian(xlim = c(-5, 5)) + 
  geom_hline(yintercept = -log10(padjcutoff)) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jbins.keep <- jfilt$bin

jfilt2 <- params.lst[[jmarkcompare]] %>%
  filter(bin %in% jbins.keep)

jfilt.merge <- left_join(jfilt, jfilt2, by = c("bin"))

ggplot(jfilt.merge, aes(x = log2fc.x, y = log2fc.y)) + 
  geom_point() + 
  ggtitle(paste(jmarkref, "vs", jmarkcompare)) + 
  theme_bw() +  
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_density_2d() + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# assign bin to nearest gene

inf.tsspretty <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.50000.again.cut.pretty.bed"
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
bins.annot.tmp <- AnnotateCoordsFromList(coords.vec = jbins.keep, inf.tss = inf.tsspretty, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

# get closest 
jgenes.keep <- bins.annot.tmp$out2.df.closest$gene


# Check Giladi 
inf.public <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/giladi_pseudobulk_exprs_data.rds"
dat.public <- readRDS(inf.public)

dat.public$gene <- as.character(dat.public$gene)


m2c <- MarkerToCelltype()
dat.public$celltype2 <- sapply(as.character(dat.public$celltype), function(x) m2c[[x]])

jgenes.keep.str <- paste(jgenes.keep, collapse = "|")

# check what's the distribution of exprs across celltypes 

dat.public.sum <- dat.public %>%
  group_by(gene) %>%
  summarise(exprs = mean(exprs)) %>%
  mutate()

ggplot(dat.public.sum, aes(x = exprs)) + 
  geom_density() +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dat.public.sum$in.set <- sapply(as.character(dat.public.sum$gene), function(g){
  jsplit <- strsplit(g, ";")[[1]]
  if (length(jsplit) == 1){
    jcheck <- g %in% jgenes.keep
  } else if (length(jsplit) > 1) {
    for (i in seq(length(jsplit))){
      jcheck <- jsplit[[i]] %in% jgenes.keep
      if (jcheck == TRUE){
        return(jcheck)
      }
    }
  }
  return(jcheck)
})


jgenes.keep.in.giladi <- subset(dat.public.sum, in.set)$gene

# dat.public.filt <- subset(dat.public, grepl(jgenes.keep.str, as.character(gene)))

dat.public.filt <- subset(dat.public, gene %in% jgenes.keep.in.giladi)
dat.public.filt$gene <- as.character(dat.public.filt$gene)

jgenes.keep.grepout <- unique(dat.public.filt$gene)

ggplot(dat.public.filt, aes(x = celltype, y = exprs)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# assign filtered genes and plot density with two groups
# dat.public.sum$in.set <- sapply(dat.public.sum$gene, function(g) g %in% jgenes.keep.grepout)

# ggplot(dat.public.sum, aes(x = celltype, ))

ggplot(dat.public.sum, aes(x = exprs, fill = in.set)) + 
  geom_density(alpha = 0.5) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jsum <- dat.public.sum %>%
  group_by(in.set) %>%
  summarise(ngenes = length(gene))

jtitle <- paste(paste(jsum$in.set, jsum$ngenes, sep = "="), paste0("Pval=", padjcutoff), collapse = ", ")

ggplot(dat.public.sum, aes(x = in.set, y = exprs)) + 
  geom_boxplot() + 
  ggtitle(jtitle) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load cluster changes  ---------------------------------------------------

# check k9me3 50kb changes


inf.k9me3 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.allmarks_spikeins/poisson_fit_bins.H3K9me3.2020-12-26.newannot2.witherrors.MoreBins.newestmeta.RData"
load(inf.k9me3, v=T)

params.long <- SummarizeParamsPvalues(jfits.lst, jmark = jmark, paramname = "Cluster") %>%
  mutate(log2fc = estimate / log(2))
params.long$padj <- p.adjust(params.long$pval.param)
jnames <- names(jfits.lst); names(jnames) <- jnames
pvals.long <- lapply(jnames, function(jname){
  x <- jfits.lst[[jname]]
  xvec <- x$pval
  data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
}) %>%
  bind_rows()

# filter bins
jsub <- subset(params.long, bin %in% jbins.keep)
# jsub <- params.long

ggplot(jsub %>% filter(abs(estimate) < 5), aes(x = estimate, y = -log10(pval.param), color = param)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  facet_wrap(~param) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

inf.k27me3 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.allmarks_spikeins/poisson_fit_bins.H3K27me3.2020-12-25.newannot2.witherrors.MoreBins.newestmeta.RData"
load(inf.k27me3, v=T)

jfits.lst.k27me3 <- jfits.lst
jmat.mark.k27me3 <- jmat.mark
dat.annots.filt.mark.k27me3 <- dat.annots.filt.mark
ncuts.for.fit.mark.k27me3 <- ncuts.for.fit.mark

params.long.k27me3 <- SummarizeParamsPvalues(jfits.lst.k27me3, jmark = jmark, paramname = "Cluster") %>%
  mutate(log2fc = estimate / log(2))
params.long.k27me3$padj <- p.adjust(params.long.k27me3$pval.param)
jnames.k27me3 <- names(jfits.lst); names(jnames.k27me3) <- jnames.k27me3
pvals.long.k27me3 <- lapply(jnames.k27me3, function(jname){
  x <- jfits.lst.k27me3[[jname]]
  xvec <- x$pval
  data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
}) %>%
  bind_rows()

# filter bins
jsub.k27me3 <- subset(params.long.k27me3, bin %in% jbins.keep)
# jsub <- params.long

ggplot(jsub.k27me3 %>% filter(abs(estimate) < 5), aes(x = estimate, y = -log10(pval.param), color = param)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  facet_wrap(~param) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# combine k27me3 with k9me3 
jsub2.k9me3 <- subset(jsub, param == "ClusterEryths.Estimate")
jsub2.k27me3 <- subset(jsub.k27me3, param == "ClusterEryths.Estimate")

jsub2.merge <- left_join(jsub2.k27me3, jsub2.k9me3, by = c("bin", "param"))

ggplot(jsub2.k9me3 %>% filter(abs(estimate) < 5), aes(x = estimate)) + geom_density() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jsub2.k9me3 %>% filter(abs(estimate) < 5), aes(y = estimate)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jsub2.k9me3 %>% filter(abs(estimate) < 5), aes(x = estimate, y = -log10(pval.param))) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(jsub2.merge %>% filter(abs(estimate.x) < 5 & abs(estimate.y) < 5), aes(x = estimate.x, y = estimate.y)) + 
  geom_point() + 
  geom_density_2d() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Check TEs in bins  ------------------------------------------------------

# load TEs 

jval.cutoff <- 2000
# read TEs
hubprefix <- "/home/jyeung/hub_oudenaarden"

inf.tes <- file.path(hubprefix, "jyeung/data/databases/repeatmaskers/mm10/mm10_repeatmasker_track_viz.gz")
dat.tes <- fread(inf.tes, header = FALSE, col.names = c("Chr", "Start", "End", "Name", "Value", "Strand"))

dat.tes$Class <- sapply(dat.tes$Name, function(x) strsplit(x, split = "#")[[1]][[2]])

dat.tes.filt <- subset(dat.tes, !grepl("\\?", Class)) %>%
  rowwise() %>%
  mutate(Class2 = strsplit(Class, "\\/")[[1]][[1]]) %>%
  ungroup() %>%
  filter(Class2 %in% c("LINE", "LTR", "SINE", "DNA"),
         Value > jval.cutoff)

set.seed(0)
jbins.background <- sample(x = unique(params.lst[[jmarkref]]$bin), size = length(jbins.keep), replace = FALSE)

dat.bins <- data.frame(Chr = sapply(jbins.keep, GetChromo), 
                       Start = sapply(jbins.keep, GetStart),
                       End = sapply(jbins.keep, GetEnd), 
                       Name = jbins.keep, 
                       stringsAsFactors = FALSE)

dat.bins.bg <- data.frame(Chr = sapply(jbins.background, GetChromo), 
                       Start = sapply(jbins.background, GetStart),
                       End = sapply(jbins.background, GetEnd), 
                       Name = jbins.background, 
                       stringsAsFactors = FALSE) %>%
  filter(Chr %in% jchromos)
print(unique(dat.bins.bg$Chr))

dat.bins.bg <- dat.bins.bg[gtools::mixedorder(dat.bins.bg$Name), ]

tes.gr <- makeGRangesFromDataFrame(dat.tes.filt, keep.extra.columns = TRUE)
bins.gr <- makeGRangesFromDataFrame(dat.bins, keep.extra.columns = TRUE)
binsbg.gr <- makeGRangesFromDataFrame(dat.bins.bg, keep.extra.columns = TRUE)

# gr.out <- findOverlaps(bins.gr, tes.gr, type = "within")
gr.out <- findOverlaps(tes.gr, bins.gr, type = "within")
grbg.out <- findOverlaps(tes.gr, binsbg.gr, type = "within")
# dat.overlap <- data.frame(bins.gr[queryHits(gr.out), ], tes.gr[subjectHits(gr.out), ], stringsAsFactors = FALSE)
dat.overlap <- data.frame(bins.gr[subjectHits(gr.out), ], tes.gr[queryHits(gr.out), ], stringsAsFactors = FALSE)
dat.overlap.bg <- data.frame(binsbg.gr[subjectHits(grbg.out), ], tes.gr[queryHits(grbg.out), ], stringsAsFactors = FALSE)

# count TEs

dat.overlap.sum <- dat.overlap %>%
  group_by(Name) %>%
  summarise(ntes = length(Name.1)) %>%
  mutate(Type = "HSPC-specific")

dat.overlap.bg.sum <- dat.overlap.bg %>%
  group_by(Name) %>%
  summarise(ntes = length(Name.1)) %>%
  mutate(Type = "Random")

dat.overlap.merge <- rbind(dat.overlap.sum, dat.overlap.bg.sum)

ggplot(dat.overlap.merge, aes(x = ntes, fill = Type)) + 
  geom_density(alpha = 0.5) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


if (make.plots){
  dev.off()
}



# Write bins to output -----------------------------------------------------


# Get H3K27me3 bins
jbins.keep <- jfilt$bin

cnames.check <- c("bin", "ClusternotHSPCs.Estimate")
params.signifs.wide.lst <- lapply(jmarks, function(jmark){
  jparam <- params.lst[[jmark]]
  jparam.wide <- data.table::dcast(subset(jparam, bin %in% jbins.keep), formula = bin ~ param, value.var = "log2fc")
  assertthat::assert_that(all(colnames(jparam.wide) %in% cnames.check))
  cnames.new <- c("bin", paste0("log2fc.ClusternotHSPCs.Estimate.", jmark))
  colnames(jparam.wide) <- cnames.new
  return(jparam.wide)
})

pvals.signifs.wide.lst <- lapply(jmarks, function(jmark){
  jparam <- params.lst[[jmark]]
  jparam.wide <- data.table::dcast(subset(jparam, bin %in% jbins.keep), formula = bin ~ param, value.var = "pval.param")
  assertthat::assert_that(all(colnames(jparam.wide) %in% cnames.check))
  cnames.new <- c("bin", paste0("pval.ClusternotHSPCs.Estimate.", jmark))
  colnames(jparam.wide) <- cnames.new
  return(jparam.wide)
})

# join
params.signifs.wide.joined <- Reduce(left_join, params.signifs.wide.lst)
pvals.signifs.wide.joined <- Reduce(left_join, pvals.signifs.wide.lst)
params.pvals.joined <- left_join(params.signifs.wide.joined, pvals.signifs.wide.joined)


# Add nearest genes -------------------------------------------------------


dat.bins.filt <- data.frame(Chr = sapply(jbins.keep, JFuncs::GetChromo),
                            Start = sapply(jbins.keep, JFuncs::GetStart, returnAsInt = TRUE),
                            End = sapply(jbins.keep, JFuncs::GetEnd, returnAsInt = TRUE), 
                            Name = jbins.keep,
                            stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(Midpt = mean(c(Start, End)))

# fourth column needs to be gene, fifth column needs to be distance
bins.mid <- paste(dat.bins.filt$Chr, paste(dat.bins.filt$Midpt - 1, dat.bins.filt$Midpt + 1, sep = "-"), sep = ":")
bins.orig <- dat.bins.filt$Name

bins.hash <- hash::hash(bins.mid, bins.orig)

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
bins.uniq <- unique(bins.mid)
dat.annot <- AnnotateCoordsFromList(coords.vec = bins.uniq, inf.tss = "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.100000.bed", txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

dat.annot2 <- dat.annot$regions.annotated

dat.annot2$bin <- sapply(bins.orig, function(x) AssignHash(x = x, jhash = bins.hash, null.fill = x))

params.pvals.joined.annot <- left_join(params.pvals.joined, dat.annot2)





# Write output ------------------------------------------------------------

if (write.tables){
  outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/tables_H3K27me3_HSPC_lost_bins"
  fname <- paste0("H3K27me3_HSPCs_lost_bins_DE_estimates.", Sys.Date(), ".txt")
  outf <- file.path(outdir, fname)
  fwrite(params.pvals.joined.annot, file = outf, sep = "\t")
}

# 
# # write background K27me3? 
# 
# dat.pbulk.lst <- lapply(jmarks, function(jmark.check){
#   print(jmark.check)
#   cnames.keep.lst <- split(x = dat.metas[[jmark.check]]$cell, f = dat.metas[[jmark.check]]$cluster)
#   pbulks.lst <- SumAcrossClusters(outs.all.lst[[jmark.check]]$count.mat, cnames.keep.lst = cnames.keep.lst)
#   
#   
#   mat.pbulk <- bind_rows(pbulks.lst) %>%
#     as.data.frame()
#   rownames(mat.pbulk) <- rownames(outs.all.lst[[jmark.check]]$count.mat)
#   mat.pbulk <- sweep(mat.pbulk, MARGIN = 2, STATS = colSums(mat.pbulk), FUN = "/")
#   
#   dat.pbulk <- as.matrix(mat.pbulk) %>%
#     melt()
#   colnames(dat.pbulk) <- c("bin", "ctype", "cuts")
#   dat.pbulk$mark <- jmark.check
#   return(dat.pbulk)
# })
# 
# 
# 
# 
# ggplot(params.lst$H3K27me3, aes(x = estimate)) + 
#   geom_density() + 
#   coord_cartesian(xlim = c(-5, 5)) + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



