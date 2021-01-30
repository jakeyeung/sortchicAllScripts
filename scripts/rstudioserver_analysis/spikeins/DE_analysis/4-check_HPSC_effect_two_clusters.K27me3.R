# Jake Yeung
# Date of Creation: 2021-01-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/DE_analysis/4-check_HPSC_effect_two_clusters.K27me3.R
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


if (make.plots){
  dev.off()
}



