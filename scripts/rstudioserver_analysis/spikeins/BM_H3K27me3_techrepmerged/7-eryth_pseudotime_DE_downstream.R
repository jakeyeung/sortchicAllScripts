# Jake Yeung
# Date of Creation: 2020-12-26
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_techrepmerged/7-eryth_pseudotime_DE_downstream.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)



# Load DE  ----------------------------------------------------------------


# inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.allmarks_spikeins/poisson_fit_bins.H3K27me3.2020-12-26.newannot2.witherrors.MoreBins.newestmeta.totalcuts.pseudotime.RData"
inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.allmarks_spikeins/poisson_fit_bins.H3K27me3.2020-12-26.newannot2.witherrors.MoreBins.newestmeta.spikeins.pseudotime.RData"
load(inf.de, v=T)
jfits.lst.spikeins.ptime <- jfits.lst
params.dat.spikeins.ptime.withpval <- SummarizeParamsPvalues(jfits.lst.spikeins.ptime, jmark = "H3K27me3", paramname = "Pseudotime")

# inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.allmarks_spikeins/poisson_fit_bins.H3K27me3.2020-12-26.newannot2.witherrors.MoreBins.newestmet"
inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.allmarks_spikeins/poisson_fit_bins.H3K27me3.2020-12-25.newannot2.witherrors.MoreBins.newestmeta.RData"
load(inf.de, v=T)
jfits.lst.spikeins <- jfits.lst
params.dat.spikeins.withpval <- SummarizeParamsPvalues(jfits.lst.spikeins, jmark = "H3K27me3", paramname = "Cluster")


inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.H3K27me3.2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData"
load(inf.de, v=T)
print(head(ncuts.for.fit.mark))
jfits.lst.total <- jfits.lst
params.dat.total.withpval <- SummarizeParamsPvalues(jfits.lst.total, jmark = "H3K27me3", paramname = "Cluster")


# inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.H3K4me1.2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData"
inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.allmarks_spikeins/poisson_fit_bins.H3K4me1.2020-12-25.newannot2.witherrors.MoreBins.newestmeta.RData"
load(inf.de, v=T)
jfits.lst.spikeins.k4me1 <- jfits.lst
params.dat.spikeins.withpval.k4me1 <- SummarizeParamsPvalues(jfits.lst.spikeins.k4me1, jmark = "H3K4me1", paramname = "Cluster")



ggplot(params.dat.spikeins.withpval, aes(x = estimate, y = -log10(pval.param))) + 
  geom_point() + 
  facet_wrap(~param) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(params.dat.total.withpval, aes(x = estimate, y = -log10(pval.param))) + 
  geom_point() + 
  facet_wrap(~param) + 
  theme_bw() + 
  coord_cartesian(xlim = c(-5, 5)) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(params.dat.total.withpval, aes(x = estimate, fill = param)) +
  geom_density(alpha = 0.25) + 
  facet_wrap(~param) + 
  theme_bw() + 
  coord_cartesian(xlim = c(-5, 5)) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Compare total vs spikeins -----------------------------------------------

params.dat.spikeins.withpval %>% 
  arrange(pval.param)

jbin <- "chr7:112700000-112750000"  # tead1
jbin <- "chr15:103000000-103050000"

jtest <- subset(params.dat.total.withpval, bin == jbin)
jtest2 <- subset(params.dat.spikeins.withpval, bin == jbin)
print(jtest)
print(jtest2)

params.merge <- left_join(params.dat.total.withpval %>% dplyr::rename(estimate.total = estimate, pval.param.total = pval.param) %>% 
                            dplyr::select(c(bin, param, estimate.total, pval.param.total)),
                         params.dat.spikeins.withpval %>% dplyr::rename(estimate.spikeins = estimate, pval.param.spikeins = pval.param) %>% 
                           dplyr::select(c(bin, param, estimate.spikeins, pval.param.spikeins)))

params.merge.filt <- params.merge %>% filter(pval.param.total < 10^-10 & pval.param.spikeins < 10^-10)
ggplot(params.merge.filt, aes(x = estimate.total, y = estimate.spikeins)) + 
  facet_wrap(~param) + 
  geom_point() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Get global loss HSPC to all  --------------------------------------------


head(params.dat.total.withpval)

params.dat.total.withpval.sum <- params.dat.total.withpval %>%
  group_by(bin) %>%
  summarise(estimate.mean = (1 / length(estimate)) * sum(estimate),
            se.mean = (1 / length(se)) * sqrt(sum(se^2))) %>%
  rowwise() %>%
  mutate(z = estimate.mean / se.mean,
         pval.param.mean = exp(-0.717 * z - 0.416 * z ^ 2)) %>%
  arrange(pval.param.mean)

print(params.dat.total.withpval.sum %>% filter(estimate.mean < 0))


# Annotate bins  ----------------------------------------------------------

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
bins.annot.total.loss <- AnnotateCoordsFromList.GeneWise(coords.vec = params.dat.total.withpval.sum$bin, 
                                                         inf.tss = "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.100000.bed", 
                                                         txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

params.dat.total.withpval.sum.annot <- left_join(params.dat.total.withpval.sum, bins.annot.total.loss$out2.df %>% dplyr::select(region_coord, gene, dist.to.tss), by = c("bin" = "region_coord"))

params.dat.total.withpval.sum.annot.filt <- params.dat.total.withpval.sum.annot %>% 
  filter(pval.param.mean < 10^-10) %>%
  arrange(estimate.mean) %>%
  group_by(bin) %>%
  filter(abs(dist.to.tss) == min(abs(dist.to.tss)))

# Get pseudotime  ---------------------------------------------------------

ggplot(params.dat.spikeins.ptime.withpval %>% filter(abs(estimate) < 1), aes(x = estimate, y = -log10(pval.param))) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


head(params.dat.spikeins.ptime.withpval %>% arrange(pval.param))

# AnnotateCoordsFromList(params.dat.spikeins.ptime.withpval$bin, )

bins.annot <- AnnotateCoordsFromList.GeneWise(coords.vec = params.dat.spikeins.ptime.withpval$bin, inf.tss = "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.100000.bed", txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

params.dat.spikeins.ptime.withpval.annot <- left_join(params.dat.spikeins.ptime.withpval, subset(bins.annot$out2.df, select = c(region_coord, gene, dist.to.tss)), by = c("bin" = "region_coord")) %>%
  arrange(pval.param)
  
params.dat.spikeins.ptime.withpval.annot.sum <- params.dat.spikeins.ptime.withpval.annot %>%
  group_by(bin) %>%
  filter(abs(dist.to.tss) == min(abs(dist.to.tss)))

ggplot(params.dat.spikeins.ptime.withpval.annot.sum %>% filter(abs(estimate) < 1), aes(x = estimate))  + 
  geom_density() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(params.dat.spikeins.ptime.withpval.annot.sum %>% filter(abs(estimate) < 1), aes(x = estimate, y = -log10(pval.param)))  + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# check raw
jgene <- "Hoxb2"
jgene <- "Hbb-bs"
jgene <- "Hoxb2"
jgene <- "Hbb-bs"

jgene <- "Pax6"
jgene <- "Sox6"
jgene <- "Sox6"
jgene <- "Hoxc4"
jgene <- "Zic2"
jgene <- "Hoxa9"
jgene <- "Sfi1"
jgene <- "Hbb-bs"
jgene <- "Hoxc4"
(jbin <- subset(params.dat.spikeins.ptime.withpval.annot.sum, gene == jgene)$bin[[1]])


cuts.raw <- data.frame(cell = colnames(jmat.mark), cuts = jmat.mark[jbin, ], stringsAsFactors = FALSE) %>%
  left_join(., ncuts.for.fit.mark) %>%
  left_join(., dat.annots.filt.mark)

ggplot(cuts.raw, aes(x = Pseudotime, y = log2(cuts / ncuts.total))) + 
  geom_point() + 
  ggtitle(jgene, jbin) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(cuts.raw, aes(x = umap1, y = umap2, color = log2(cuts / ncuts.total))) + 
  geom_point() + 
  ggtitle(jgene, jbin) + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(cuts.raw, aes(x = umap1, y = umap2, color = log2(cuts_in_peak / ncuts.total))) + 
  geom_point() + 
  ggtitle(jgene, jbin) + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(cuts.raw, aes(x = Pseudotime, y = log2(cuts / ncuts.total))) + 
  geom_point() + 
  ggtitle(jgene, jbin) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(cuts.raw, aes(x = Pseudotime, y = log2(cuts_in_peak / ncuts.total))) + 
  geom_point() + 
  ggtitle(jgene, jbin) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Check active mark -------------------------------------------------------


inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.H3K4me1.2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData"
load(inf, v=T)

jmark <- "H3K4me1"
jnames <- names(jfits.lst); names(jnames) <- jnames

# https://stats.stackexchange.com/questions/315311/how-to-find-p-value-using-estimate-and-standard-error
params.dat.all.withpval.k4me1 <- lapply(jnames, function(jname){
  x <- jfits.lst[[jname]]
  xkeep <- grepl("^Cluster.*.Estimate$", x = names(x))
  jparams <- x[xkeep]
  xkeep.se <- grepl("^Cluster.*.StdError$", x = names(x))
  jparams.se <- x[xkeep.se]
  jout <- data.frame(bin = jname, param = names(jparams), estimate = unlist(jparams), se = unlist(jparams.se), mark = jmark, stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(z = estimate / se, 
           pval.param = exp(-0.717 * z - 0.416 * z ^ 2))
}) %>%
  bind_rows() 

print(params.dat.all.withpval.annot.sum)

jbin <- "chr15:103000000-103050000"
jbin <- "chr15:102950000-103000000"
jbin <- "chr15:102900000-102950000"
jbin <- "chr2:74700000-74750000"
jbin <- "chr14:122450000-122500000"
jbin <- "chr17:85600000-85650000"
jbin <- "chr11:96250000-96300000"
jbin <- "chr7:79500000-79550000"
jbin <- "chr7:103800000-103850000"

# top eryth k4me1 high
jbin <- "chr7:90050000-90100000"
jbin <- "chr12:76650000-76700000"

subset(params.dat.all.withpval.k4me1, bin == jbin)
subset(params.dat.all.withpval, bin == jbin)


jsub <- subset(params.dat.all.withpval.k4me1, param == "ClusterEryths.Estimate") %>%
  arrange(pval.param)

print(jsub)
 


 
# Giladi ------------------------------------------------------------------

inf.public <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/giladi_pseudobulk_exprs_data.rds"
dat.public <- readRDS(inf.public)

m2c <- MarkerToCelltype()
dat.public$celltype2 <- sapply(as.character(dat.public$celltype), function(x) m2c[[x]])

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

dat.public$gene <- as.character(dat.public$gene)
jsub.pub <- subset(dat.public, grepl("Hoxc6", gene))
jsub.pub <- subset(dat.public, grepl("Hox", gene))

ggplot(jsub.pub, aes(x = celltype2, y = zscore)) + geom_boxplot() + geom_point()

  # bins.uniq <- unique(bins.keep)
# dat.annot <- AnnotateCoordsFromList(coords.vec = bins.uniq, inf.tss = "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.100000.bed", txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

