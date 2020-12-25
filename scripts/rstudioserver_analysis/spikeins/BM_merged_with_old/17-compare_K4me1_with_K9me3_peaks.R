# Jake Yeung
# Date of Creation: 2020-12-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/17-compare_K4me1_with_K9me3_peaks.R
# Look at k4me1 peaks to find anticorrelation

rm(list=ls())
   
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(hash)
library(igraph)
library(umap)

library(JFuncs)

library(DescTools)


GetCutFracsFromBinSet <- function(bset){
  m1 <- ggplot(bset.dat, aes(x = param , y = estimate2)) + 
    ggtitle(paste("Nbins:", length(bset))) + 
    geom_boxplot() 
  
  m2 <- ggplot(params.dat.merge %>% filter(bin %in% bset), aes(x = estimate1, y = estimate2)) + 
    geom_point(alpha = 0.25, color = "grey75")  + 
    geom_density_2d() + 
    facet_wrap(~param) + 
    ggtitle(paste("Nbins:", length(bset))) + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  print(m1)
  print(m2)
  
  # show on UMAP 
  cellsums1 <- colSums(mats1)
  cellsums2 <- colSums(mats2)
  
  cutfrac1 <- data.frame(cell = colnames(mats1), cutfrac = colSums(sweep(mats1[bset, ], MARGIN = 2, STATS = cellsums1, FUN = "/")))
  cutfrac2 <- data.frame(cell = colnames(mats2), cutfrac = colSums(sweep(mats2[bset, ], MARGIN = 2, STATS = cellsums2, FUN = "/")))
  
  # add to UMAP 
  dat.annot1.merge <- left_join(dat.annot1, cutfrac1)
  dat.annot2.merge <- left_join(dat.annot2, cutfrac2)
  return(list(dat.annot1.merge = dat.annot1.merge, dat.annot2.merge = dat.annot2.merge))
}

# Load UMAPs --------------------------------------------------------------
hubprefix <- "/home/jyeung/hub_oudenaarden"
jmark1 <- "H3K4me1"
jmark2 <- "H3K9me3"

outpdf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/anticorrelations_", jmark1, "_vs_", jmark2))




inf.annot1 <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.", jmark1, ".txt"))
# inf.annot2 <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins/cell_cluster_table_with_spikeins.", jmark2, ".2020-11-18.dupfilt.txt"))
inf.annot2 <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/BM_celltypes.bincutoff_0.binskeep_1000.byplate.szname_none.niter_500.reorder_rownames.dupfilt.2020-11-23.H3K9me3.txt"))

dat.annot1 <- fread(inf.annot1)
dat.annot2 <- fread(inf.annot2)


# Load count tables -------------------------------------------------------

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/count_tables_from_peaks.from_sitecount_mat.H3K4me1_peaks"
inf.counts1 <- list.files(indir, pattern = paste0("BM_round1_round2_merged_", jmark1, "_"), full.names = TRUE)
inf.counts2 <- list.files(indir, pattern = paste0("BM_round1_round2_merged_", jmark2, "_"), full.names = TRUE)


mats1 <- lapply(inf.counts1, function(inf){
  mat <- ReadMatTSSFormat(inf)
  mat <- mat[!duplicated(rownames(mat)), ]
})

mats2 <- lapply(inf.counts2, function(inf){
  mat <- ReadMatTSSFormat(inf)
  mat <- mat[!duplicated(rownames(mat)), ]
})
all.rnames1 <- unique(unlist(lapply(mats1, rownames)))
all.rnames2 <- unique(unlist(lapply(mats2, rownames)))
all.rnames <- unique(intersect(all.rnames1, all.rnames2))

mats1 <- JFuncs::cbind.fill.lst(mats1, all.rnames, fill = 0)
print(dim(mats1))

mats2 <- JFuncs::cbind.fill.lst(mats2, all.rnames, fill = 0)
print(dim(mats2))


# Do pseudobulk -----------------------------------------------------------

plot(density(log10(colSums(mats1))))
plot(density(log10(colSums(mats2))))

cells.bycluster1 <- split(f = dat.annot1$cluster, x = dat.annot1$cell)
cells.bycluster2 <- split(f = dat.annot2$cluster, x = dat.annot2$cell)

mats1.pbulk <- SumAcrossClusters(mats1, cells.bycluster1)
mats1.pbulk <- do.call(cbind, mats1.pbulk)
mats1.pbulk <- sweep(mats1.pbulk, MARGIN = 2, STATS = colSums(mats1.pbulk), FUN = "/")

mats2.pbulk <- SumAcrossClusters(mats2, cells.bycluster2)
mats2.pbulk <- do.call(cbind, mats2.pbulk)
mats2.pbulk <- sweep(mats2.pbulk, MARGIN = 2, STATS = colSums(mats2.pbulk), FUN = "/")
# 
# # Compare corrlations -----------------------------------------------------
# 
# plot(log2(mats1.pbulk[, "Bcells"]), log2(mats2.pbulk[, "Lymphoid"]), pch = 20)
# plot(mats1.pbulk[, "Bcells"], mats2.pbulk[, "Lymphoid"], pch = 20)
# 
# plot(log2(mats1.pbulk[, "Bcells"]), log2(mats2.pbulk[, "Granulocytes"]), pch = 20)
# plot(mats1.pbulk[, "Bcells"], mats2.pbulk[, "Granulocytes"], pch = 20)
# 
# plot(log2(mats1.pbulk[, "Bcells"]), log2(mats2.pbulk[, "HSPCs"]), pch = 20)
# plot(mats1.pbulk[, "Bcells"], mats2.pbulk[, "HSPCs"], pch = 20)
#   
# plot(mats1.pbulk[, "HSPCs"], mats2.pbulk[, "HSPCs"], pch = 20)
# 
# plot(mats1.pbulk[, "Eryths"], mats2.pbulk[, "Eryth"], pch = 20)
# 
# plot(mats1.pbulk[, "Bcells"], mats2.pbulk[, "Lymphoid"], pch = 20, xlab = "K4me1_Bcells", ylab = "K9me3_Bcells")
# 
# ggplot(data = data.frame(x = mats1.pbulk[, "Bcells"], y = mats2.pbulk[, "Lymphoid"]), aes(x = x, y = y)) + 
#   geom_point(alpha = 0.25) +  
#   geom_density_2d(bins = 100) + 
#   scale_x_log10() + scale_y_log10() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# plot(mats1.pbulk[, "Eryths"], mats2.pbulk[, "Lymphoid"], pch = 20)
# plot(mats1.pbulk[, "Granulocytes"], mats2.pbulk[, "Lymphoid"], pch = 20)
# 

# Get diff genes ----------------------------------------------------------


inf.de1 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3/poisson_fit_hiddendomains.H3K4me1_peaks.", jmark1, ".2020-12-03.newannot2.witherrors.RData")
inf.de2 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3/poisson_fit_hiddendomains.H3K4me1_peaks.", jmark2, ".2020-12-04.newannot2.witherrors.RData")

load(inf.de1, v=T)

jfits.lst1 <- jfits.lst
params.lst1 <- lapply(jfits.lst1, function(x){
  xkeep <- grepl("^Cluster.*.Estimate$", x = names(x))
  x[xkeep]
})

pvals.lst1 <- lapply(jfits.lst1, function(x) x$pval)

load(inf.de2, v=T)

jfits.lst2 <- jfits.lst
params.lst2 <- lapply(jfits.lst2, function(x){
  xkeep <- grepl("^Cluster.*.Estimate$", x = names(x))
  x[xkeep]
})
# jnames2 <- names(jfits.lst2)
# names(jnames2) <- jnames2
# params.dat2 <- lapply(jnames2, f)
pvals.lst2 <- lapply(jfits.lst2, function(x) x$pval)

plot(density(-log10(unlist(pvals.lst1))))
plot(density(-log10(unlist(pvals.lst2))))

# find diff genes in k9me3
k9.bins <- which(pvals.lst2 < 1e-10)

bins.keep <- unique(names(pvals.lst2)[k9.bins])

jnames <- bins.keep
names(jnames) <- jnames
params.dat1 <- lapply(jnames, function(jname){
  jparams <- params.lst1[[jname]]
  data.frame(bin = jname, param = names(jparams), estimate1 = unlist(jparams), stringsAsFactors = FALSE)
}) %>%
  bind_rows()

params.dat2 <- lapply(jnames, function(jname){
  jparams <- params.lst2[[jname]]
  data.frame(bin = jname, param = names(jparams), estimate2 = unlist(jparams), stringsAsFactors = FALSE)
}) %>%
  bind_rows()

unique(params.dat2$param)

# match parm names then merge
# rename eryth to eryths, lymphoid to bcells

params.dat2 <- params.dat2 %>%
  mutate(param = gsub("Eryth", "Eryths", param),
         param = gsub("Lymphoid", "Bcells", param))

# # scramble?
# params.dat2 <- params.dat2 %>%
#   rowwise() %>%
#   mutate(param = gsub("Eryth", "NeutroPlaceholder", param),
#          param = gsub("Lymphoid", "Eryths", param),
#          param = gsub("Granulocytes", "Bcells", param),
#          param = gsub("NeutroPlaceholder", "Granulocytes", param))

(params.keep <- unique(params.dat2$param))

params.dat1.sub <- subset(params.dat1, param %in% params.keep)

# merge
params.dat.merge <- left_join(params.dat1.sub, params.dat2)

jbin.keep <- bins.keep[[1]]

pdf(file = outpdf, useDingbats = FALSE)

m.clean <- ggplot(params.dat.merge, aes(x = estimate1, y = estimate2)) + 
  geom_point(alpha = 0.25, color = "grey75")  + 
  geom_density_2d() + 
  facet_wrap(~param) + 
  ggtitle("Clean", paste("Nbins:", length(bins.keep))) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.clean)

# Scramble param1  --------------------------------------------------------

set.seed(2)
params.dat1.scrambled <- params.dat1 %>%
  group_by(bin) %>%
  mutate(param = sample(x = param, size = length(param), replace = FALSE)) %>%
  filter(param %in% params.keep) %>%
  filter(abs(estimate1) < 5) 
  
m.scrambled <- ggplot(left_join(params.dat1.scrambled, params.dat2), aes(x = estimate1, y = estimate2)) + 
  geom_point(alpha = 0.25, color = "grey75")  + 
  geom_density_2d() + 
  facet_wrap(~param) + 
  ggtitle("Scrambled", paste("Nbins:", length(bins.keep))) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.scrambled)

multiplot(m.clean, m.scrambled)


# Show anticorrelation on UMAP  -------------------------------------------





# HSPCs -------------------------------------------------------------------


bset.dat <- params.dat2 %>%
  group_by(bin) %>%
  filter(min(estimate2) > 0)
bset.hspcs <- unique(bset.dat$bin)

jout.hspcs <- GetCutFracsFromBinSet(bset.hspcs)

ggplot(jout.hspcs$dat.annot1.merge, aes(x = umap1, y = umap2, color = log2(Winsorize(cutfrac, probs = c(0.02, 0.98))))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jmark1, paste("HSPC-specific K9me3 signal", length(bset.hspcs), "bins")) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(jout.hspcs$dat.annot2.merge, aes(x = umap1, y = umap2, color = log2(Winsorize(cutfrac, probs = c(0.02, 0.98))))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jmark2, paste("HSPC-specific K9me3 signal", length(bset.hspcs), "bins")) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


# Eryths ------------------------------------------------------------------

bset.dat <- params.dat2 %>%
  group_by(bin) %>%
  filter(estimate2[[1]] > 0 & median(estimate2) < 0)
bset.eryth <- unique(bset.dat$bin)

jout.eryth <- GetCutFracsFromBinSet(bset.eryth)

ggplot(jout.eryth$dat.annot1.merge, aes(x = umap1, y = umap2, color = log2(Winsorize(cutfrac, probs = c(0.02, 0.98))))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jmark1, paste("Eryth-specific K9me3 signal", length(bset.eryth), "bins")) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(jout.eryth$dat.annot2.merge, aes(x = umap1, y = umap2, color = log2(Winsorize(cutfrac, probs = c(0.02, 0.98))))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jmark2, paste("Eryth-specific K9me3 signal", length(bset.eryth), "bins")) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

 # Granus ------------------------------------------------------------------ 

print(head(params.dat2))

bset.dat <- params.dat2 %>%
  group_by(bin) %>%
  filter(estimate2[[2]] > 0 & median(estimate2) < 0)

bset.granu <- unique(bset.dat$bin)

jout.granu <- GetCutFracsFromBinSet(bset.granu)

ggplot(jout.granu$dat.annot1.merge, aes(x = umap1, y = umap2, color = log2(Winsorize(cutfrac, probs = c(0.02, 0.98))))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jmark1, paste("Granu-specific K9me3 signal", length(bset.granu), "bins")) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(jout.granu$dat.annot2.merge, aes(x = umap1, y = umap2, color = log2(Winsorize(cutfrac, probs = c(0.02, 0.98))))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jmark2, paste("Granu-specific K9me3 signal", length(bset.granu), "bins")) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")



# Bcells ------------------------------------------------------------------

print(head(params.dat2))

bset.dat <- params.dat2 %>%
  group_by(bin) %>%
  filter(estimate2[[3]] > 0 & median(estimate2) < 0)

bset.bcell <- unique(bset.dat$bin)

jout.bcell <- GetCutFracsFromBinSet(bset.bcell)

ggplot(jout.bcell$dat.annot1.merge, aes(x = umap1, y = umap2, color = log2(Winsorize(cutfrac, probs = c(0.02, 0.98))))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jmark1, paste("Bcell-specific K9me3 signal", length(bset.bcell), "bins")) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(jout.bcell$dat.annot2.merge, aes(x = umap1, y = umap2, color = log2(Winsorize(cutfrac, probs = c(0.02, 0.98))))) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  ggtitle(jmark2, paste("Bcell-specific K9me3 signal", length(bset.bcell), "bins")) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

dev.off()

