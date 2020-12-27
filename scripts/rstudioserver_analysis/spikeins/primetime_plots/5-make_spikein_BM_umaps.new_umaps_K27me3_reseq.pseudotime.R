# Jake Yeung
# Date of Creation: 2020-12-26
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/5-make_spikein_BM_umaps.new_umaps_K27me3_reseq.pseudotime.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/global_hist_mod_BM"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2"
outpdf <- file.path(outdir, paste0("bone_marrow_hist_mod_differences.", Sys.Date(), ".normtorep3.same_annot_file.K27me3reseq.spread.pseudotime.pdf"))

make.plots <- TRUE

if (make.plots){
  pdf(outpdf, useDingbats = FALSE)
} else {
  print("Skip making plots")
}

# Load UMAPs --------------------------------------------------------------

indir.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned"
jdate <- "2020-12-23"
dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0("cell_cluster_table_with_spikeins.", jmark, ".", jdate, ".umap_spread.final.txt")
  inf.meta <- file.path(indir.meta, fname)
  dat.tmp <- fread(inf.meta) %>%
    rowwise() %>% 
    mutate(stype = batch,
           platenbr = plate,
           plate = ClipLast(cell, jsep = "_"))
})



# Show K27me3 and UMAP  ---------------------------------------------------

ggplot(dat.metas$H3K27me3, aes(x = umap1, y = umap2, color = log2(cuts_in_peak / spikein_cuts))) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

clstrs <- sort(unique(dat.metas$H3K27me3$cluster))
for (jclst in clstrs){
  m <- ggplot(dat.metas$H3K27me3 %>% filter(cluster == jclst), aes(x = umap1, y = umap2, color = log2(cuts_in_peak / spikein_cuts))) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c() + 
    ggtitle(jclst) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
}


ggplot(dat.metas$H3K27me3 %>% filter(cluster == "Eryths"), aes(x = umap1, y = umap2, color = log2(cuts_in_peak / spikein_cuts))) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Check from glmpcaout ----------------------------------------------------

inf.glmpca <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final_glmpca_tweak_umap.H3K27me3_rep2rep3reseq/glmpca_factors.cleaned.H3K27me3.2020-12-23.txt"

# Predict?

# dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
# dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)

glmpca.k27me3 <- fread(inf.glmpca) %>%
  dplyr::rename(cell = V1) %>%
  left_join(., dat.metas$H3K27me3)

ggplot(glmpca.k27me3 %>% filter(cluster == "Eryths"), aes(x = dim5, y = dim6, color = log2(cuts_in_peak / spikein_cuts))) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jfilt.full <- dat.metas$H3K27me3 %>% filter(cluster == "Eryths")
jfilt <- jfilt.full %>%
  dplyr::select(c(cell, umap1, umap2)) %>%
  as.data.frame() 
rownames(jfilt) <- jfilt$cell
jfilt$cell <- NULL

pcurveout <- princurve::principal_curve(x = as.matrix(jfilt))

pcurve.dat <- data.frame(pcurveout$s, cell = rownames(pcurveout$s), ptime = pcurveout$lambda, stringsAsFactors = FALSE)

jfilt.full.join <- left_join(jfilt.full, subset(pcurve.dat, select = -c(umap1, umap2)), by = "cell")

ggplot(pcurve.dat, aes(x = umap1, y = umap2, color = ptime)) + 
  # geom_line(mapping = aes(color = ptime)) + 
  geom_line() + 
  #  geom_point(data = jfilt.full, mapping = aes(color = log2(cuts_in_peak / spikein_cuts))) + 
  geom_point(data = jfilt.full.join) + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jfilt.full.join, aes(x = log2(cuts_in_peak / spikein_cuts), y = ptime)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

inf.fits <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/bone_marrow_hist_mod_differences.2020-12-26.normtorep3.same_annot_file.K27me3reseq.spread.rds"
jfits.ds.lst.bymark <- readRDS(inf.fits)

head(jfits.ds.lst.bymark$H3K27me3$jsub.effects)

ggplot(jfits.ds.lst.bymark$H3K27me3$jsub.effects, aes(x = umap1, y = umap2, color = log2(cuts_in_peak / spikein_cuts) - Estimate)) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jfilt.full.join2 <- left_join(jfilt.full.join, subset(jfits.ds.lst.bymark$H3K27me3$jsub.effects, select = c(cell, Estimate)))

ggplot(jfilt.full.join2, aes(x = log2(cuts_in_peak / spikein_cuts) - Estimate, y = ptime)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jfilt.full.join2, aes(x = log2(cuts_in_peak / spikein_cuts), y = ptime)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jfilt.full.join2, aes(y = log2(cuts_in_peak / spikein_cuts) - Estimate, x = ptime / max(ptime), color = batch)) + 
  geom_point() + 
  facet_wrap(~batch) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jfilt.full.join2, aes(y = log2(cuts_in_peak / spikein_cuts) - Estimate, x = ptime / max(ptime))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("H3K27me3 erythroblast pseudotime") + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Check DE genes ----------------------------------------------------------


inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.allmarks_spikeins/poisson_fit_bins.H3K27me3.2020-12-25.newannot2.witherrors.MoreBins.newestmeta.RData"
# inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again/poisson_fit_bins.H3K27me3.2020-12-12.newannot2.witherrors.MoreBins.RData"
load(inf.de, v=T)

jmark <- "H3K27me3"
jnames <- names(jfits.lst); names(jnames) <- jnames

pvals.dat <- lapply(jnames, function(jname){
  data.frame(bin = jname, pval = jfits.lst[[jname]]$pval, stringsAsFactors = FALSE)
}) %>%
  bind_rows()

params.dat.all <- lapply(jnames, function(jname){
  x <- jfits.lst[[jname]]
  xkeep <- grepl("^Cluster.*.Estimate$", x = names(x))
  jparams <- x[xkeep]
  data.frame(bin = jname, param = names(jparams), estimate = unlist(jparams), mark = jmark, stringsAsFactors = FALSE)
}) %>%
  bind_rows() 

# https://stats.stackexchange.com/questions/315311/how-to-find-p-value-using-estimate-and-standard-error
params.dat.all.withpval <- lapply(jnames, function(jname){
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




# make params more readable
params.dat.all$ctype <- params.dat.all$param
params.dat.all$ctype <- gsub("Cluster", "", params.dat.all$ctype)
params.dat.all$ctype <- gsub(".Estimate", "", params.dat.all$ctype)


params.dat.all.withpval$ctype <- params.dat.all.withpval$param
params.dat.all.withpval$ctype <- gsub("Cluster", "", params.dat.all.withpval$ctype)
params.dat.all.withpval$ctype <- gsub(".Estimate", "", params.dat.all.withpval$ctype)

ggplot(params.dat.all %>% filter(abs(estimate) < 5), aes(x = estimate, fill = ctype)) + 
  geom_density() + 
  facet_wrap(~ctype) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(params.dat.all.withpval %>% filter(abs(estimate) < 5), aes(x = estimate, y = -log10(pval.param))) + 
  geom_point() + 
  facet_wrap(~ctype) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Assign bins to gene  ----------------------------------------------------

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

jbins <- unique(params.dat.all$bin)
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

bins.annot <- AnnotateCoordsFromList.GeneWise(coords.vec = jbins, inf.tss = "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.100000.bed", txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

# chekc HSPC-specifc genes

famous.genes.lst <- 
  list("Eryths" = c("Tal1", "Aqp1", "Gata1", "Comt", "Sphk1", "Hbb-bs", "Hbb-bt", "Sox6"),
       "Bcells" = c("Pax5", "Ebf1", "Cd79a", "Cd79b", "Snn", "Blnk", "Cd72", "Blk", "Kzf3", "Cd19"), 
       "NKs" = c("Stat4", "Tcf7", "Tbx21", "Cd3d", "Gimap4"), 
       "Granulocytes" = c("Cxcr2", "Mmp9", "Cd177", "Ncf4", "S100a9", "S100a8", "Ncf1", "Cebpe"),
       "Basophils" = c("Il4", "Il6", "Il1r1", "Arap3", "Cpa3"), 
       "pDCs" = c("Irf8", "Selpg", "Itgb7", "Ccr9", "Unc93b1"), 
       "DCs" = c("Ccl9", "Apoe", "Nlrp3", "Csf1r"), 
       "HSPCs" = c("Hoxa9", "Hoxa7", "Hoxa3", "Meis1", "Runx2", "Kit", "Hlf", "Hoxa10", "Erg", "Cd34", "Hoxa6", "Gata3", "Hoxb4")) 

jgene <- "Gata1"
jsub <- subset(bins.annot$out2.df, gene == jgene) %>%
  arrange(abs(dist.to.tss))
jbin <- jsub$region_coord[[1]]

jfits.sub <- subset(params.dat.all, bin == jbin)


# Find bins that are lost in all celltyupes -------------------------------

jfits.sum <- params.dat.all %>%
  group_by(bin) %>%
  summarise(estimate.mean = mean(estimate))  %>%
  arrange(estimate.mean) %>%
  filter(abs(estimate.mean) < 5) %>%
  left_join(., bins.annot$out2.df %>% dplyr::select(region_coord, gene, dist.to.tss), by = c("bin" = "region_coord")) %>%
  left_join(., pvals.dat, by = "bin")

jfits.sum.filt <- subset(jfits.sum, pval < 10^-50)

ggplot(jfits.sum, aes(x = estimate.mean, y = -log10(pval))) + 
  geom_point(alpha = 0.25)  + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# ggplot(jfits.sum %>% filter(abs(estimate.mean) < 5), aes(x = estimate.mean)) + 
ggplot(jfits.sum, aes(x = estimate.mean)) + 
  geom_density()  + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# pvals 
print(head(data.frame(jfits.sum.filt), n = 50))


# Find eryth-specific loss  -----------------------------------------------

jfits.wide <- params.dat.all %>%
  data.table::dcast(., formula = bin ~ ctype, value.var = "estimate") 
jfits.wide2 <- left_join(jfits.wide, jfits.sum)

jfits.eryth <- jfits.wide2 %>%
  mutate(estimate.eryth = Eryths - estimate.mean)  %>%
  arrange(estimate.eryth) 

jfits.eryth.sub <- subset(jfits.eryth, pval < 10^-10)


# Write pseudotime to output ----------------------------------------------

  
# indir.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned"
# outf.meta <- file.path(indir.meta, paste0("cell_cluster_table_with_spikeins.", jmark, ".", Sys.Date(), ".umap_spread.final.with_pseudotime.txt"))
# fwrite(jfilt.full.join2, file = outf.meta, sep = "\t")




# Dev off  ----------------------------------------------------------------



if (make.plots){
  dev.off()
} else {
  print("Skipping making plots")
}


