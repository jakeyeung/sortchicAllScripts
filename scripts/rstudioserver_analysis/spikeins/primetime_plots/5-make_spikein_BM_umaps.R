# Jake Yeung
# Date of Creation: 2020-11-29
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/5-make_spikein_BM_umaps.R
# Fig 4 init


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

DoFitsDownstream <- function(jfit, jfit.null, jgrep = "^experi|Intercept"){
  jcompare <- anova(jfit.null, jfit)
  jfit.effects <- data.frame(param = rownames(summary(jfit)$coefficients), summary(jfit)$coefficients, stringsAsFactors = FALSE) %>%
    filter(!grepl(jgrep, param))
  jfit.int <- data.frame(param = rownames(summary(jfit)$coefficients), summary(jfit)$coefficients, stringsAsFactors = FALSE) %>%
    filter(param == "(Intercept)") %>%
    rowwise()
  jint <- jfit.int$Estimate
  jint.se <- jfit.int$Std..Error
  
  jfit.int <- data.frame(param = rownames(summary(jfit)$coefficients), summary(jfit)$coefficients, stringsAsFactors = FALSE) %>%
    filter(param == "(Intercept)") %>%
    rowwise() %>%
    mutate(est = jint,
           est.se = jint.se)
  
  jfit.dat <- jfit.effects %>%
    rowwise() %>%
    mutate(est = Estimate + jint,
           est.se = sqrt(Std..Error ^ 2 + jint.se ^ 2))
  
  jfit.merge <- rbind(jfit.dat, jfit.int)
  return(list(jfit = jfit, jfit.null = jfit.null, jcompare = jcompare, jfit.merge = jfit.merge))
}


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/global_hist_mod_BM"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2"
outpdf <- file.path(outdir, paste0("bone_marrow_hist_mod_differences.", Sys.Date(), ".normtorep3.same_annot_file.pdf"))

make.plots <- TRUE

if (make.plots){
  pdf(outpdf, useDingbats = FALSE)
} else {
  print("Skip making plots")
}

# Load UMAPs --------------------------------------------------------------


# indir.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap"
# indir.meta.k9 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins"
indir.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins"
jdate <- "2020-11-18"
dat.metas.init <- lapply(jmarks, function(jmark){
  fname <- paste0("cell_cluster_table_with_spikeins.", jmark, ".", jdate, ".dupfilt.txt")
  inf.meta <- file.path(indir.meta, fname)
  # if (jmark == "H3K9me3"){
  # } else {
  #   fname <- paste0("geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.", jmark, ".txt")
  #   inf.meta <- file.path(indir.meta, fname)
  # }
  dat.tmp <- fread(inf.meta)
  dat.tmp <- subset(dat.tmp, select = -c(umap1, umap2))
})

# update UMAPs
indir.umaps.k9 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2"
indir.umaps <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap"
dat.umaps <- lapply(jmarks, function(jmark){
  if (jmark != "H3K9me3"){
    fname.tmp <- paste0("geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.", jmark, ".txt")
    inf.umaps <- file.path(indir.umaps, fname.tmp)
  } else {
    fname.tmp <- paste0("BM_celltypes.bincutoff_0.binskeep_1000.byplate.szname_none.niter_500.reorder_rownames.dupfilt.2020-11-23.", jmark, ".txt")
    inf.umaps <- file.path(indir.umaps.k9, fname.tmp)
  }
  #   inf.meta <- file.path(indir.meta, fname)
  dat.tmp <- fread(inf.umaps)
  dat.tmp$mark <- jmark
  dat.tmp <- subset(dat.tmp, select = c(cell, umap1, umap2))
  return(dat.tmp)
})

dat.metas <- lapply(jmarks, function(jmark){
  left_join(dat.metas.init[[jmark]], dat.umaps[[jmark]], by = "cell")
})

# check basophils
basos <- subset(dat.metas.init$H3K4me1, cluster == "Basophils")$cell
dat.sub <- subset(dat.umaps$H3K4me1, cell %in% basos)

# check glmpca output
inf.check <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins.from_sitecount_mat/lda_outputs.count_mat_from_sitecount_mat.H3K4me1.filtNAcells_allbins.K-30.binarize.FALSE/ldaOut.count_mat_from_sitecount_mat.H3K4me1.filtNAcells_allbins.K-30.Robj")
load(inf.check, v=T)

colnames(count.mat)[(which(colnames(count.mat) %in% basos))]


# Plot plate effects ------------------------------------------------------

# jmark <- "H3K27me3"
for (jmark in jmarks){
  
  m1 <- ggplot(dat.metas[[jmark]] %>% filter(!is.na(spikein_cuts)), aes(x = umap1, y = umap2, color = cuts_in_peak / spikein_cuts)) + 
    geom_point() + 
    theme_bw() +
    scale_color_viridis_c() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  m2 <- ggplot(dat.metas[[jmark]] %>% filter(!is.na(spikein_cuts)), aes(x = cluster, y = log2(cuts_in_peak / spikein_cuts))) + 
    geom_jitter(width = 0.1) + 
    theme_bw() +
    facet_wrap(~plate) + 
    ggtitle(jmark) + 
    geom_boxplot(outlier.shape = NA) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  m3 <- ggplot(dat.metas[[jmark]] %>% filter(!is.na(spikein_cuts)), aes(x = cluster, y = log2(cuts_in_peak))) + 
    geom_jitter(width = 0.1) + 
    theme_bw() +
    facet_wrap(~plate) + 
    ggtitle(jmark) + 
    geom_boxplot(outlier.shape = NA) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  m4 <- ggplot(dat.metas[[jmark]] %>% filter(!is.na(spikein_cuts)), aes(x = cluster, y = cuts_in_peak / spikein_cuts)) + 
    geom_jitter(width = 0.1) + 
    theme_bw() +
    facet_wrap(~plate) + 
    ggtitle(jmark) + 
    geom_boxplot(outlier.shape = NA) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m1)
  print(m2)
  print(m3)
  print(m4)
}


# Calculate plate effects -------------------------------------------------

# fit cell and plate effects

jfits.ds.lst.bymark <- lapply(jmarks, function(jmark){
  print(jmark)
  jsub <- dat.metas[[jmark]] %>%
    filter(!is.na(spikein_cuts))
  
  # make relative to granulocytes
  jchange <- "HSPCs"
  # jsub$clst <- gsub("Granulocytes", "aGranulocytes", jsub$cluster)
  jsub$clst <- gsub(jchange, paste0("a", jchange), jsub$cluster)
  jsub$plate <- gsub("PZ-BM-rep3", "PZ-BM-repa3", jsub$plate)
  
  jfit <- lm(formula = log2(cuts_in_peak / spikein_cuts) ~ 1 + plate + clst, data = jsub)
  jfit.null <- lm(formula = log2(cuts_in_peak / spikein_cuts) ~ 1 + plate, data = jsub)
  jsum <- anova(jfit.null, jfit)
  pval <- jsum$`Pr(>F)`[[2]]
  jfit.ds.lst <- DoFitsDownstream(jfit, jfit.null)
  jfit.ds.lst$jfit.merge$mark <- jmark
  jfit.ds.lst$jfit.merge$param <- gsub(pattern = "clst", replacement = "", x = jfit.ds.lst$jfit.merge$param)
  # jfit.ds.lst$jfit.merge$param <- gsub(pattern = "\\(Intercept\\)", replacement = "Granulocytes", x = jfit.ds.lst$jfit.merge$param)
  jfit.ds.lst$jfit.merge$param <- gsub(pattern = "\\(Intercept\\)", replacement = jchange, x = jfit.ds.lst$jfit.merge$param)
  
  jbaseline <- subset(jfit.ds.lst$jfit.merge, param == jchange)$Estimate
  jfit.ds.lst$jfit.merge$est <- jfit.ds.lst$jfit.merge$est - jbaseline
  
  # subtract plate effects from jsub
  
  plate.effects <- jfit.ds.lst$jfit.merge %>%
    dplyr::filter(grepl("^platePZ", param)) %>%
    dplyr::mutate(plate = gsub("plate", "", param))
  
  # if NA then its intercept plate 
  jsub.effects <- left_join(jsub, plate.effects, by = "plate") %>%
    rowwise() %>%
    mutate(Estimate = ifelse(is.na(Estimate), 0, Estimate),
           Std..Error = ifelse(is.na(Std..Error), 0, Std..Error))
  
  return(list(jfits.ds.lst = jfit.ds.lst, jsub.effects = jsub.effects))
})


# Show spikeins in UMAPs --------------------------------------------------

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
m.lst.ctype <- lapply(jmarks, function(jmark){
  print(jmark)
  dat.umap <- jfits.ds.lst.bymark[[jmark]]$jsub.effects
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + 
    theme_bw() + 
    facet_wrap(~cluster) + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})
print(m.lst.ctype)

m.lst.umap <- lapply(jmarks, function(jmark){
  print(jmark)
  dat.umap <- jfits.ds.lst.bymark[[jmark]]$jsub.effects
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = log2(cuts_in_peak / spikein_cuts) - Estimate)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})
print(m.lst.umap)

m.lst.boxplots <- lapply(jmarks, function(jmark){
  print(jmark)
  dat.umap <- jfits.ds.lst.bymark[[jmark]]$jsub.effects
  m <- ggplot(dat.umap, aes(x = forcats::fct_reorder(.f = clst, .x = log2(cuts_in_peak/spikein_cuts) - Estimate, .desc = TRUE, .fun = median), y = log2(cuts_in_peak / spikein_cuts) - Estimate)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(alpha = 0.25, width = 0.1) + 
    theme_bw(24) + 
    ggtitle(jmark) + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  return(m)
})

print(m.lst.boxplots)

ymax <- 1
ymin <- -4.5
m.lst.effects <- lapply(jmarks, function(jmark){
  print(jmark)
  
  jsub.tmp <- jfits.ds.lst.bymark[[jmark]]$jfits.ds.lst$jfit.merge %>%
    filter(!grepl("platePZ", param))
  
  m <- ggplot(jsub.tmp, aes(x = forcats::fct_reorder(.f = param, .x = est, .desc = TRUE, .fun = median), y = est, ymin = est - 2 * est.se, ymax = est + 2 * est.se)) + 
    geom_errorbar() + 
    theme_bw(24) + 
    ggtitle(jmark) + 
    scale_color_viridis_c() + 
    ylim(c(ymin, ymax)) + 
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  return(m)
})

print(m.lst.effects)



# add pseudotime of eryths? 
jmark <- "H3K27me3"
dat.umap <- jfits.ds.lst.bymark[[jmark]]$jsub.effects

jcols <- c("grey", "red", "blue")
ggplot(dat.umap %>% filter(cluster == "Eryths" & umap2 < -6), aes(x = -umap1, y = log2(cuts_in_peak/spikein_cuts) - Estimate, color = stype)) + 
  geom_point(size = 2) + 
  ggtitle("Eryths, umap2 < -6, umap1 as pseudotime") + 
  stat_smooth(geom='smooth', alpha=0.25, se=FALSE, color = 'blue', size = 2) + 
  # geom_smooth(method = "loess", se = FALSE, color = 'blue', size = 2, alpha = 0.5) + 
  theme_bw(24) +  
  scale_color_manual(values = jcols) + 
  xlab("Pseudotime") + 
  ylab("log2(cuts/spikeins)") + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap %>% filter(cluster == "Eryths" & umap2 < -6), aes(x = -umap1, y = log2(cuts_in_peak/spikein_cuts) - Estimate)) + 
  geom_point(size = 2, color = "#0072B2") + 
  stat_smooth(geom='line', alpha=0.85, se=FALSE, color = 'grey', size = 2) + 
  ggtitle("Eryths, umap2 < -6, umap1 as pseudotime") + 
  # geom_smooth(method = "loess", se = FALSE, color = 'blue', size = 2, alpha = 0.5) + 
  theme_bw(24) +  
  scale_color_manual(values = jcols) + 
  xlab("Pseudotime") + 
  ylab("log2(cuts/spikeins)") + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap %>% filter(cluster == "HSPCs"), aes(x = -umap1, y = log2(cuts_in_peak/spikein_cuts) - Estimate)) + 
  geom_point() + 
  theme_bw(24) +  
  xlab("Pseudotime") + 
  ylab("Corrected log2(cuts/spikeins)") + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



if (make.plots){
  dev.off()
} else {
  print("Skipping making plots")
}


