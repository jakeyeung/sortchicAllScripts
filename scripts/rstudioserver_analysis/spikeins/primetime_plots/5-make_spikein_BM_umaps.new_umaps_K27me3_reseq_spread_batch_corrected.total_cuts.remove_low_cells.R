# Jake Yeung
# Date of Creation: 2021-02-02
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/5-make_spikein_BM_umaps.new_umaps_K27me3_reseq_spread_batch_corrected.total_cuts.R
# 

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

# zval <- 2  # 95% con interval
zval <- 2.58  # 99% con interval

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/global_hist_mod_BM"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/batch_corrected"
outpdf <- file.path(outdir, paste0("bone_marrow_hist_mod_differences.", Sys.Date(), ".normtorep3.same_annot_file.K27me3reseq.spread.batch_corrected.totalcuts.mincuts_1000.pdf"))
outrds <- file.path(outdir, paste0("bone_marrow_hist_mod_differences.", Sys.Date(), ".normtorep3.same_annot_file.K27me3reseq.spread.batch_corrected.totalcuts.mincuts_1000.rds"))

make.plots <- TRUE

if (make.plots){
  pdf(outpdf, useDingbats = FALSE)
} else {
  print("Skip making plots")
}

# Load UMAPs --------------------------------------------------------------

indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28")
dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.meta <- file.path(indir.meta, paste0("metadata_batch_corrected.", jmark, ".2020-12-28.txt"))
  print(inf.meta)
  dat.meta <- fread(inf.meta)
  # adjuast plate 
  dat.meta$plate <- sapply(dat.meta$cell, function(x) ClipLast(x, jsep = "_"))
  return(dat.meta)
})


# Plot plate effects ------------------------------------------------------

# jmark <- "H3K27me3"
for (jmark in jmarks){
  
  m1 <- ggplot(dat.metas[[jmark]] %>% filter(!is.na(spikein_cuts)), aes(x = umap1, y = umap2, color = cuts_total / spikein_cuts)) + 
    geom_point() + 
    theme_bw() +
    scale_color_viridis_c() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  m2 <- ggplot(dat.metas[[jmark]] %>% filter(!is.na(spikein_cuts)), aes(x = cluster, y = log2(cuts_total / spikein_cuts))) + 
    geom_jitter(width = 0.1) + 
    theme_bw() +
    facet_wrap(~plate) + 
    ggtitle(jmark) + 
    geom_boxplot(outlier.shape = NA) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  m3 <- ggplot(dat.metas[[jmark]] %>% filter(!is.na(spikein_cuts)), aes(x = cluster, y = log2(cuts_total))) + 
    geom_jitter(width = 0.1) + 
    theme_bw() +
    facet_wrap(~plate) + 
    ggtitle(jmark) + 
    geom_boxplot(outlier.shape = NA) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  m4 <- ggplot(dat.metas[[jmark]] %>% filter(!is.na(spikein_cuts)), aes(x = cluster, y = cuts_total / spikein_cuts)) + 
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
    filter(!is.na(spikein_cuts) & cuts_total > 1000)
    # filter(!is.na(spikein_cuts) & jrep == "rep3")
  
  # make relative to granulocytes
  jchange <- "HSPCs"
  # jsub$clst <- gsub("Granulocytes", "aGranulocytes", jsub$cluster)
  jsub$clst <- gsub(jchange, paste0("a", jchange), jsub$cluster)
  jsub$plate <- gsub("PZ-BM-rep3", "PZ-BM-repa3", jsub$plate)
  
  jfit <- lm(formula = log2(cuts_total / spikein_cuts) ~ 1 + plate + clst, data = jsub)
  jfit.null <- lm(formula = log2(cuts_total / spikein_cuts) ~ 1 + plate, data = jsub)
  # jfit <- lm(formula = log2(cuts_in_peak / spikein_cuts) ~ 1 + plate + clst, data = jsub)
  # jfit.null <- lm(formula = log2(cuts_in_peak / spikein_cuts) ~ 1 + plate, data = jsub)
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

jmark <- "H3K27me3"
jfits.ds.lst.bymark[[jmark]]$jfits.ds.lst$jfit
jfits.ds.lst.bymark$H3K4me1$jfits.ds.lst$jfit

infcheck <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq/varfilt_2/BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.metadata.txt")
dat.check <- fread(infcheck)

infcheck2 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/metadata_batch_corrected.arranged_by_lineage.H3K27me3.2021-01-02.txt")
dat.check2 <- fread(infcheck2)

infcheck3 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned/cell_cluster_table_with_spikeins.H3K27me3.2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
dat.check3 <- fread(infcheck3)

infcheck4 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq/varfilt_2/BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.metadata.txt")
dat.check4 <- fread(infcheck4)

range(dat.check3$cuts_total)


# check total counts

indir.cuts.other <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/count_tables_from_bins/binsize_10000")
indir.cuts.k27me3 <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/bams_split_by_cluster/bams_merged_by_cluster/counts_tables_10000")

bin.mat.lst <- lapply("H3K27me3", function(jmark){
  print(jmark)
  indir.cuts <- ifelse(jmark == "H3K27me3", indir.cuts.k27me3, indir.cuts.other)
  infs.cuts <- list.files(path = indir.cuts, pattern = paste0(jmark, ".*.txt"), all.files = TRUE, full.names = TRUE)
  
  dat.csv.lst <- lapply(infs.cuts, function(inf.tmp){
    print(inf.tmp)
    # mat.tmp <- ReadMatTSSFormat(inf.tmp, as.sparse = TRUE, add.coord = TRUE, sort.rnames = FALSE)
    mat.tmp <- ReadMatSlideWinFormat(inf.tmp, as.sparse = TRUE, sort.rnames = FALSE, add.chromo = TRUE)
    rownames(mat.tmp) <- paste("chr", rownames(mat.tmp), sep = "")
    return(mat.tmp)
  })
  rnames.all <- sort(unique(unlist(lapply(dat.csv.lst, function(jmat) rownames(jmat)))))
  mat.tmp <- cbind.fill.lst(dat.csv.lst, all.rnames = rnames.all)
  return(mat.tmp)
})

# 
# dat.totalcuts <- data.frame(cell = colnames(bin.mat.lst[[jmark]]), totalcuts = colSums(bin.mat.lst[[jmark]]), stringsAsFactors = FALSE) %>%
#   left_join(., dat.metas[[jmark]])
# 
# cells.keep <- subset(dat.check4, is.good)$samp
# 
# dat.check4.filt <- subset(dat.check4, is.good)
# 
# dat.check4.merge <- left_join(dat.check4.filt, subset(dat.totalcuts, select = c(cell, totalcuts)))
# 
# ggplot(dat.check4.merge, aes(x = chromocounts, totalcuts)) + 
#   geom_point() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# fwrite(subset(dat.totalcuts, cell %in% cells.keep) %>% dplyr::select(-totalcuts), file = "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pseudobulk_tables_TSS_and_bins/single_cell_total_cuts_spikeins_by_plate_H3K27me3.txt", sep = "\t")
# 
# ggplot(dat.totalcuts, aes(x = ctype, y = log2(cuts_total / spikein_cuts))) + 
#   geom_point() + 
#   facet_wrap(~plate) + 
#   geom_boxplot() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# # sum across celltypes
# 
# dat.totalcuts.sum <- dat.totalcuts %>%
#   filter(ctype %in% c("HSPCs", "Eryths")) %>%
#   group_by(plate, ctype) %>%
#   summarise(totalcuts = sum(totalcuts),
#             cuts_total = sum(cuts_total),
#             spikein_cuts = sum(spikein_cuts)) %>%
#   rowwise() %>%
#   mutate(ratio1 = log2(totalcuts / spikein_cuts),
#          ratio2 = log2(cuts_total / spikein_cuts))
# print(dat.totalcuts.sum)
# 
# dat.totalcuts.sum2 <- dat.totalcuts.sum %>%
#   group_by(plate) %>%
#   summarise(ratio1.diff = ratio1[2] - ratio1[1],
#             ratio2.diff = ratio2[2] - ratio2[1]) %>%
#   filter(!is.na(ratio1.diff))
# 
# 
# jjchange <- "HSPCs"
# jjplate <- unique(dat.totalcuts$plate[2:3])
# jjsub <- subset(dat.totalcuts, plate %in% jjplate)
# jjsub$clst <- gsub(jjchange, paste0("a", jjchange), jjsub$cluster)
# jjsub$plate <- gsub("PZ-BM-rep3", "PZ-BM-repa3", jjsub$plate)
# # jjfit <- lm(formula = log2(cuts_total / spikein_cuts) ~ 1 + plate + clst, data = jjsub)
# jjfit <- lm(formula = log2(cuts_total / spikein_cuts) ~ 1  + clst, data = jjsub)
# 
# ggplot(jjsub, aes(x = ctype, y = log2(totalcuts / spikein_cuts))) + 
#   geom_point() + 
#   facet_wrap(~plate) + 
#   geom_boxplot() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# jjsub.sum <- jjsub %>%
#   filter(ctype %in% c("HSPCs", "Eryths")) %>%
#   group_by(ctype, plate) %>%
#   summarise(log2ratio = log2(sum(cuts_total) / sum(spikein_cuts)))
# 
# jjsub.sum <- jjsub %>%
#   filter(ctype %in% c("HSPCs", "Eryths")) %>%
#   group_by(ctype, plate) %>%
#   summarise(log2ratio = mean(log2(cuts_total / spikein_cuts)))
# 
# print(jjfit)
# 
# 
# 
# ggplot(dat.totalcuts.sum, aes(x = ctype, y = ratio1)) + 
#   geom_point() + 
#   geom_boxplot()  + 
#   facet_wrap(~plate) + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# 


# Show spikeins in UMAPs --------------------------------------------------

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

jprobs <- c(0.05, 0.95)

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
  dat.umap <- jfits.ds.lst.bymark[[jmark]]$jsub.effects %>%
    rowwise() %>%
    mutate(l2r = log2(cuts_total / spikein_cuts) - Estimate) %>%
    # ungroup() %>%
    group_by(cluster) %>%
    mutate(l2r.win = DescTools::Winsorize(l2r, probs = jprobs))
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = l2r.win)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})
print(m.lst.umap)

# set same UMAP
jmerged <- lapply(jmarks, function(jmark){
  jtmp <- jfits.ds.lst.bymark[[jmark]]$jsub.effects
  jtmp$mark <- jmark
  return(jtmp)
}) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(l2r = log2(cuts_total / spikein_cuts) - Estimate) %>%
  group_by(mark, clst) %>%
  mutate(l2r.win = DescTools::Winsorize(l2r, probs = jprobs))
  # mutate(l2r.win = l2r)
  # group_by(mark) %>%
  # mutate(l2r.win = scale(l2r.win, center = TRUE, scale = FALSE))

# shift so HSPCs are centered at zero? 
jshifts <- jmerged %>%
  group_by(mark) %>%
  filter(clst == "aHSPCs") %>%
  # summarise(l2r.med = ifelse(mark == "H3K4me1", median(l2r.win), median(l2r.win)))
  summarise(l2r.med = median(l2r.win)) %>%
  mutate(l2r.med = ifelse(mark == "H3K4me1", -0.75, 0))


jmerged2 <- left_join(jmerged, jshifts)

m <- ggplot(jmerged2, aes(x = umap1, y = umap2, color = l2r.win - l2r.med)) + 
  geom_point(size = 0.7) + 
  theme_minimal(8) + 
  scale_color_viridis_c() + 
  ggtitle(jmark) + 
  facet_wrap(~mark, scales = "free") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m)

m.lst.boxplots.check <- ggplot(jmerged2, 
              aes(x = forcats::fct_reorder(.f = clst, .x = l2r.win, .desc = TRUE, .fun = median), y = l2r.win - l2r.med)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(alpha = 0.25, width = 0.1) + 
    theme_bw(24) + 
    facet_wrap(~mark) + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

print(m.lst.boxplots.check)

m.lst.boxplots <- lapply(jmarks, function(jmark){
  print(jmark)
  dat.umap <- jfits.ds.lst.bymark[[jmark]]$jsub.effects
  m <- ggplot(dat.umap, 
              aes(x = forcats::fct_reorder(.f = clst, .x = log2(cuts_total/spikein_cuts) - Estimate, .desc = TRUE, .fun = median), y = log2(cuts_total / spikein_cuts) - Estimate)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(alpha = 0.25, width = 0.1) + 
    theme_bw(24) + 
    ggtitle(jmark) + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  return(m)
})

print(m.lst.boxplots)

JFuncs::multiplot(m.lst.boxplots[[1]], m.lst.boxplots[[2]], m.lst.boxplots[[3]], m.lst.boxplots[[4]], cols = 4)

ymax <- 1
ymin <- -4.5
m.lst.effects <- lapply(jmarks, function(jmark){
  print(jmark)
  
  jsub.tmp <- jfits.ds.lst.bymark[[jmark]]$jfits.ds.lst$jfit.merge %>%
    filter(!grepl("platePZ", param))
  
  m <- ggplot(jsub.tmp, aes(x = forcats::fct_reorder(.f = param, .x = est, .desc = TRUE, .fun = median), y = est, ymin = est - zval * est.se, ymax = est + zval * est.se)) + 
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
ggplot(dat.umap %>% filter(cluster == "Eryths" & umap2 < -6), aes(x = -umap1, y = log2(cuts_total/spikein_cuts) - Estimate, color = batch)) + 
  geom_point(size = 2) + 
  ggtitle("Eryths, umap2 < -6, umap1 as pseudotime") + 
  stat_smooth(geom='smooth', alpha=0.25, se=FALSE, color = 'blue', size = 2) + 
  # geom_smooth(method = "loess", se = FALSE, color = 'blue', size = 2, alpha = 0.5) + 
  theme_bw(24) +  
  scale_color_manual(values = jcols) + 
  xlab("Pseudotime") + 
  ylab("log2(cuts/spikeins)") + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap %>% filter(cluster == "Eryths" & umap2 < -6), aes(x = -umap1, y = log2(cuts_total/spikein_cuts) - Estimate)) + 
  geom_point(size = 2, color = "#0072B2") + 
  stat_smooth(geom='line', alpha=0.85, se=FALSE, color = 'grey', size = 2) + 
  ggtitle("Eryths, umap2 < -6, umap1 as pseudotime") + 
  # geom_smooth(method = "loess", se = FALSE, color = 'blue', size = 2, alpha = 0.5) + 
  theme_bw(24) +  
  scale_color_manual(values = jcols) + 
  xlab("Pseudotime") + 
  ylab("log2(cuts/spikeins)") + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap %>% filter(cluster == "HSPCs"), aes(x = -umap1, y = log2(cuts_total/spikein_cuts) - Estimate)) + 
  geom_point() + 
  theme_bw(24) +  
  xlab("Pseudotime") + 
  ylab("Corrected log2(cuts/spikeins)") + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())






# Write batch correction to output ----------------------------------------

# saveRDS(jfits.ds.lst.bymark, file = outrds)



# Do pseudotime  ----------------------------------------------------------

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
  geom_line(size = 2, alpha = 0.5) + 
  #  geom_point(data = jfilt.full, mapping = aes(color = log2(cuts_total / spikein_cuts))) + 
  geom_point(data = jfilt.full.join) + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(pcurve.dat, aes(x = umap1, y = umap2)) + 
  # geom_line(mapping = aes(color = ptime)) + 
  geom_line(size = 2, alpha = 0.5) + 
  geom_point(data = jfilt.full, mapping = aes(color = log2(cuts_total / spikein_cuts))) + 
  # geom_point(data = jfilt.full.join) + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(pcurve.dat, aes(x = umap1, y = umap2)) + 
  # geom_line(mapping = aes(color = ptime)) + 
  geom_line(size = 2, alpha = 0.5) + 
  geom_point(data = jfilt.full, mapping = aes(color = batch)) + 
  # geom_point(data = jfilt.full.join) + 
  theme_bw() + 
  scale_color_manual(values = c("grey", "red", "blue")) + 
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
  #  geom_point(data = jfilt.full, mapping = aes(color = log2(cuts_total / spikein_cuts))) + 
  geom_point(data = jfilt.full.join) + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggplot(jfilt.full.join, aes(x = log2(cuts_total / spikein_cuts), y = ptime)) + 
#   geom_point() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(jfilt.full.join, aes(x = log2(cuts_total / spikein_cuts), y = ptime)) + 
#   geom_point() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())





if (make.plots){
  dev.off()
} else {
  print("Skipping making plots")
}

