# Jake Yeung
# Date of Creation: 2020-10-31
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged/1-LDA_downstream.WithRatios.again.K9me3_clean_with_celltypes.ForPeter.R
# 



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)
library(scchicFuncs)

library(topicmodels)

library(ggrepel)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)


SplitGetLast <- function(x, jsplit = "-"){
  # get last element after splitting
  xsplit <- strsplit(x, jsplit)[[1]]
  xnew <- xsplit[[length(xsplit)]]
  return(xnew)
}



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

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_mouse.BMround2_umaps_and_ratios"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios"
dir.create(outdir)
outprefix <- file.path(outdir, paste0("spikeins_mouse.BMround2_umaps_and_ratios.colfix.celltyping.", Sys.Date(), ".WithRelLevels"))
outpdf <- paste0(outprefix, ".pdf")
pdf(file = outpdf, useDingbats = FALSE)


# Constnats ---------------------------------------------------------------

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#41f0d9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#28f9ff", "#88497e", "#bcf5c3", "#86f115")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"



# Add log2ratio -----------------------------------------------------------

inf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all/spikein_info_BM_round2_all.txt"

dat.spikeins.mat <- fread(inf.spikeins) %>%
  rowwise() %>%
  mutate(plate = as.numeric(SplitGetLast(experi, jsplit = "-")),
         jrep = GetRepBM(experiname = experi), 
         stype = AnnotateSortFromLayoutBMall(plate = plate, rowcoord = rowcoord, colcoord = colcoord, jrep = jrep, jmark = mark))


dat.spikeins.mat.sub <- subset(dat.spikeins.mat, select = c(samp, chromocounts, spikeincounts)) %>%
  rowwise() %>%
  mutate(l2r = log2(chromocounts / spikeincounts))
# 
# dat.merge <- left_join(dat.merge, jsub, by = c("cell" = "samp"))



# Load objects, k9me3 will be replaced ------------------------------------



inf.objs <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_mouse.BMround2all.varfilt.blfix/LDA_downstream_objects.RData"
load(inf.objs, v=T)

dat.umaps.lst <- split(dat.merge, dat.merge$mark)

# Load K9me3 separately cleaned up  ---------------------------------------

jmarkrep <- "H3K9me3"; names(jmarkrep) <- "H3K9me3"
inf.k9me3 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_H3K9me3_badclustremoved/lda_outputs.countmat_H3K9me3_newonly_badclustremoved.K-30.binarize.FALSE/ldaOut.countmat_H3K9me3_newonly_badclustremoved.K-30.Robj"
load(inf.k9me3, v=T)

tm.result.k9me3 <- posterior(out.lda)
dat.umap.k9me3 <- DoUmapAndLouvain(tm.result.k9me3$topics, jsettings) %>%
  left_join(., dat.spikeins.mat.sub, by = c("cell" = "samp"))
dat.umap.k9me3$mark <- jmarkrep


ggplot(dat.umap.k9me3, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# add var
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.impute.log <- log2(t(tm.result.k9me3$topics %*% tm.result.k9me3$terms))
dat.var.k9me3 <- CalculateVarAll(dat.impute.log, jchromos)

dat.umap.k9me3.merge <- left_join(dat.umap.k9me3, dat.var.k9me3, by = "cell")

tm.result.lst[[jmarkrep]] <- tm.result.k9me3
dat.umaps.lst[[jmarkrep]] <- dat.umap.k9me3.merge


jmarksall <- c(jmarks, jmarkrep)

for (jmark in jmarksall){
  m <- ggplot(dat.umaps.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jmark)
  print(m)
}


dat.umaps.long <- bind_rows(dat.umaps.lst) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"), 
         plate = as.numeric(SplitGetLast(experi, jsplit = "-")),
         rowcoord = AddPlateCoordinates(cell)$rowcoord,
         colcoord = AddPlateCoordinates(cell)$colcoord,
         jrep = GetRepBM(experiname = experi), 
         stype = AnnotateSortFromLayoutBMall(plate = plate, rowcoord = rowcoord, colcoord = colcoord, jrep = jrep, jmark = mark))

dat.merge <- dat.umaps.long

dat.meta <- subset(dat.umaps.long, select = c(cell, experi, plate, rowcoord, colcoord, stype))

# # Variance stuff ----------------------------------------------------------
# 
# jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
# dat.var.lst <- lapply(jmarks, function(jmark){
#   tm.result <- tm.result.lst[[jmark]]
#   dat.impute <- t(log2(tm.result$topics %*% tm.result$terms))
#   dat.var <- CalculateVarAll(dat.impute, jchromos)
# }) 
# 
# dat.var <- dat.var.lst %>%
#   bind_rows()
# 
# dat.merge <- left_join(dat.umaps.long, dat.var)
# 


# Plot umap with log2ratio  -----------------------------------------------

ggplot(dat.merge, aes(x = umap1, y = umap2, color = l2r)) + 
  geom_point() + 
  facet_wrap(~mark) + 
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# show boxplots
ggplot(dat.merge, aes(x = stype, y = l2r)) + 
  geom_boxplot() + 
  theme_bw() + 
  facet_wrap(~mark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# do K9me3 only 

louv2clst <- list("1" = "zHSPCs", 
                  "2" = "aGranulocytes", 
                  "3" = "aGranulocytes",
                  "4" = "Eryth",
                  "5" = "Lymphoid")

jhash <- hash::hash(louv2clst)

dat.k9me3 <- subset(dat.merge, mark == "H3K9me3") %>%
  rowwise() %>%
  mutate(cluster = AssignHash(x = as.character(louvain), louv2clst, null.fill = NA))


ggplot(dat.k9me3, aes(x = cluster, y = l2r)) + 
  geom_boxplot() + 
  theme_bw() + 
  facet_wrap(~experi) + 
  ylab("log2(chromocounts/spikeins)") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Get common celltypes for all marks ---------------------------------------------

ggplot(dat.merge %>% filter(mark == "H3K4me1"), aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# jcheck <- c("13", "12", "7")
# jcheck <- c("13")
jcheck <- c("3")
ggplot(dat.merge %>% filter(mark == "H3K4me1") %>% mutate(louvain = louvain %in% jcheck), aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + ggtitle(jcheck) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

louv2clst.k4me1 <- list("14" = "NKs", 
                        "10" = "pDCs", 
                        "2" = "Eryths", 
                        "6" = "Eryths", 
                        "8" = "Monocytes",
                        "9" = "BcellsNaive",
                        "4" = "BcellsPlasma",
                        "5" = "BcellsPlasma",
                        "11" = "Basophils",
                        "3" = "zHSPCs",
                        "13" = "aGranulocytes",
                        "12" = "aGranulocytes",
                        "7" = "aGranulocytes",
                        "1" = "Eosinophils")

jhash.k4me1 <- hash::hash(louv2clst.k4me1)

dat.k4me1 <- subset(dat.merge, mark == "H3K4me1")  %>%
  rowwise() %>%
  mutate(cluster = AssignHash(x = as.character(louvain), louv2clst.k4me1, null.fill = NA),
         cluster = ifelse(umap2 < -6.5 & umap2 >= -10, NA, cluster))

ggplot(dat.k4me1, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + ggtitle("H3K4me1") + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# H3K4me3 -----------------------------------------------------------------

ggplot(dat.merge %>% filter(mark == "H3K4me3"), aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# jcheck <- c("13", "12", "7")
# jcheck <- c("13")
# jcheck <- c("8")  # 15
# jcheck <- c("1")  # 15
# jcheck <- c("2", "5")  # 15
# jcheck <- c("7", "11", "3", "6")  # 15
jcheck <- "15"
ggplot(dat.merge %>% filter(mark == "H3K4me3") %>% mutate(louvain = louvain == jcheck), aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + ggtitle(jcheck) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

louv2clst.k4me3 <- list("14" = "NKs", 
                        "13" = "pDCs",
                        "9" = "DCs",
                        "4" = "DCs",
                        "8" = "BcellsNaive",
                        "15" = "BcellsPlasma",
                        "12" = "Eryths",
                        "1" = "Eryths",
                        "2" = "zHSPCs",
                        "5" = "zHSPCs",
                        "7" = "aGranulocytes",
                        "11" = "aGranulocytes",
                        "3" = "aGranulocytes",
                        "6" = "aGranulocytes")

jhash.k4me3 <- hash::hash(louv2clst.k4me3)

dat.k4me3 <- subset(dat.merge, mark == "H3K4me3")  %>%
  rowwise() %>%
  mutate(cluster = AssignHash(x = as.character(louvain), louv2clst.k4me3, null.fill = NA))

ggplot(dat.k4me3, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + ggtitle("H3K4me3") + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Fit K27me3 --------------------------------------------------------------


ggplot(dat.merge %>% filter(mark == "H3K27me3"), aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jcheck <- c("10", "2")
jcheck <- c("9")
# jcheck <- c("12")
jcheck <- c("6")
jcheck <- c("3")
jcheck <- c("7")
jcheck <- c("8")
jcheck <- c("5")
jcheck <- c("4")
jcheck <- c("11", "12", "1")

ggplot(dat.merge %>% filter(mark == "H3K27me3") %>% mutate(louvain = louvain %in% jcheck), aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + ggtitle(jcheck) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

louv2clst.k27me3 <- list("10" = "Eryths", 
                        "2" = "Eryths",
                        "9" = "NKs",
                        "6" = "pDCs",
                        "3" = "DCs",
                        "13" = "Basophils",
                        "8" = "BcellsNaive",
                        "7" = "BcellsPlasma",
                        "5" = "BcellsPlasma",
                        "4" = "zHSPCs",
                        "11" = "aGranulocytes",
                        "12" = "aGranulocytes",
                        "1" = "aGranulocytes"
                        )

jhash.k27me3 <- hash::hash(louv2clst.k27me3)

dat.k27me3 <- subset(dat.merge, mark == "H3K27me3")  %>%
  rowwise() %>%
  mutate(cluster = AssignHash(x = as.character(louvain), louv2clst.k27me3, null.fill = NA))

ggplot(dat.k27me3, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + ggtitle("H3K27me3") + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Plot celltypes ----------------------------------------------------------



dat.merge.annot <- list("H3K4me1" = dat.k4me1, "H3K4me3" = dat.k4me3, "H3K27me3" = dat.k27me3, "H3K9me3" = dat.k9me3)

for (jmark in jmarksall){
  m <- ggplot(dat.merge.annot[[jmark]], aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + 
    ggtitle(jmark) + 
    scale_color_manual(values = cbPalette, na.value = "grey95") + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}

for (jmark in jmarksall){
  m <- ggplot(dat.merge.annot[[jmark]], aes(x = umap1, y = umap2, color = stype)) + 
    geom_point() + 
    ggtitle(jmark) + 
    scale_color_manual(values = cbPalette, na.value = "grey95") + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}

for (jmark in jmarksall){
  m <- ggplot(dat.merge.annot[[jmark]], aes(x = umap1, y = umap2, color = log2(chromocounts / spikeincounts))) + 
    geom_point() + 
    ggtitle(jmark) + 
    scale_color_viridis_c() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}



# Fit ratio across celltypes  ---------------------------------------------

for (jmark in jmarksall){
  m <- ggplot(dat.merge.annot[[jmark]] %>% filter(!is.na(cluster)), aes(x = cluster, y = log2(chromocounts / spikeincounts))) + 
    geom_point() + 
    geom_boxplot() + 
    theme_bw() + 
    facet_wrap(~experi) + 
    ggtitle(jmark) + 
    xlab("") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  strip.text.x = element_text(size = 5))
  print(m)
}

for (jmark in jmarksall){
  m <- ggplot(dat.merge.annot[[jmark]] %>% filter(!is.na(cluster)) %>% mutate(stype = gsub(pattern = "LSK", replacement = "zLSK", x = stype)), aes(x = stype, y = log2(chromocounts / spikeincounts))) + 
    geom_point() + 
    geom_boxplot() + 
    theme_bw() + 
    facet_wrap(~experi) + 
    ggtitle(jmark) + 
    xlab("") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  strip.text.x = element_text(size = 5))
  print(m)
}

# Fit across plates --------------------------------------------------------



dat.fit.lst <- lapply(jmarksall, function(jmark){
  jsub <- subset(dat.merge.annot[[jmark]], !is.na(cluster))
  jsub$stype <- gsub(pattern = "LSK", replacement = "zLSK", x = jsub$stype)
  jsub$stype <- gsub(pattern = "Unenriched", replacement = "aUnenriched", x = jsub$stype)
  # jsub$experi <- as.character(frank(jsub$experi, ties.method = "dense"))
  jfit <- lm(formula = log2(chromocounts / spikeincounts) ~ 1 + experi + stype, data = jsub)
  jfit.null <- lm(formula = log2(chromocounts / spikeincounts) ~ 1 + experi, data = jsub)
  jsum <- anova(jfit.null, jfit)
  pval <- jsum$`Pr(>F)`[[2]]
  jfit.ds.lst <- DoFitsDownstream(jfit, jfit.null)
  jfit.ds.lst$jfit.merge$mark <- jmark
  jfit.ds.lst$jfit.merge$param <- gsub(pattern = "stype", replacement = "", x = jfit.ds.lst$jfit.merge$param)
  jfit.ds.lst$jfit.merge$param <- gsub(pattern = "\\(Intercept\\)", replacement = "aUnenriched", x = jfit.ds.lst$jfit.merge$param)
  # dat.fit <- data.frame(mark = jmark, pval = pval, t(data.frame(jfit$coefficients, stringsAsFactors = FALSE)), stringsAsFactors = FALSE)
  return(jfit.ds.lst$jfit.merge)
})


dat.fit.cluster.lst <- lapply(jmarksall, function(jmark){
  print(jmark)
  jsub <- subset(dat.merge.annot[[jmark]], !is.na(cluster))
  # jsub$experi <- as.character(frank(jsub$experi, ties.method = "dense"))
  jfit <- lm(formula = log2(chromocounts / spikeincounts) ~ 1 + experi + cluster, data = jsub)
  jfit.null <- lm(formula = log2(chromocounts / spikeincounts) ~ 1 + experi, data = jsub)
  jsum <- anova(jfit.null, jfit)
  pval <- jsum$`Pr(>F)`[[2]]
  jfit.ds.lst <- DoFitsDownstream(jfit, jfit.null)
  # dat.fit <- data.frame(mark = jmark, pval = pval, t(data.frame(jfit$coefficients, stringsAsFactors = FALSE)), stringsAsFactors = FALSE)
  jfit.ds.lst$jfit.merge$mark <- jmark
  jfit.ds.lst$jfit.merge$param <- gsub(pattern = "cluster", replacement = "", x = jfit.ds.lst$jfit.merge$param)
  jfit.ds.lst$jfit.merge$param <- gsub(pattern = "\\(Intercept\\)", replacement = "Granulocytes", x = jfit.ds.lst$jfit.merge$param)
  return(jfit.ds.lst$jfit.merge)
})


signif.factor <- 1.96
for (jmark in jmarksall){
  # print(jmark)
  # jsub.tmp <- dat.fit.lst[[jmark]]
  # jbaseline <- subset(jsub.tmp, stype == "Unenriched")``
  m <- ggplot(dat.fit.lst[[jmark]], aes(x = param, y = est, ymin = est - signif.factor * est.se, ymax = est + signif.factor * est.se)) + 
    geom_point() +
    geom_errorbar() + 
    theme_bw() + 
    ggtitle(jmark) + 
    ylab("log2(Cuts/Spikeins)") + 
    xlab("") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
}

for (jmark in jmarksall){
  print(jmark)
  m <- ggplot(dat.fit.cluster.lst[[jmark]], aes(x = forcats::fct_reorder(.f = param, .x = est, .desc = TRUE), y = est, ymin = est - signif.factor * est.se, ymax = est + signif.factor * est.se)) + 
    geom_point() +
    geom_errorbar() + 
    theme_bw() + 
    ggtitle(jmark) + 
    ylab("log2(Cuts/Spikeins)") + 
    xlab("") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
}


ymax <- 2.5
# plot relative to toi Granulocytes
for (jmark in jmarksall){
  print(jmark)
  jsub.tmp <- dat.fit.lst[[jmark]]
  jbaseline <- subset(jsub.tmp, param == "aUnenriched")$Estimate
  jsub.tmp$est <- jsub.tmp$est - jbaseline
  # jsub.tmp$est.se <- jsub.tmp$est.se - jbaseline
  m <- ggplot(jsub.tmp, aes(x = forcats::fct_reorder(.f = param, .x = est, .desc = TRUE), y = est, ymin = est - signif.factor * est.se, ymax = est + signif.factor * est.se)) + 
    geom_point() +
    geom_errorbar() + 
    theme_bw() + 
    ggtitle(jmark) + 
    ylab("log2(FC)") + 
    xlab("") + 
    ylim(c(-ymax, ymax)) + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    theme(aspect.ratio=1.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
}

# ymin <- 3
for (jmark in jmarksall){
  print(jmark)
  jsub.tmp <- dat.fit.cluster.lst[[jmark]]
  jbaseline <- subset(jsub.tmp, param == "Granulocytes")$Estimate
  jsub.tmp$est <- jsub.tmp$est - jbaseline
  m <- ggplot(jsub.tmp, aes(x = forcats::fct_reorder(.f = param, .x = est, .desc = TRUE), y = est, ymin = est - signif.factor * est.se, ymax = est + signif.factor * est.se)) + 
    geom_point() +
    geom_errorbar() + 
    theme_bw() + 
    ggtitle(jmark) + 
    ylab("log2(FC)") + 
    ylim(c(-ymax, ymax)) + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    xlab("") + 
    theme(aspect.ratio=1.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
}

dev.off()

# save objects
dat.merge.annot.out <- lapply(dat.merge.annot, function(jdat){
  jdat <- jdat %>%
    rowwise() %>%
    mutate(cluster = gsub(pattern = "zHSPCs", replacement = "HSPCs", cluster)) %>%
    mutate(cluster = gsub(pattern = "aGranulocytes", replacement = "Granulocytes", cluster)) %>%
    filter(!is.na(cluster))
})

for (jmark in jmarksall){
  print(jmark)
  outf <- paste0(outprefix, ".mark_", jmark, ".cell_cluster_tables.txt")
  jout.tmp <- dat.merge.annot.out[[jmark]]
  # make cluter second column
  print(length(jout.tmp))
  indx <- which(colnames(jout.tmp) == "cluster")
  assertthat::assert_that(length(indx) > 0)
  col.indx.vec <- seq(ncol(jout.tmp))
  col.indx.vec.append1 <- which(col.indx.vec %in% c(1, indx))
  col.indx.vec.append2 <- which(!col.indx.vec %in% c(1, indx))
  col.indx.vec.new <- c(col.indx.vec.append1, col.indx.vec.append2)
  jout.tmp.new <- jout.tmp[, col.indx.vec.new]
  # saveRDS(dat.merge.annot.out[[jmark]], file = outrds)
  print(dim(jout.tmp))
  print(dim(jout.tmp.new))
  fwrite(jout.tmp.new, file = outf, sep = "\t", na = "NA", quote = FALSE)
}
