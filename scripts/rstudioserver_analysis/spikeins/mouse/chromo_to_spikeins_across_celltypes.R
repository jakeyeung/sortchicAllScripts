# Jake Yeung
# Date of Creation: 2020-09-13
# File: 
# 

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


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Do boxplots of chromo to spikeincounts ----------------------------------


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_mouse.round2/chromo_spikein_counts"
# dir.create(outdir)
pdf(file.path(outdir, paste0("chromo_spikein_counts.allmark.", Sys.Date(), ".pdf")))



# Load cells keep ---------------------------------------------------------

# hubprefix <- "/home/jyeung/hub_oudenaarden"

inmain <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_VAN5046_varfilt"
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

jvarfilt <- 1
jsuffix <- paste0("_l2r_filt.2020-09-13.minl2r_-1.varfilt_", jvarfilt, ".K-30")


ldas.lst <- lapply(jmarks, function(jmark){
  dname <- paste0("lda_outputs.count_mat_", jmark, jsuffix, ".binarize.FALSE")
  indir <- file.path(inmain, dname)
  fname <- paste0("ldaOut.count_mat_", jmark, jsuffix, ".Robj")
  inf.lda <- file.path(indir, fname)
  assertthat::assert_that(file.exists(inf.lda))
  print(inf.lda)
  load(inf.lda, v=T)
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)
  dat.umap$mark <- jmark
  return(list(dat.umap = dat.umap, tm.result = tm.result))
})


dat.umap.long <- lapply(jmarks, function(jmark){
  ldas.lst[[jmark]]$dat.umap
}) %>%
  bind_rows()

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var.long <- lapply(ldas.lst, function(jdat){
  dat.impute <- t(log2(jdat$tm.result$topics %*% jdat$tm.result$terms))
  dat.var <- CalculateVarAll(dat.impute, jchromos)
  return(dat.var)
}) %>%
  bind_rows()

dat.umap.long <- left_join(dat.umap.long, dat.var.long, by = "cell")

# Load sipkeins -----------------------------------------------------------

inf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.chromo2spikeinfilt.VAN5046/spikein_info.txt"
dat.spikeins.mat <- fread(inf.spikeins)

# dat.merge <- dat.spikeins.mat %>%
dat.merge <- left_join(dat.umap.long, subset(dat.spikeins.mat, select = c(samp, spikeincounts, chromocounts)), by = c("cell" = "samp")) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"), 
         plate = as.numeric(strsplit(experi, "-")[[1]][[6]]),
         rowcoord = AddPlateCoordinates(cell)$rowcoord,
         colcoord = AddPlateCoordinates(cell)$colcoord,
         stype = AnnotateSortFromLayout(plate, rowcoord, colcoord))

# statistically significant? 


# Load cells keep ---------------------------------------------------------



dat.fit.lst <- lapply(jmarks, function(jmark){
  jsub <- subset(dat.merge, mark == jmark)
  jsub$stype <- gsub(pattern = "LSK", replacement = "aLSK", x = jsub$stype)
  jsub$plate <- as.character(frank(jsub$plate, ties.method = "dense"))
  jfit <- lm(formula = log2(chromocounts / spikeincounts) ~ 1 + plate + stype, data = jsub)
  jfit.null <- lm(formula = log2(chromocounts / spikeincounts) ~ 1 + plate, data = jsub)
  jsum <- anova(jfit.null, jfit)
  pval <- jsum$`Pr(>F)`[[2]]
  dat.fit <- data.frame(mark = jmark, pval = pval, t(data.frame(jfit$coefficients, stringsAsFactors = FALSE)), stringsAsFactors = FALSE)
  return(dat.fit)
})


for (jmark in jmarks){
  jsub <- dat.merge %>% filter(mark == jmark)
  jsub$stype <- factor(jsub$stype, levels = c("LSK", "LinNeg", "Unenriched"))
  jpval <- signif(subset(dat.fit.lst[[jmark]])$pval, digits = 2)
  m <- ggplot(jsub, aes(x = stype, y = log2(chromocounts / spikeincounts))) + 
    geom_boxplot() + 
    geom_jitter(alpha = 0.33, width = 0.12)  +  
    ggtitle(jmark, paste("pval:", jpval)) + 
    facet_wrap(~plate) + 
    scale_color_viridis_c() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}


for (jmark in jmarks){
  jsub <- dat.merge %>% filter(mark == jmark)
  jsub$stype <- factor(jsub$stype, levels = c("LSK", "LinNeg", "Unenriched"))
  jpval <- signif(subset(dat.fit.lst[[jmark]])$pval, digits = 2)
  m <- ggplot(jsub, aes(x = stype, y = log2(chromocounts))) + 
    geom_boxplot() + 
    geom_jitter(alpha = 0.33, width = 0.12)  +  
    ggtitle(jmark) + 
    facet_wrap(~plate) + 
    scale_color_viridis_c() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}


for (jmark in jmarks){
  jsub <- dat.merge %>% filter(mark == jmark)
  jsub$stype <- factor(jsub$stype, levels = c("LSK", "LinNeg", "Unenriched"))
  jpval <- signif(subset(dat.fit.lst[[jmark]])$pval, digits = 2)
  m <- ggplot(jsub, aes(x = stype, y = cell.var.within.sum.norm)) + 
    geom_boxplot() + 
    geom_jitter(alpha = 0.33, width = 0.12)  +  
    ggtitle(jmark) + 
    facet_wrap(~plate) + 
    scale_color_viridis_c() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}

dev.off()
