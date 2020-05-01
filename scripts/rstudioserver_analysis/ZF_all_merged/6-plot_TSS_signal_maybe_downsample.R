# Jake Yeung
# Date of Creation: 2020-04-29
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/6-plot_TSS_signal_maybe_downsample.R
# Plot TSS signal in zebrafish 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(DropletUtils)
library(JFuncs)

library(topicmodels)
library(hash)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)
library(GenomicRanges)
library(forcats)

# Marks -------------------------------------------------------------------

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# jvarcutoffs <- c(0.75, 2, 1, 0.5)
# names(jvarcutoffs) <- jmarks
# hubprefix <- "/home/jyeung/hub_oudenaarden"
# inmain <- paste0(hubprefix, "/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZF_AllMerged.winsize_50000.imputevarfilt.lessstringent")
# assertthat::assert_that(dir.exists(inmain))

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/TSS_signal.stringentfilt2"
dir.create(outdir)


for (jmark in jmarks){
  outname <- paste0("TSS_signal_from_5kb_bins.", jmark, ".pdf")
  outpdf <- file.path(outdir, outname)
  pdf(outpdf, useDingbats = FALSE)
  
  
  # jmark <- "H3K4me3"
  # jmark <- "H3K27me3"
  print(jmark)
  
  # Load pseudobulk ---------------------------------------------------------
  
  inf.pbulk <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/rdata_robjs/WKM_scChICseq_pseudobulk_downsampled.ens.", jmark, ".RData")
  load(inf.pbulk, v=T)  # loads annot.louv an dannot.glmpca.filt
  
  
  # filter by previously defined celltypes? 
  inf.annot.louv <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/LDA_downstream/LDA_downstream_ZF.2020-04-23.imputevarfilt.lessstringent/ZF_LDA_output.", jmark, ".keepn_150.final.ClusterTables.txt")
  assertthat::assert_that(file.exists(inf.annot.louv))
  
  
  inf.annot.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt")
  assertthat::assert_that(file.exists(inf.annot.glmpca))
  
  annot.louv <- fread(inf.annot.louv)
  annot.louv$clusterplate <- paste(annot.louv$cluster, annot.louv$plate, "_")
  
  annot.glmpca <- fread(inf.annot.glmpca)
  annot.glmpca.filt <- subset(annot.glmpca, cell %in% annot.louv$cell) %>%
    rowwise() %>%
    mutate(clusterplate = paste(cluster, plate, sep = "_")) %>%
    mutate(cluster = ifelse(cluster %in% c("lymph1", "lymph2"), "lymph", cluster)) %>%   # rename lymph1 and lymph2 into lymph
    ungroup() %>%
    filter(cluster != "Unknown") %>%
    left_join(., subset(annot.louv, select = c(cell, var.imputed)))
  print("annot glmpca filt")
  print(annot.glmpca.filt)
  
  
  
  # plot uMAP 
  
  cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
  m <- ggplot(annot.glmpca.filt, aes(x = umap1, y = umap2, color = cluster)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette)
  print(m)
  
  # # load LDA
  # jvarcutoff <- jvarcutoffs[[jmark]]
  # infname <- paste0("lda_outputs.counts_table_var_filt.", jmark, ".imputevar_", jvarcutoff, ".K-30.binarize.FALSE/ldaOut.counts_table_var_filt.", jmark, ".imputevar_", jvarcutoff, ".K-30.Robj")
  # inf <- file.path(inmain, infname)
  # assertthat::assert_that(file.exists(inf))
  # 
  # load(inf, v=T)
  # count.mat <- as.matrix(count.mat)
  # tm.result <- posterior(out.lda)
  # colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = "_")
  # rownames(tm.result$terms) <- paste("topic", rownames(tm.result$terms), sep = "_")
  # topics.mat <- tm.result$topics
  
  
  # Load count TSS data  ----------------------------------------------------
  
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.TSS.CodingOnly.imputevarfilt.lessstringent.mapq_40.winsize_10000/", jmark, ".imputevarfilt.lessstringent.mapq_40.remerged.countTable.TSS.csv")
  mat <- ReadMatTSSFormat(inf)
  
  
  
  # head(pbulk.long)
  
  # Calculate entropy  ------------------------------------------------------
  
  # CalculateEntropy()
  pbulk.long <- pbulk.long %>%
    group_by(pbulk) %>%
    mutate(ncuts.frac = ncuts / sum(ncuts))
  
  pbulk.H <- pbulk.long %>%
    group_by(pbulk) %>%
    summarise(H = CalculateEntropy(ncuts, normalize.p = TRUE))
  
  m <- ggplot(pbulk.long, aes(x = ncuts.frac)) + geom_density() + facet_wrap(~pbulk, ncol = 1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  # Do it single cells ------------------------------------------------------
  
  genes.keep <- unique(pbulk.long$genefull)  
  cells.keep <- annot.glmpca.filt$cell
  
  mat.filt <- mat[genes.keep, cells.keep]
  
  mat.filt.long <- data.frame(genefull = rownames(mat.filt), as.data.frame(as.matrix(mat.filt)), stringsAsFactors = FALSE) %>%
    reshape2::melt(., variable.name = "cell", value.name = "ncuts") %>%
    group_by(cell) %>%
    mutate(ncuts.frac = ncuts / sum(ncuts))
  
  mat.filt.sum <- mat.filt.long %>%
    group_by(cell) %>%
    summarise(H = CalculateEntropy(ncuts.frac, normalize.p=TRUE)) 
  
  mat.filt.sum <- left_join(mat.filt.sum, annot.glmpca.filt %>% mutate(cell = make.names(cell)))
  
  m <- ggplot(mat.filt.sum, aes(x = umap1, y = umap2, color = H)) + facet_wrap(~plate) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point() + scale_color_viridis_c() + ggtitle("Entropy from chic TSS")
  print(m)
  
  
  # Try reading from a bin file ?  ------------------------------------------
  
  
  inf.bins <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.winsize_5000/PZ-ChIC-ZF_", jmark, "_2020-04-07.countTable.csv")
  mat.bins <- ReadMatSlideWinFormat(inf.bins)
  
  # filter cells
  mat.bins <- mat.bins[, cells.keep]
  
  # # downsample?
  # jprop <- 5000 / colSums(mat.bins)
  # jprop <- sapply(jprop, function(x) ifelse(x > 1, 1, x))
  # mat.bins <- DropletUtils::downsampleMatrix(mat.bins, prop = jprop)
  
  bins.5kb <- rownames(mat.bins)
  
  
  # jwin <- 10000
  # hubprefix <- "/home/jyeung/hub_oudenaarden"
  # inf.annot <- file.path(hubprefix, paste0("jyeung/data/databases/gene_tss/gene_tss.CodingOnly.winsize_", jwin, ".species_drerio.bed"))
  inf.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/zebrafish/gene_tss.CodingOnly.winsize_10000.species_drerio.bed"
  assertthat::assert_that(file.exists(inf.annot))
  
  jchromos <- paste("chr", seq(25), sep = "")
  
  annot.out <- AnnotateCoordsFromList(bins.5kb, inf.tss = inf.annot, txdb = TxDb.Drerio.UCSC.danRer11.refGene, annodb = "org.Dr.eg.db", chromos.keep = jchromos)
  
  # assign regions to either promoter or not 
  
  # creaete hash or coord to distance??
  
  # dhash <- hash(annot.out$regions.annotated$region_coord, annot.out$regions.annotated$distanceToTSS)
  annot.sub <- subset(annot.out$regions.annotated, select = c(distanceToTSS, region_coord, ENSEMBL, SYMBOL, GENENAME))
  
  # merge with count matrix?
  # create a TSS and a nonTSS matrix
  bins.tss <- subset(annot.sub, abs(distanceToTSS) <= 1000)
  bins.others <- subset(annot.sub, abs(distanceToTSS) > 100000)
  
  mat.bins.tss <- mat.bins[bins.tss$region_coord, ]
  mat.bins.others <- mat.bins[bins.others$region_coord, ]
  
  # get total counts
  ncuts.total <- data.frame(cell = colnames(mat.bins), ncuts.total = colSums(mat.bins), stringsAsFactors = FALSE)
  ncuts.tss <- data.frame(cell = colnames(mat.bins.tss), ncuts.tss = colSums(mat.bins.tss), stringsAsFactors = FALSE)
  ncuts.others <- data.frame(cell = colnames(mat.bins.others), ncuts.others = colSums(mat.bins.others), stringsAsFactors = FALSE)
  
  ncuts.merge <- Reduce(left_join, x = list(ncuts.tss, ncuts.others), init = ncuts.total) %>%
    rowwise() %>%
    mutate(ncuts.tss.frac = ncuts.tss / ncuts.total) %>%
    left_join(., annot.glmpca.filt)
  print(ncuts.merge)
  
  
  m <- ggplot(ncuts.merge, aes(x = ncuts.tss, y = ncuts.others)) + 
    geom_point() +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~plate) + geom_abline(slope = 1) + ggtitle(jmark)
  print(m)
  
  m <- ggplot(ncuts.merge, aes(x = umap1, y = umap2, color = ncuts.tss.frac)) + 
    geom_point() + scale_color_viridis_c() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~plate) + ggtitle(jmark)
  print(m)
  
  m <- ggplot(ncuts.merge, aes(x = umap1, y = umap2, color = ncuts.tss.frac)) + 
    geom_point() + scale_color_viridis_c() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~cluster) + ggtitle(jmark)
  print(m)
  
  m <- ggplot(ncuts.merge, aes(x = umap1, y = umap2, color = ncuts.tss.frac)) + 
    geom_point() + scale_color_viridis_c() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark)
  print(m)
  
  m <- ggplot(ncuts.merge, aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + scale_color_manual(values = cbPalette) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~plate) + ggtitle(jmark)
  print(m)
  
  m <- ggplot(ncuts.merge %>% filter(cluster != "eryth"), aes(x = umap1, y = umap2, color = ncuts.tss.frac)) + 
    geom_point() + scale_color_viridis_c() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~plate) + ggtitle(jmark)
  print(m)
  
  m <- ggplot(ncuts.merge, aes(x = umap1, y = umap2, color = log10(ncuts.tss))) + 
    geom_point() + scale_color_viridis_c() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~plate) + ggtitle(jmark)
  print(m)
  
  m <- ggplot(ncuts.merge, aes(x = umap1, y = umap2, color = log10(ncuts.total))) + 
    geom_point() + scale_color_viridis_c() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~plate) + ggtitle(jmark)
  print(m)
  
  m <- ggplot(ncuts.merge, aes(x = umap1, y = umap2, color = log10(ncuts.others))) + 
    geom_point() + scale_color_viridis_c() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~plate) + ggtitle(jmark)
  print(m)
  
  # dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
  # jchromos <- paste("chr", seq(25), sep = "")
  # dat.var <- CalculateVarAll(dat.impute.log, jchromos) %>%
  #   left_join(., annot.glmpca.filt)
  m.var <- ggplot(annot.glmpca.filt, aes(x = umap1, y = umap2, color = var.imputed)) + 
    scale_color_viridis_c(direction = -1) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.var)
  
  
  # Average across celltypes?  ----------------------------------------------
  
  ncuts.merge.ctype <- ncuts.merge %>%
    group_by(cluster) %>%
    summarise(ncuts.tss.frac.mean = mean(ncuts.tss.frac),
              ncuts.tss.frac.sd = sd(ncuts.tss.frac))
  
  dat.var.sum <- annot.glmpca.filt %>%
    group_by(cluster) %>%
    summarise(var.impute.mean = mean(var.imputed),
              var.impute.sd = sd(var.imputed))
  
  m <- ggplot(ncuts.merge.ctype, aes(x = cluster, y = ncuts.tss.frac.mean, ymin = ncuts.tss.frac.mean - ncuts.tss.frac.sd, ymax = ncuts.tss.frac.mean + ncuts.tss.frac.sd)) + 
    geom_col() + geom_errorbar(width = 0.25) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark)
  print(m)
  
  m.var.sum <- ggplot(dat.var.sum, aes(x = cluster, y = var.impute.mean, ymin = var.impute.mean - var.impute.sd, ymax = var.impute.mean + var.impute.sd)) + 
    geom_col() + geom_errorbar(width = 0.25) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark)
  print(m.var.sum)
  
  # plot TSS fraction versus variance 
  m.compare <- ggplot(ncuts.merge, aes(x = ncuts.tss.frac, y = var.imputed, color = cluster)) + geom_point() +
    scale_color_manual(values = cbPalette) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(jmark)
  print(m.compare)
  m.compare <- ggplot(ncuts.merge, aes(x = ncuts.others / ncuts.total, y = var.imputed, color = cluster)) + geom_point() +
    scale_color_manual(values = cbPalette) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(jmark)
  print(m.compare)
  m.compare <- ggplot(ncuts.merge, aes(x = ncuts.tss, y = var.imputed, color = cluster)) + geom_point() +
    scale_color_manual(values = cbPalette) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(jmark)
  print(m.compare)
  m.compare <- ggplot(ncuts.merge, aes(x = ncuts.others, y = var.imputed, color = cluster)) + geom_point() +
    scale_color_manual(values = cbPalette) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(jmark)
  print(m.compare)
  
  
  dev.off()
}








