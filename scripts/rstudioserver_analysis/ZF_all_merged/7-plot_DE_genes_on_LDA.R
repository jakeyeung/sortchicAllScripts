# Jake Yeung
# Date of Creation: 2020-04-29
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/7-plot_DE_genes_on_LDA.R
# Plot DE genes on LDA 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

# Load DE genes -----------------------------------------------------------

inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/from_data/zebrafish.poisson.2020-04-27/diff_exprs_Chloe_seurat.full.ctypefilt.rds"
dat.de <- readRDS(inf.de)


# Contants ----------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jvarcutoffs <- c(0.75, 2, 1, 0.5)
names(jvarcutoffs) <- jmarks
jmark <- "H3K4me1"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/DE_genes_on_UMAP"


for (jmark in jmarks){
  
  print(jmark)

  outpdf <- file.path(outdir, paste0(jmark, ".DE_genes_on_umap3.pdf")) 
  pdf(outpdf, useDingbats = FALSE)
  
  
  # # Load GLMPCA -------------------------------------------------------------
  # 
  # inf.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/GLMPCA_outputs.ZF_AllMerged.imputevarfilt.lessstringent/ZF_", jmark, ".AllMerged.ZF_AllMerged.imputevarfilt.lessstringent.GLMPCA_var_correction.mergebinsize_1000.binskeep_500.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.5.winsorize_TRUE.2020-04-29.RData")
  # load(inf.glmpca, v=T)
  # # do it on the GLMPCA? 
  # jsettings <- umap.defaults
  # jsettings$n_neighbors <- 30
  # jsettings$min_dist <- 0.1
  # jsettings$random_state <- 123
  # dat.umap.glm <- DoUmapAndLouvain(glm.out$factors, jsettings = jsettings)
  # 
  # glm.impute <- t(as.matrix(glm.out$factors) %*% as.matrix(t(glm.out$loadings)))
  
  # Load LDA output ---------------------------------------------------------
  
  
  print(jmark)
  
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
    filter(cluster != "Unknown")  # remove the small cluster Unknown
  
  annot.glmpca.filt <- left_join(annot.glmpca.filt, subset(annot.louv, select = c(cell, var.imputed)))
  
  # plot uMAP 
  
  cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
  m <- ggplot(annot.glmpca.filt, aes(x = umap1, y = umap2, color = cluster)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette)
  print(m)
  
  
  # Load LDA  ---------------------------------------------------------------
  
  hubprefix <- "/home/jyeung/hub_oudenaarden"
  inmain <- paste0(hubprefix, "/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZF_AllMerged.winsize_50000.imputevarfilt.lessstringent")
  assertthat::assert_that(dir.exists(inmain))
  
  jvarcutoff <- jvarcutoffs[[jmark]]
  infname <- paste0("lda_outputs.counts_table_var_filt.", jmark, ".imputevar_", jvarcutoff, ".K-30.binarize.FALSE/ldaOut.counts_table_var_filt.", jmark, ".imputevar_", jvarcutoff, ".K-30.Robj")
  # infname <- paste0("lda_outputs.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.Robj")
  inf <- file.path(inmain, infname)
  assertthat::assert_that(file.exists(inf))
  
  load(inf, v=T)
  count.mat <- as.matrix(count.mat)
  tm.result <- posterior(out.lda)
  colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = "_")
  rownames(tm.result$terms) <- paste("topic", rownames(tm.result$terms), sep = "_")
  topics.mat <- tm.result$topics
  
  
  # Plot variance? ----------------------------------------------------------
  
  dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
  jchromos <- paste("chr", seq(25), sep = "")
  dat.var <- CalculateVarAll(dat.impute.log, jchromos) %>%
    left_join(., annot.glmpca.filt)
  
  m <- ggplot(dat.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    scale_color_viridis_c(direction = -1) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark)
  print(m)
  
  m <- ggplot(annot.glmpca.filt, aes(x = umap1, y = umap2, color = var.imputed)) + 
    scale_color_viridis_c(direction = -1) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(jmark)
  print(m)
  
  dat.var.sum <- annot.glmpca.filt %>%
    group_by(cluster) %>%
    summarise(var.impute.mean = mean(var.imputed),
              var.impute.sd = sd(var.imputed))
  
  m.var.sum <- ggplot(dat.var.sum, aes(x = cluster, y = var.impute.mean, ymin = var.impute.mean - var.impute.sd, ymax = var.impute.mean + var.impute.sd)) + 
    geom_col() + geom_errorbar(width = 0.25) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark)
  print(m.var.sum)
  
  
  # Plot DE genes on the UMAP  ----------------------------------------------
  
  
  # annotate bins
  
  inf.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/zebrafish/gene_tss.CodingOnly.winsize_20000.species_drerio.bed"
  assertthat::assert_that(file.exists(inf.annot))
  
  jchromos <- paste("chr", seq(25), sep = "")
  
  annot.out <- AnnotateCoordsFromList(rownames(dat.impute.log), inf.tss = inf.annot, txdb = TxDb.Drerio.UCSC.danRer11.refGene, annodb = "org.Dr.eg.db", chromos.keep = jchromos)
  
  # assign regions to either promoter or not 
  
  # creaete hash or coord to distance??
  
  # dhash <- hash(annot.out$regions.annotated$region_coord, annot.out$regions.annotated$distanceToTSS)
  annot.sub <- subset(annot.out$regions.annotated, select = c(distanceToTSS, region_coord, ENSEMBL, SYMBOL, GENENAME))
  
  # get a DE gene
  
  print(unique(dat.de$cluster))
  
  for (jclst in as.character(unique(dat.de$cluster))){
    print(jclst)
    
    jsub <- subset(dat.de, cluster == jclst) %>% arrange(desc(avg_logFC)) %>% filter(p_val_adj < 0.01)
    print(jsub)
    jgenes <- subset(dat.de, cluster == jclst)$gene[1:50]
    ens <- sapply(jgenes, function(jgenefull) strsplit(jgenefull, "-")[[1]][[1]], simplify = TRUE, USE.NAMES = FALSE)
    print(ens)
    
    # get appropriate coordinate
    jsub.annot <- subset(annot.sub, ENSEMBL %in% ens & abs(distanceToTSS) < 10000)
    ens.matched <- unique(jsub.annot$ENSEMBL)
    regions_coord.matched <- unique(jsub.annot$region_coord)
    
    # get expression in UMAP 
    assertthat::assert_that(length(regions_coord.matched) != 0)
    if (length(regions_coord.matched) > 1){
      exprs.dat.lda <- data.frame(cell = colnames(dat.impute.log), exprs = colMeans(dat.impute.log[regions_coord.matched, ]))
      # rnames.keep <- which(rownames(glm.impute) %in% regions_coord.matched)
      # assertthat::assert_that(length(rnames.keep) > 0)
      # exprs.dat <- data.frame(cell = colnames(glm.impute), exprs = colMeans(glm.impute[rnames.keep, ]))
    } else {
      # exprs.dat <- data.frame(cell = colnames(dat.impute.log), exprs = colMeans(dat.impute.log[regions_coord.matched, ]))
      exprs.dat <- data.frame(cell = colnames(glm.impute), exprs = glm.impute[regions_coord.matched, ])
    }
    ngenes <- length(regions_coord.matched)
    
    # plojt on UMAP 
    jmerge.lda <- left_join(annot.glmpca.filt, exprs.dat.lda)
    m <- ggplot(jmerge.lda, aes(x = umap1, y = umap2, color = exprs)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      scale_color_viridis_c() + ggtitle(jmark, paste(jclst, "ngenes:", ngenes))
    print(m)
    
    # jmerge.glmpca <- left_join(dat.umap.glm, exprs.dat)
    # ggplot(jmerge.glmpca, aes(x = umap1, y = umap2, color = exprs)) + geom_point() +
    #   scale_color_viridis_c() + 
    #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
  }
  dev.off()
}


