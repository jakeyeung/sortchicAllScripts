# Jake Yeung
# Date of Creation: 2020-08-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/downstream_K562_fits.discrete_covariates.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)
library(GenomicRanges)



# Load inputs -------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

jsuffix <- "topn_5000.glmpcaout.penalty_5.by_plate.RData"
jprefix <- ".G1_G2_S."

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_K562_cellcycle_gene_sets"


# hubprefix 
jmark.test <- "H3K27me3"
jinf.tss <- file.path(hubprefix, "jyeung/data/databases/refseq/GRCh38/GCF_000001405.39_GRCh38.p13_genomic.parsed.TSS.txt.merged.chromorenamed.4columns.bed")
jinf.test <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/K562_count_tables_50000.", jmark.test, jprefix, "rds"))
mat.test <- readRDS(jinf.test)


jchromos.keep <- paste("chr", c(seq(22), "X", "Y"), sep = "")
annot.out <- AnnotateCoordsFromList.GeneWise(coords.vec = rownames(mat.test), inf.tss = jinf.tss, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene, annodb = "org.Hs.eg.db", chromos.keep = jchromos.keep)
annot.out$out2.df$geneclean <- sapply(annot.out$out2.df$gene, function(x) strsplit(x, "\\.\\.")[[1]][[2]])

slopemax <- 5

make.plots <- TRUE

jgenes <- c("PLK1", "UNG", "PCNA", "HPRT1", "CENPE", "KIF11", "UBE2C", "TOP2A", "KIF2C", "CDK1", "TTK")

for (jmark in jmarks){
  
  print(jmark)
  
  outf <- file.path(outdir, paste0(jmark, ".bin_by_bin_fits_downstream4.pdf"))
  
  if (make.plots){
    pdf(file = outf, useDingbats = FALSE)
  }
  
  # jmark <- "H3K27me3"
  # jmark <- "H3K9me3"
  
  
  # inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein/K562_count_tables_50000.", jmark, ".G1_G2_S.glmpcaout", jsuffix))
  inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein/K562_count_tables_50000.", jmark, jprefix, jsuffix))
  assertthat::assert_that(file.exists(inf.glmpca))
  
  
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/K562_count_tables_50000.", jmark, jprefix, "rds"))
  inf.spike <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/spikein_counts/spikein_counts_all.RData")
  
  assertthat::assert_that(file.exists(inf))
  assertthat::assert_that(file.exists(inf.spike))
  
  mat <- readRDS(inf)
  load(inf.spike, v=T)
  load(inf.glmpca, v=T)
  
  dat.spikeins.mat <- AddCellCycleLabel.bydat(dat.spikeins.mat)
  
  
  # Load fits ---------------------------------------------------------------
  
  # inf.fits <- paste0("/home/jyeung/data/from_rstudioserver/spikein_fits_cellcycle/fit_cellcycle_pseudotime.", jmark, ".rdata")
  # inf.fits <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/spikein_fits_cellcycle/fit_cellcycle_pseudotime.", jmark, ".2020-08-14.RData"))
  inf.fits <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/spikein_fits_cellcycle.discrete_covariate/fit_cellcycle_pseudotime.", jmark, ".2020-08-16.simplified.RData"))
  load(inf.fits, v=T)
  
  # 
  # # Get slope  --------------------------------------------------------------
  # 
  # slopes.S <- lapply(jfits, function(jfit){
  #   return(coef(jfit)[["cellcycle.str1_S"]])
  # })
  # 
  # slopes.G2 <- lapply(jfits, function(jfit){
  #   return(coef(jfit)[["cellcycle.str2_G2/M"]])
  # })
  # 
  # pvals.S <- lapply(jfits, function(jfit){
  #   return(summary(jfit)$coefficients["cellcycle.str1_S", "Pr(>|z|)"])
  # })
  # 
  # pvals.G2 <- lapply(jfits, function(jfit){
  #   return(summary(jfit)$coefficients["cellcycle.str2_G2/M", "Pr(>|z|)"])
  # })
  # 
  # coefs.dat <- data.frame(coord = names(slopes.S), 
  #                         slope.S = unlist(slopes.S), pval.S = unlist(pvals.S), 
  #                         slope.G2 = unlist(slopes.G2), pval.G2 = unlist(pvals.G2), 
  #                         stringsAsFactors = FALSE)
  
  m <- ggplot(coefs.dat, aes(x = slope.S / log(2), y = -log10(pval.S))) + geom_vline(xintercept = 1, color = 'blue') + 
    geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    coord_cartesian(xlim = c(-slopemax, slopemax)) + 
    ggtitle(jmark) 
  print(m)
  
  m <- ggplot(coefs.dat, aes(x = slope.G2 / log(2), y = -log10(pval.G2))) + geom_vline(xintercept = 1, color = 'blue') + 
    geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    coord_cartesian(xlim = c(-slopemax, slopemax)) + 
    ggtitle(jmark)  
  print(m)
  
  coefs.dat.long <- melt(coefs.dat %>% dplyr::select(c(coord, slope.S, slope.G2)), id.vars = "coord")
  colnames(coefs.dat.long) <- c("coord", "fit.param", "slope")
  
  m <- ggplot(coefs.dat.long, aes(x = slope / log(2), fill = fit.param)) + 
    geom_vline(xintercept = 1, color = 'blue') + 
    geom_vline(xintercept = 0, color = 'black', linetype = "dotted") + 
    geom_density(alpha = 0.3) + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    coord_cartesian(xlim = c(-slopemax, slopemax)) + 
    ggtitle(jmark) 
  print(m)
  
  m <- ggplot(coefs.dat.long, aes(x = slope / log(2), fill = fit.param)) + 
    geom_vline(xintercept = 1, color = 'blue') + 
    geom_vline(xintercept = 0, color = 'black', linetype = "dotted") + 
    geom_density(alpha = 0.3) + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    coord_cartesian(xlim = c(-slopemax, slopemax)) + 
    ggtitle(jmark) 
  print(m)
  
  # add gene name
  
  
  
  # Add RNA  ----------------------------------------------------------------
  
  
  # Load RNA data  ----------------------------------------------------------
  
  inf.rna <- file.path(hubprefix, "jyeung/data/public_data/TableS1.csv")
  
  dat.rna <- fread(inf.rna, skip = 1)
  
  row2gene <- hash::hash(annot.out$out2.df$region_coord, annot.out$out2.df$geneclean)
  
  
  genes.keep <- subset(dat.rna, strategy_group == "B")$gene_symbol
  print(length(genes.keep))
  
  gene2strategy <- hash::hash(dat.rna$gene_symbol, dat.rna$strategy_group)
  row2strategy <- hash::hash(coefs.dat$coord, sapply(coefs.dat$coord, function(x) {
    jgene <- AssignHash(x, row2gene, null.fill = NA)
    jstrat <- AssignHash(x, gene2strategy, null.fill = NA)
    return(jstrat)
  })) 
  
  coefs.dat$geneclean <- sapply(coefs.dat$coord, function(x) AssignHash(x = x, row2gene, null.fill = x))
  coefs.dat$strategy <- sapply(coefs.dat$geneclean, function(x) AssignHash(x, gene2strategy, null.fill = NA))
  
  
  m <- ggplot(coefs.dat %>% mutate(is.cc = geneclean %in% genes.keep), aes(x = slope.S / log(2), y = -log10(pval.S))) + geom_vline(xintercept = 1, color = 'blue') + 
    geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
    coord_cartesian(xlim = c(-slopemax, slopemax)) + 
    facet_wrap(~is.cc) + 
    ggtitle(jmark)
  print(m)
  
  m <- ggplot(coefs.dat, aes(x = slope.S / log(2), y = -log10(pval.S))) + geom_vline(xintercept = 1, color = 'blue') + 
    geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
    coord_cartesian(xlim = c(-slopemax, slopemax)) + 
    facet_wrap(~strategy)  + 
    ggtitle(jmark)
  print(m)
  
  m <- ggplot(coefs.dat, aes(x = slope.G2 / log(2), y = -log10(pval.G2))) + geom_vline(xintercept = 1, color = 'blue') + 
    geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
    coord_cartesian(xlim = c(-slopemax, slopemax)) + 
    facet_wrap(~strategy) + 
    ggtitle(jmark)
  print(m)
  
  m <- ggplot(coefs.dat %>% filter(!is.na(strategy)), aes(y = slope.S / log(2), x = strategy)) +
    geom_boxplot() + geom_hline(yintercept = 0, linetype = "dotted") + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    coord_cartesian(ylim = c(-slopemax, slopemax))  + 
    ggtitle(jmark)
  print(m)
  
  m <- ggplot(coefs.dat %>% filter(!is.na(strategy)), aes(y = slope.G2 / log(2), x = strategy)) +
    geom_boxplot() + geom_hline(yintercept = 0, linetype = "dotted") + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    coord_cartesian(ylim = c(-slopemax, slopemax))  + 
    ggtitle(jmark)
  print(m)
  

  # Fit single gene  --------------------------------------------------------------------
  
  
  # jgene <- "PLK1"
  for (jgene in jgenes){
    
    jcoefs <- subset(coefs.dat, geneclean == jgene)
    if (nrow(jcoefs) == 0){
      next
    }
    jrow <- jcoefs$coord
    jvec <- mat[jrow, ]
    
    jdat <- data.frame(cuts = jvec, cell = names(jvec), stringsAsFactors = FALSE) %>%
      left_join(., dat.spikeins.mat, by = c("cell" = "samp"))
    # refit
    jrefit <- glm(formula = cuts ~ cellcycle.str + offset(log(spikeincounts)), family = "poisson", data = jdat)
    jci <- confint(jrefit)
    
    int.G1 <- jrefit$coefficients[[1]]
    int.S <- int.G1 + jrefit$coefficients[[2]]
    int.G2 <- int.G1 + jrefit$coefficients[[3]]
    
    int.G1.lower <- jci[1, 1]
    int.S.lower <- int.G1.lower + jci[2, 1]
    int.G2.lower <- int.G1.lower + jci[3, 1] 
    
    int.G1.higher <- jci[1, 2]
    int.S.higher <- int.G1.higher + jci[2, 2]
    int.G2.higher <- int.G1.higher + jci[3, 2] 
    
    jfit.dat <- data.frame(y = c(int.G1, int.S, int.G2), 
                           ymin = c(int.G1.lower, int.S.lower, int.G2.lower),
                           ymax = c(int.G1.higher, int.S.higher, int.G2.higher),
                           cellcycle.str = c("0_G1", "1_S", "2_G2/M"), stringsAsFactors = FALSE)
    
    m <- ggplot(jdat, aes(x = cellcycle.str, y = log(cuts / spikeincounts))) + 
      geom_point(alpha = 0.25) + 
      theme_bw(18) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      xlab("") + 
      geom_errorbar(mapping = aes(x = cellcycle.str, ymin = ymin, ymax = ymax), data = jfit.dat, inherit.aes = FALSE, width = 0.5) + 
      ggtitle(paste(jmark, jgene, jrow))
    print(m)
    
    
  }
  
  
  
  # Is it stat signif?  -----------------------------------------------------
  
  
  # maybe Fisher's exact teset 
  
  coefs.dat <- coefs.dat %>%
    ungroup() %>%
    mutate(qval.S = p.adjust(pval.S, method = "BH"),
           qval.G2 = p.adjust(pval.G2, method = "BH"))
  
  jthres <- 0.05
  coefs.dat$is.signif.S <- sapply(coefs.dat$qval.S, function(x) x <= jthres)
  coefs.dat$is.signif.G2 <- sapply(coefs.dat$qval.G2, function(x) x <= jthres)
  
  m <- ggplot(coefs.dat, aes(x = -log10(pval.S), y = -log10(qval.S))) + geom_point() +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  print(m)
  
  m <- ggplot(coefs.dat, aes(x = -log10(pval.G2), y = -log10(qval.G2))) + geom_point() +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark)
  print(m)
  
  m <- ggplot(coefs.dat, aes(x = slope.S / log(2), y = -log10(pval.S))) + geom_vline(xintercept = 1, color = 'blue') + 
    geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
    coord_cartesian(xlim = c(-slopemax, slopemax)) + 
    facet_wrap(~strategy) + 
    ggtitle(jmark)
  print(m)
  
  m <- ggplot(coefs.dat, aes(x = slope.S / log(2), y = -log10(pval.S), color = is.signif.S)) + geom_vline(xintercept = 1, color = 'blue') + 
    geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
    coord_cartesian(xlim = c(-slopemax, slopemax)) + 
    facet_wrap(~strategy) + 
    ggtitle(paste("qval thres:", jthres, jmark))
  print(m)
  
  m <- ggplot(coefs.dat, aes(x = slope.S / log(2), y = -log10(pval.S), color = is.signif.G2)) + geom_vline(xintercept = 1, color = 'blue') + 
    geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
    coord_cartesian(xlim = c(-slopemax, slopemax)) + 
    facet_wrap(~strategy) + 
    ggtitle(paste("qval thres:", jthres, jmark))
  print(m)
  
  m <- ggplot(coefs.dat %>% filter(strategy == "B"), aes(x = slope.S / log(2), y = -log10(pval.S), color = is.signif.S)) + geom_vline(xintercept = 1, color = 'blue') + 
    geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
    coord_cartesian(xlim = c(-slopemax, slopemax)) + 
    facet_wrap(~strategy) + 
    ggtitle(paste("qval thres:", jthres, jmark, "strategy B"))
  print(m)
  
  m <- ggplot(coefs.dat %>% filter(strategy == "C"), aes(x = slope.S / log(2), y = -log10(pval.S), color = is.signif.S)) + geom_vline(xintercept = 1, color = 'blue') + 
    geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
    coord_cartesian(xlim = c(-slopemax, slopemax)) + 
    facet_wrap(~strategy) + 
    ggtitle(paste("qval thres:", jthres, jmark))
  print(m)
  
  # count signifs 
  
  jstrats <- unique((coefs.dat %>% filter(!is.na(strategy)))$strategy)
  names(jstrats) <- jstrats
  
  fisher.S.out <- lapply(jstrats, function(jstrat){
    tp <- nrow(subset(coefs.dat, strategy == jstrat & is.signif.S))
    fp <- nrow(subset(coefs.dat, strategy == jstrat & !is.signif.S))
    tn <- nrow(subset(coefs.dat, is.na(strategy) & !is.signif.S))
    fn <- nrow(subset(coefs.dat, is.na(strategy) & is.signif.S))
    x <- matrix(c(tp, fp, fn, tn), nrow = 2)
    fisher.test(x)
  })
  
  fisher.G2.out <- lapply(jstrats, function(jstrat){
    tp <- nrow(subset(coefs.dat, strategy == jstrat & is.signif.G2))
    fp <- nrow(subset(coefs.dat, strategy == jstrat & !is.signif.G2))
    tn <- nrow(subset(coefs.dat, is.na(strategy) & !is.signif.G2))
    fn <- nrow(subset(coefs.dat, is.na(strategy) & is.signif.G2))
    x <- matrix(c(tp, fp, fn, tn), nrow = 2)
    fisher.test(x)
  })
  
  fisher.S.dat <- lapply(jstrats, function(jstrat){
    jout <- fisher.S.out[[jstrat]]
    pval <- jout$p.value
    or <- jout$estimate
    data.frame(pval = pval, or = or, strat = jstrat, cellcycle = "S", stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  
  fisher.G2.dat <- lapply(jstrats, function(jstrat){
    jout <- fisher.G2.out[[jstrat]]
    pval <- jout$p.value
    or <- jout$estimate
    data.frame(pval = pval, or = or, strat = jstrat, cellcycle = "G2/M", stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  
  
  m <- ggplot(fisher.G2.dat, aes(x = strat, y = -log10(pval))) + geom_col() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle("G2 enrichment", jmark)
  print(m)
  
  m <- ggplot(fisher.S.dat, aes(x = strat, y = -log10(pval))) + geom_col() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle("S enrichment", jmark)
  print(m)
  
  if (make.plots){
    dev.off()
  }
  
}




