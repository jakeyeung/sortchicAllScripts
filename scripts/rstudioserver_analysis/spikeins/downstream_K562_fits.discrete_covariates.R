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



# Load inputs -------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmark <- "H3K27me3"

jsuffix <- "topn_5000.glmpcaout.penalty_5.by_plate.RData"
jprefix <- ".G1_G2_S."

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


# Load fits ---------------------------------------------------------------

# inf.fits <- paste0("/home/jyeung/data/from_rstudioserver/spikein_fits_cellcycle/fit_cellcycle_pseudotime.", jmark, ".rdata")
# inf.fits <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/spikein_fits_cellcycle/fit_cellcycle_pseudotime.", jmark, ".2020-08-14.RData"))
inf.fits <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/spikein_fits_cellcycle.discrete_covariate/fit_cellcycle_pseudotime.", jmark, ".2020-08-16.RData"))
load(inf.fits, v=T)


# Get slope  --------------------------------------------------------------

slopes.S <- lapply(jfits, function(jfit){
  return(coef(jfit)[["cellcycle.str1_S"]])
})

slopes.G2 <- lapply(jfits, function(jfit){
  return(coef(jfit)[["cellcycle.str2_G2/M"]])
})

pvals.S <- lapply(jfits, function(jfit){
  return(summary(jfit)$coefficients["cellcycle.str1_S", "Pr(>|z|)"])
})

pvals.G2 <- lapply(jfits, function(jfit){
  return(summary(jfit)$coefficients["cellcycle.str2_G2/M", "Pr(>|z|)"])
})

coefs.dat <- data.frame(coord = names(slopes.S), 
                        slope.S = unlist(slopes.S), pval.S = unlist(pvals.S), 
                        slope.G2 = unlist(slopes.G2), pval.G2 = unlist(pvals.G2), 
                        stringsAsFactors = FALSE)

ggplot(coefs.dat, aes(x = slope.S / log(2), y = -log10(pval.S))) + geom_vline(xintercept = 1) + 
  geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  coord_cartesian(xlim = c(-20, 20))

# add gene name

# hubprefix 
jmark.test <- "H3K27me3"
jinf.tss <- file.path(hubprefix, "jyeung/data/databases/refseq/GRCh38/GCF_000001405.39_GRCh38.p13_genomic.parsed.TSS.txt.merged.chromorenamed.4columns.bed")
jinf.test <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/K562_count_tables_50000.", jmark.test, jprefix, "rds"))
mat.test <- readRDS(jinf.test)


jchromos.keep <- paste("chr", c(seq(22), "X", "Y"), sep = "")
annot.out <- AnnotateCoordsFromList.GeneWise(coords.vec = rownames(mat.test), inf.tss = jinf.tss, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene, annodb = "org.Hs.eg.db", chromos.keep = jchromos.keep)
annot.out$out2.df$geneclean <- sapply(annot.out$out2.df$gene, function(x) strsplit(x, "\\.\\.")[[1]][[2]])



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


ggplot(coefs.dat %>% mutate(is.cc = geneclean %in% genes.keep), aes(x = slope.S / log(2), y = -log10(pval.S))) + geom_vline(xintercept = 1) + 
  geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  coord_cartesian(xlim = c(-20, 20)) + 
  facet_wrap(~is.cc)

ggplot(coefs.dat, aes(x = slope.S / log(2), y = -log10(pval.S))) + geom_vline(xintercept = 1) + 
  geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  coord_cartesian(xlim = c(-20, 20)) + 
  facet_wrap(~strategy)


ggplot(coefs.dat %>% filter(!is.na(strategy)), aes(y = slope.S / log(2), x = strategy)) +
  geom_boxplot() + geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  coord_cartesian(ylim = c(-20, 20)) 



# Is it stat signif?  -----------------------------------------------------


# maybe Fisher's exact teset 

coefs.dat <- coefs.dat %>%
  ungroup() %>%
  mutate(qval.S = p.adjust(pval.S, method = "BH"),
         qval.G2 = p.adjust(pval.G2, method = "BH"))

jthres <- 0.05
coefs.dat$is.signif.S <- sapply(coefs.dat$qval.S, function(x) x <= jthres)
coefs.dat$is.signif.G2 <- sapply(coefs.dat$qval.G2, function(x) x <= jthres)


ggplot(coefs.dat, aes(x = slope.S / log(2), y = -log10(pval.S))) + geom_vline(xintercept = 1) + 
  geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  coord_cartesian(xlim = c(-20, 20)) + 
  facet_wrap(~strategy)

ggplot(coefs.dat, aes(x = slope.S / log(2), y = -log10(pval.S), color = is.signif.S)) + geom_vline(xintercept = 1) + 
  geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  coord_cartesian(xlim = c(-20, 20)) + 
  facet_wrap(~strategy)

ggplot(coefs.dat %>% filter(strategy == "B"), aes(x = slope.S / log(2), y = -log10(pval.S), color = is.signif.S)) + geom_vline(xintercept = 1) + 
  geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  coord_cartesian(xlim = c(-20, 20)) + 
  facet_wrap(~strategy)

ggplot(coefs.dat %>% filter(strategy == "C"), aes(x = slope.S / log(2), y = -log10(pval.S), color = is.signif.S)) + geom_vline(xintercept = 1) + 
  geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  coord_cartesian(xlim = c(-20, 20)) + 
  facet_wrap(~strategy)

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


ggplot(fisher.G2.dat, aes(x = strat, y = -log10(pval + 0.01))) + geom_col() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("G2 enrichment")

ggplot(fisher.S.dat, aes(x = strat, y = -log10(pval + 0.01))) + geom_col() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("S enrichment")




