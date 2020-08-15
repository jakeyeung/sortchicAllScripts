# Jake Yeung
# Date of Creation: 2020-08-14
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/downstream_K562_fits.R
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
inf.fits <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/spikein_fits_cellcycle/fit_cellcycle_pseudotime.", jmark, ".2020-08-14.RData"))
load(inf.fits, v=T)


# Get slope  --------------------------------------------------------------

slopes <- lapply(jfits, function(jfit){
  return(coef(jfit)[["pseudotime"]])
})

pvals <- lapply(jfits, function(jfit){
  return(summary(jfit)$coefficients["pseudotime", "Pr(>|z|)"])
})

coefs.dat <- data.frame(coord = names(slopes), slope = unlist(slopes), pval = unlist(pvals), stringsAsFactors = FALSE)

ggplot(coefs.dat, aes(x = slope / log(2), y = -log10(pval))) + geom_vline(xintercept = 1) + 
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


ggplot(coefs.dat %>% mutate(is.cc = geneclean %in% genes.keep), aes(x = slope / log(2), y = -log10(pval))) + geom_vline(xintercept = 1) + 
  geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  coord_cartesian(xlim = c(-20, 20)) + 
  facet_wrap(~is.cc)

ggplot(coefs.dat, aes(x = slope / log(2), y = -log10(pval))) + geom_vline(xintercept = 1) + 
  geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  coord_cartesian(xlim = c(-20, 20)) + 
  facet_wrap(~strategy)


ggplot(coefs.dat %>% filter(!is.na(strategy)), aes(y = slope / log(2), x = strategy)) +
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  coord_cartesian(ylim = c(-20, 20)) 
