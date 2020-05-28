# Jake Yeung
# Date of Creation: 2020-05-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/get_DE_genes_like_in_BM_no_lymphocytes.R
# Lymphocytes are small population in chic, try removing them 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(hash)

exprs.max <- 5
logfc.min <- 2
logfc.max <- 1
exprs.thres.hk <- exprs.max
exprs.thres.ne <- exprs.max * 0.5

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/pseudobulk_analysis.get_DE_genes_from_pbulk_scrnaseq.nolymph"
dir.create(outdir)
outname <- paste0("de_genes_sorted_and_giladi.NoLymph.WithHouseKeepAndNotExpressed.FixExprsOther.exprsmax_", exprs.max, ".logfcmin_", logfc.min, ".logfcmax_", logfc.max)
outf <- file.path(outdir, paste0(outname, ".RData"))
outpdf <- file.path(outdir, paste0(outname, ".pdf"))


pdf(outpdf, width = 1020/72, height = 815/72, useDingbats = FALSE)


# Load DE pseudobulk ------------------------------------------------------

# inf.de <- "/home/jyeung/data/from_rstudioserver/zebrafish.poisson.SameNbrCells2.2020-04-30/WKM_pseudobulk_scrnaseq_downsampled.2020-04-30.SameNbrCells.RData"
inf.de <- "/home/jyeung/data/from_rstudioserver/zebrafish.poisson.SameNbrCells2.NoLymph.2020-05-05/WKM_pseudobulk_scrnaseq_downsampled.2020-05-05.SameNbrCells.NoLymph.RData"
load(inf.de, v=T)

pbulk.ctypefilt.long <- pbulk.ctypefilt.long %>%
  rowwise() %>%
  mutate(ens = strsplit(as.character(gene), "_")[[1]][[1]])

jgenesfull <- rownames(mat.pbulk.ds.ctypefilt)
jens <- sapply(jgenesfull, function(x) strsplit(x, "_")[[1]][[1]], USE.NAMES = FALSE)
jgenes <- sapply(jgenesfull, function(x) strsplit(x, "_")[[1]][[2]], USE.NAMES = FALSE)

g2e <- hash(jgenes, jens)


# Set cutoffs -------------------------------------------------------------

# ggplot(pbulk.ctypefilt.long, aes(x = pbulk, y = log2p1counts)) + 
#   geom_boxplot() + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(pbulk.ctypefilt.long, aes(x = log2p1counts, fill = pbulk)) + 
  facet_wrap(~pbulk, ncol = 1) + 
  geom_density() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = c(exprs.thres.ne, exprs.max))

jclsts <- as.character(sort(unique(pbulk.ctypefilt.long$pbulk)))
names(jclsts) <- jclsts

# set cutoffs for each jclust? 

de.ens.sorted.stringent <- lapply(jclsts, function(jclst){
  print(jclst)
  jsub <- pbulk.ctypefilt.long %>% 
    rowwise() %>%
    mutate(is.ctype = pbulk == jclst) %>%
    group_by(ens, is.ctype) %>%
    summarise(jexprs = mean(log2p1counts)) %>%
    group_by(ens) %>%
    mutate(logFC = jexprs[[2]] - jexprs[[1]],
           jexprs.other = jexprs[[1]])  %>%
    ungroup() %>%
    arrange(desc(logFC)) %>%
    filter(jexprs.other <= exprs.max & logFC >= logfc.min)
  return(unique(jsub$ens))
})

lapply(de.ens.sorted.stringent, length)


# plot outputs

# define house keeping?
jsub.hk <- pbulk.ctypefilt.long %>% 
  rowwise() %>%
  group_by(ens) %>%
  summarise(exprs.mean = mean(log2p1counts),
            exprs.sd = sd(log2p1counts)) %>%
  arrange(exprs.mean) %>%
  filter(exprs.mean >= exprs.thres.hk & exprs.sd <= logfc.max * 0.2)
genes.hk <- unique(jsub.hk$ens)
print(length(genes.hk))

jsub.ne <- pbulk.ctypefilt.long %>% 
  rowwise() %>%
  group_by(ens) %>%
  summarise(exprs.mean = mean(log2p1counts),
            exprs.sd = sd(log2p1counts)) %>%
  arrange(desc(exprs.mean)) %>%
  filter(exprs.mean <= exprs.thres.ne & exprs.sd <= logfc.max * 0.2)
genes.ne <- unique(jsub.ne$ens)
print(length(genes.ne))

de.ens.sorted.stringent$HighExprs <- genes.hk
de.ens.sorted.stringent$LowExprs <- genes.ne

lapply(de.ens.sorted.stringent, length)

jclsts <- as.list(jclsts)

jclsts$HighExprs <- "HighExprs"
jclsts$LowExprs <- "LowExprs"

# Plot outputs ------------------------------------------------------------

m.lst <- lapply(jclsts, function(jclst){
  jens <- de.ens.sorted.stringent[[jclst]]
  jsub <- subset(pbulk.ctypefilt.long, ens %in% jens)
  m <- ggplot(jsub, aes(y = zscore, x = pbulk, fill = pbulk)) + 
    geom_boxplot() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jclst)
  return(m)
})
print(m.lst)

dev.off()

# Save outputs ------------------------------------------------------------

# rename
de.ens.zf.stringent <- de.ens.sorted.stringent
pbulk.zf.ctypefilt.long <- pbulk.ctypefilt.long
g2e.zf <- g2e

save(de.ens.zf.stringent, pbulk.zf.ctypefilt.long, g2e.zf, file = outf)



