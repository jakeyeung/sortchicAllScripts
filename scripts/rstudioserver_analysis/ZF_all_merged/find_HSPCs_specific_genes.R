# Jake Yeung
# Date of Creation: 2020-04-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/find_HSPCs_specific_genes.R
# Find HSPC-specific genes that we trust

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

# Load DE pseudobulk ------------------------------------------------------

# inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/rdata_robjs/WKM_pseudobulk_scrnaseq_downsampled.2020-04-27.EosinophilsKeep.AllGenes.RData"
inf.de <- "/home/jyeung/data/from_rstudioserver/zebrafish.poisson.SameNbrCells2.2020-04-30/WKM_pseudobulk_scrnaseq_downsampled.2020-04-30.SameNbrCells.RData"
load(inf.de, v=T)

pbulk.ctypefilt.long <- pbulk.ctypefilt.long %>%
  rowwise() %>%
  mutate(ens = strsplit(as.character(gene), "_")[[1]][[1]])

jgenesfull <- rownames(mat.pbulk.ds.ctypefilt)
jens <- sapply(jgenesfull, function(x) strsplit(x, "_")[[1]][[1]], USE.NAMES = FALSE)
jgenes <- sapply(jgenesfull, function(x) strsplit(x, "_")[[1]][[2]], USE.NAMES = FALSE)

g2e <- hash(jgenes, jens)

# load from Seurat
inf.seurat <- "/home/jyeung/data/from_rstudioserver/zebrafish.poisson.SameNbrCells2.2020-04-30/diff_exprs_Chloe_seurat.full.ctypefilt.SameNbrCells.rds"
dat.seurat <- readRDS(inf.seurat) %>%
  rowwise() %>%
  mutate(ens = strsplit(gene, "-")[[1]][[1]])


# Load pseudobulks ChIC  -------------------------------------------------------


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

pbulk.long.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.pbulk.chic <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/rdata_robjs/WKM_scChICseq_pseudobulk_downsampled.ens.", jmark, ".RData")
  load(inf.pbulk.chic, v=T)
  return(pbulk.long)
})

pbulk.chic.mat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.pbulk.chic <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/rdata_robjs/WKM_scChICseq_pseudobulk_downsampled.ens.", jmark, ".RData")
  load(inf.pbulk.chic, v=T)
  return(pbulk.filt.ds)
})


# Create one dataframe ----------------------------------------------------

ctypes.keep <- c("eryth", "HSC", "lymph", "monocyte")

# ctypes: eryth, HSC, lymph, monocytes
pbulk.k4me1 <- subset(pbulk.long.lst$H3K4me1) %>%
  ungroup() %>%
  mutate(mark = "H3K4me1") %>%
  filter(pbulk %in% ctypes.keep)

pbulk.k4me3 <- subset(pbulk.long.lst$H3K4me3) %>%
  ungroup() %>%
  mutate(mark = "H3K4me3") %>%
  filter(pbulk %in% ctypes.keep)

pbulk.k27me3 <- subset(pbulk.long.lst$H3K27me3) %>%
  ungroup() %>%
  mutate(mark = "H3K27me3",
         pbulk = gsub("eryth2", "eryth", pbulk),
         pbulk = gsub("HSC2", "HSC", pbulk)) %>%
  filter(pbulk %in% ctypes.keep)

# combine it all?
pbulk.merge <- bind_rows(pbulk.k4me1, pbulk.k4me3, pbulk.k27me3) %>%
  mutate(pbulk = gsub("monocyte", "granulocyte", pbulk))

# recalculate log2FC and zscores...
pbulk.merge <- pbulk.merge %>%
  group_by(gene, ens, mark) %>%
  mutate(log2FC = log2cuts - mean(log2cuts),
         log2zscore = log2FC / sd(log2cuts))

pbulk.merge$mark <- factor(pbulk.merge$mark, levels = c("H3K4me1", "H3K4me3", "H3K27me3"))

ggplot(subset(pbulk.merge), aes(x = pbulk, y = log2cuts, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(pbulk.merge), aes(x = pbulk, y = log2FC, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(pbulk.merge), aes(x = pbulk, y = log2FC, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(pbulk.merge), aes(x = pbulk, y = log2zscore, fill = mark)) + 
  geom_boxplot() + geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Check HSPC genes from literature ----------------------------------------
# from gateid

jgenes.hspc <- c("meis1b", "pmp22b", "ahnak", "krt8", "anxa5b", "mrc1b", "pa2g4a", "npm1a", "ahcy", "adh5", "fabp3", "myb")
jens.hspc <- sapply(jgenes.hspc, AssignHash, g2e, null.fill = NA)


# from Kobayashi
# inf.kob <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Kobayashi_2019_zf_genelists/lymphoid_genes_list.txt"
# inf.kob <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Kobayashi_2019_zf_genelists/erythmyeloid_genes_list.txt"
inf.kob <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Kobayashi_2019_zf_genelists/hspc_genes_list.txt"
kobayashi.genes <- fread(inf.kob, header = FALSE)$V1
jens.hspc <- kobayashi.genes

# take random 
jens.hspc <- sample(unique(pbulk.merge$ens), size = 1000)

# from seurat

jens.hspc <- subset(dat.seurat, avg_logFC > 1 & p_val_adj < 0.01 & cluster == "erythrocytes")$ens
jens.hspc <- subset(dat.seurat, avg_logFC > 0.5 & p_val_adj < 0.001 & cluster == "HSPCs")$ens
jens.hspc <- subset(dat.seurat, avg_logFC > 1 & p_val_adj < 0.01 & cluster == "granulocytes")$ens
jens.hspc <- subset(dat.seurat, avg_logFC > 1 & p_val_adj < 0.01 & cluster == "lymphocytes")$ens

# Plot boxplots -----------------------------------------------------------

ggplot(subset(pbulk.merge, ens %in% jens.hspc), aes(x = pbulk, y = log2cuts, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(pbulk.merge, ens %in% jens.hspc), aes(x = pbulk, y = log2FC, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0, linetype = "dotted")

ggplot(subset(pbulk.merge, ens %in% jens.hspc), aes(x = pbulk, y = log2zscore, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0, linetype = "dotted")


# take random set of genes
# jens.random <- 


# Check other genes -------------------------------------------------------


# HSPCs
jgenes.choose <- c("meis1b", "pmp22b", "ahnak", "krt8", "anxa5b", "mrc1b", "pa2g4a", "npm1a", "ahcy", "adh5", "fabp3", "myb")
# add kobayashi?
jgenes.choose <- c(jgenes.choose, kobayashi.genes)
jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)

# lymphs
jgenes.choose <- c("pax5", "cd79a", "bhlhe40", "cd83", "cxcr4a", "cd74b", "cd74a", "CD37", "zfp36l1a")  # lymphs
jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)

# monocytes
jgenes.choose1 <- c("adam8a", "odc1", "lta4h", "thy1", "scpp8", "illr4", "timp2b", "mmp9", "mmp13a", "scinlb")  # monocytes
jgenes.choose2 <- c("cpa5", "lyz", "lect2l", "npsn", "sms", "abcb9", "ch25hl2", "papss2b", "hsd3b7", "cfd")  # neutros
jgenes.choose <- c(jgenes.choose1, jgenes.choose2)
jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)

# eryths
jgenes.choose <- c("rhag", "prdx2", "epor", "gata1a", "tspo", "slc4a1a", "sptb", "cahz", "ba1", "alas2", "epb41b", "nt5c2l1")
jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)

ggplot(subset(pbulk.merge, ens %in% jens.choose), aes(x = pbulk, y = log2cuts, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(pbulk.merge, ens %in% jens.choose), aes(x = pbulk, y = log2FC, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0, linetype = "dotted")

ggplot(subset(pbulk.merge, ens %in% jens.choose), aes(x = pbulk, y = log2zscore, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0, linetype = "dotted")


# Use DE genes from seurat, use pvalue but also the mean expression -------

print(as.character(unique(dat.seurat$cluster)))
jclst <- "erythrocytes"
jclst <- "granulocytes"
# jclst <- ""
jclst <- "HSPCs"

jens.hspc <- subset(dat.seurat, avg_logFC > 0 & p_val_adj < 0.0001 & cluster == jclst)

# keep only genes where the mean of non-HSPCs are below threshold

print(head(pbulk.ctypefilt.long))

pbulk.ctypefilt.long <- pbulk.ctypefilt.long %>%
  rowwise() %>%
  mutate(is.hspc = pbulk == jclst)

pbulk.hspc.or.not <- pbulk.ctypefilt.long %>%
  group_by(gene, ens, is.hspc) %>%
  summarise(log2p1counts.mean = mean(log2p1counts))

pbulk.not <- pbulk.hspc.or.not %>%
  group_by(gene, ens) %>%
  mutate(diff = log2p1counts.mean[[2]] - log2p1counts.mean[[1]]) %>%
  ungroup() %>%
  filter(!is.hspc) %>%
  mutate(gene = as.character(gene))
  


subset(pbulk.hspc.or.not, grepl("meis1b", gene))
subset(pbulk.hspc.or.not, grepl("krt8", gene))

subset(pbulk.not, grepl("meis1b", gene))
subset(pbulk.not, grepl("krt8", gene))

plot(density(pbulk.ctypefilt.long$log2p1counts))
abline(v = 3)

exprs.min <- 2.2
jdiff <- 1

jens.hspc.meanfilt <- subset(pbulk.not, ens %in% jens.hspc$ens & log2p1counts.mean < exprs.min & diff > jdiff)

print(dim(jens.hspc.meanfilt))

# plot exprs?

jens.choose <- jens.hspc.meanfilt$ens

ggplot(subset(pbulk.merge, ens %in% jens.choose), aes(x = pbulk, y = log2cuts, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jclst)

ggplot(subset(pbulk.merge, ens %in% jens.choose), aes(x = pbulk, y = log2FC, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  ggtitle(jclst)

ggplot(subset(pbulk.merge, ens %in% jens.choose), aes(x = pbulk, y = log2zscore, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  ggtitle(jclst)

# housekeeping genes? 
pbulk.ctypefilt.long.summary <- pbulk.ctypefilt.long %>%
  group_by(gene) %>%
  summarise(log2p1counts.mean = mean(log2p1counts)) %>%
  rowwise() %>%
  mutate(ens = strsplit(as.character(gene), "_")[[1]][[1]])

plot(density(pbulk.ctypefilt.long.summary$log2p1counts.mean))

jens.choose <- sample(pbulk.ctypefilt.long.summary$ens, 100)

jens.choose <- subset(pbulk.ctypefilt.long.summary, log2p1counts.mean < 3)$ens

jsub <- subset(pbulk.merge, ens %in% jens.choose)
ggplot(jsub, aes(x = pbulk, y = log2zscore, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0, linetype = "dotted")


# Granulocytes have high K4me1... is that normalization? ------------------

ggplot(subset(pbulk.merge), aes(x = pbulk, y = log2FC, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0, linetype = "dotted")

head(mat.pbulk.ds.ctypefilt)

boxplot(log2(mat.pbulk.ds.ctypefilt + 1))



mat.pbulk.ds.ctypefilt.log <- log2(pbulk.chic.mat.lst$H3K4me1 + 1)

cnames.keep <- which(colnames(mat.pbulk.ds.ctypefilt.log) %in% ctypes.keep)

mat.pbulk.ds.ctypefilt.log <- mat.pbulk.ds.ctypefilt.log[, cnames.keep]

library(preprocessCore)

mat.pbulk.ds.ctypefilt.log.qn <- preprocessCore::normalize.quantiles(mat.pbulk.ds.ctypefilt.log)
# mat.pbulk.ds.ctypefilt.log.qn <- mat.pbulk.ds.ctypefilt.log

colnames(mat.pbulk.ds.ctypefilt.log.qn) <- colnames(mat.pbulk.ds.ctypefilt.log)
rownames(mat.pbulk.ds.ctypefilt.log.qn) <- rownames(mat.pbulk.ds.ctypefilt.log)

jgenes <- sapply(rownames(mat.pbulk.ds.ctypefilt.log), function(x) strsplit(x, split = ";")[[1]][[2]], USE.NAMES = FALSE)
jens <- sapply(jgenes, AssignHash, g2e, null.fill = NA)

# now do log2fc
pbulk.ctypefilt.long.qn <- data.frame(gene = jgenes, ens = jens, mat.pbulk.ds.ctypefilt.log.qn, stringsAsFactors = FALSE) %>%
  reshape2::melt(., id.vars = c("gene", "ens"), variable.name = "pbulk", value.name = "logexprs") %>%
  group_by(gene)  %>%
  mutate(logfc = logexprs - mean(logexprs),
         zscore = logfc / sd(logexprs))

ggplot(pbulk.ctypefilt.long.qn, aes(x = pbulk, y = zscore, fill = pbulk)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# HSPCs
jgenes.choose <- c("meis1b", "pmp22b", "ahnak", "krt8", "anxa5b", "mrc1b", "pa2g4a", "npm1a", "ahcy", "adh5", "fabp3", "myb")
# add kobayashi?
jgenes.choose <- c(jgenes.choose, kobayashi.genes)
jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)


# monocytes
jgenes.choose1 <- c("adam8a", "odc1", "lta4h", "thy1", "scpp8", "illr4", "timp2b", "mmp9", "mmp13a", "scinlb")  # monocytes
jgenes.choose2 <- c("cpa5", "lyz", "lect2l", "npsn", "sms", "abcb9", "ch25hl2", "papss2b", "hsd3b7", "cfd")  # neutros
jgenes.choose <- c(jgenes.choose1, jgenes.choose2)
jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)



ggplot(pbulk.ctypefilt.long.qn %>% filter(ens %in% jens.choose), aes(x = pbulk, y = logfc, fill = pbulk)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(pbulk.ctypefilt.long.qn %>% filter(ens %in% jens.choose), aes(x = pbulk, y = zscore, fill = pbulk)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



