# Jake Yeung
# Date of Creation: 2021-03-28
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/misc/check_H3K27me3_TSS_absolute_levels.R
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

library(topicmodels)
library(scchicFuncs)

jstart <- Sys.Date()

# Load objs  --------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

# jmark.test <- "H3K4me3"
jmark.test <- "H3K27me3"


inf.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep1rep2rep3reseq_varfilt_2020-12-15/lda_outputs.mat_H3K27me3_rep1rep2rep3reseq.TSS.K-30.binarize.FALSE/ldaOut.mat_H3K27me3_rep1rep2rep3reseq.TSS.K-30.Robj"
load(inf.lda, v=T)

# get imputed
tm.result <- posterior(out.lda)
dat.imputed.log2 <- log2(t(tm.result$topics %*% tm.result$terms))

inf.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq.with_old/mat_H3K27me3_rep1rep2rep3reseq.metadata.txt"
dat.meta <- fread(inf.meta)
pdcs.old <- subset(dat.meta, jrep == "rep1old" & cluster == "pDCs")
pdcs2basos <- sample(pdcs.old$cell, size = round(length(pdcs.old$cell) / 2))

dat.meta <- dat.meta %>%
  rowwise() %>%
  mutate(cluster = ifelse(cell %in% pdcs2basos, "Basophils", cluster))


# pdf(outpdf, useDingbats = FALSE)

tm.result.lst <- list()
tm.result.lst[[jmark.test]] <- tm.result

count.mat.raw.lst <- list()
count.mat.raw.lst[[jmark.test]] <- count.mat

inf.genes.k4me1 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/geneset_from_H3K4me1/geneset_H3K4me1_TSS_topics2.filt.colorcoded2.txt"
dat.genes.k4me1 <- fread(inf.genes.k4me1)

inf.genes <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/celltype_specific_genes_defined_by_K4me3_TSS.txt")
dat.genes <- fread(inf.genes)

dat.meta.lst <- list()
dat.meta.lst[[jmark.test]] <- dat.meta

# dat.meta.lst <- lapply(jmarks, function(jmark){
#   inf.meta <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/metadata_umap_celltype_cuts.", jmark, ".txt")
#   fread(inf.meta)
# })

dat.imputed.log.lst <- lapply(tm.result.lst, function(tm.result){
  log2(t(tm.result$topics %*% tm.result$terms))
})

dat.imputed.log.lst <- list()
dat.imputed.log.lst[[jmark.test]] <- dat.imputed.log2



# Make heatmap, fix batch effects -----------------------------------------


cells.keep <- dat.meta.lst[[jmark.test]]$cell

genes.keep.k4me1 <- dat.genes.k4me1$gene
genes.keep.k4me3 <- dat.genes$gene
genes.keep <- unique(c(genes.keep.k4me3, genes.keep.k4me1))

print(paste("Number of genes keeping", length(genes.keep)))

# rkeep <- genes.keep %in% rownames()

mat4hm <- dat.imputed.log.lst[[jmark.test]][genes.keep, cells.keep]
mat4hm.uniq <- dat.imputed.log.lst[[jmark.test]][unique(genes.keep), cells.keep]

# heatmap3::heatmap3(mat4hm, Rowv = NA, Colv = NA, ColSideColors = dat.meta.lst[[jmark.test]]$colorcode, scale = "row", RowSideColors = dat.genes$colorcode, revC = TRUE, main = paste0(jmark.test))



# Check heatmap inputs?  --------------------------------------------------

inf.heatmap <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/heatmap_pdfs_and_ordered_matrices/heatmap_ordered_with_labels.H3K27me3.2021-01-08.rearranged.RData"
load(inf.heatmap, v=T)

ctypes.arranged <-  c("Eryths", "Bcells", "NKs", "pDCs", "Granulocytes", "DCs", "Basophils", "HSPCs")
colors.arranged <- unique(rowsidecolors)
col2ctype <- hash::hash(colors.arranged, ctypes.arranged)

rnames.keep <- rownames(mat.adj.tmp)
rname2col <- hash::hash(rnames.keep, rowsidecolors)
rowsidectype <- sapply(rowsidecolors, function(x) AssignHash(x = x, jhash = col2ctype))
rname2ctype <- hash::hash(rnames.keep, rowsidectype)


# Check H3K27me3 levels in HSPCs ------------------------------------------

cells.hspcs <- subset(dat.meta, cluster == "HSPCs")$cell

count.hspcs <- rowSums(count.mat[, cells.hspcs])

dat.count.hspcs <- data.frame(tss = names(count.hspcs), count.raw = count.hspcs, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(is.blood.gene = tss %in% rnames.keep,
         ctype = AssignHash(x = tss, jhash = rname2ctype, null.fill = "Anonblood"))

dat.count.hspcs$count.norm <- dat.count.hspcs$count.raw / sum(dat.count.hspcs$count.raw)
dat.count.hspcs$count.norm.log <- log2(dat.count.hspcs$count.norm * 10^6 + 1)

jcutoff <- 4.5

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K27me3_absolute_levels_check/H3K27me3_abs_levels_output.", Sys.Date(), ".pdf")

pdf(outpdf, file = outpdf)

ggplot(dat.count.hspcs, aes(x = log2(count.norm * 10^6 + 1), fill = is.blood.gene)) + 
  geom_density(alpha = 0.33) + 
  theme_bw() + 
  geom_vline(xintercept = jcutoff, linetype = "dotted") + 
  xlab("log2(CPM + 1) H3K27me3") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  ggtitle("H3K27me3 levels in HSPCs: split genes by blood genes vs nonblood genes")

ggplot(dat.count.hspcs, aes(x = log2(count.norm * 10^6 + 1), fill = is.blood.gene)) + 
  facet_wrap(~ctype) + 
  geom_density(alpha = 0.33) + 
  theme_bw() + 
  xlab("log2(CPM + 1) H3K27me3") + 
  geom_vline(xintercept = jcutoff, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("H3K27me3 levels in HSPCs: annotate genes by blood cell type")

ggplot(dat.count.hspcs, aes(y = log2(count.norm * 10^6 + 1), x = is.blood.gene)) + 
  geom_boxplot() + 
  theme_bw() + 
  geom_hline(yintercept = jcutoff, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ylab("log2(CPM + 1) H3K27me3") + 
  ggtitle("H3K27me3 levels in HSPCs: split genes by blood genes vs nonblood genes")

ggplot(dat.count.hspcs, aes(y = log2(count.norm * 10^6 + 1), fill = is.blood.gene, x = ctype)) + 
  geom_boxplot() + 
  theme_bw() + 
  geom_hline(yintercept = jcutoff, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ylab("log2(CPM + 1) H3K27me3") + 
  ggtitle("H3K27me3 levels in HSPCs: annotate genes by blood cell type")

dev.off()

# 
# dat.count.hspcs.summary <- dat.count.hspcs %>%
#   rowwise() %>%
#   mutate(is.on = count.norm.log > jcutoff) %>%
#   group_by(is.on, is.blood.gene) %>%
#   summarise(ngenes = length(tss))
# 
# conting.tab <- matrix(data = dat.count.hspcs.summary$ngenes, nrow = 2, ncol = 2, byrow = FALSE)
# conting.tab
# 
# fisher.test(conting.tab)
