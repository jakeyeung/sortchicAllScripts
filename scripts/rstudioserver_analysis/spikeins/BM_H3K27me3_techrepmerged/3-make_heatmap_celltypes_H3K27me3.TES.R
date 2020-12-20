# Jake Yeung
# Date of Creation: 2020-12-10
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_techrepmerged/3-make_heatmap_celltypes_H3K27me3.TES.R
# Use TES of H3K27me3 it's probably better

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)



hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load meta  --------------------------------------------------------------

# inf.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/topics_loadings_H3K27me3_techrepmerged/topics_H3K27me3_peaks_topics.txt"
inf.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.H3K27me3_techrepmerged/BM_rep2_rep3reseq_H3K27me3.2020-12-10.txt"
dat.meta <- fread(inf.meta)



# Load genes --------------------------------------------------------------

inf.genes <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-09.from_LDA_topics.condensed.heatmap.famousgenes.keepn_400.refmark_H3K4me3.2020-12-09.txt")
genes.dat <- fread(inf.genes) %>%
  rowwise() %>%
  mutate(gsymbol = strsplit(gene, "\\.")[[1]][[4]],
         coord = paste("chr", strsplit(gene, ";")[[1]][[1]], sep = ""),
         tssname = strsplit(gene, ";")[[1]][[2]])


cells.keep <- dat.meta$cell
# Load TSS ----------------------------------------------------------------

indir.counts.tss <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/count_tables_from_TSS_rep2_rep3reseq_bams_together/dist_10000")

infs.counts.tss <- list.files(indir.counts.tss, pattern = "*.txt", full.names = TRUE)

count.mat.lst.tss <- lapply(infs.counts.tss, function(inf){
  print(inf)
  mat <- ReadMatTSSFormat(inf, add.coord = TRUE)
  ckeep <- colnames(mat) %in% cells.keep
  return(mat[, ckeep])
}) 

all.rnames.tss <- unique(unlist(lapply(count.mat.lst.tss, function(x) rownames(x))))
count.mat.tss <- cbind.fill.lst(count.mat.lst.tss, all.rnames = all.rnames.tss)

count.mat.norm.tss <- sweep(count.mat.tss, MARGIN = 2, STATS = colSums(count.mat.tss), FUN = "/")


jcoords <- paste("chr", sapply(all.rnames.tss, function(x) strsplit(x, ";")[[1]][[1]]), sep = "")
coord2gene <- hash::hash(jcoords, all.rnames.tss)

# Load count mat ----------------------------------------------------------


indir.counts.tes <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/count_tables_from_TES_rep2_rep3reseq_bams_together")

infs.counts.tes <- list.files(indir.counts.tes, pattern = "*.txt", full.names = TRUE)

count.mat.lst.tes <- lapply(infs.counts.tes, function(inf){
  print(inf)
  mat <- ReadMatTSSFormat(inf, add.coord = TRUE)
  ckeep <- colnames(mat) %in% cells.keep
  return(mat[, ckeep])
}) 

all.rnames.tes <- unique(unlist(lapply(count.mat.lst.tes, function(x) rownames(x))))
count.mat.tes <- cbind.fill.lst(count.mat.lst.tes, all.rnames = all.rnames.tes)


# gene legnth norm
gene.coords <- sapply(rownames(count.mat.tes), function(x) strsplit(x, ";")[[1]][[1]])
gene.starts <- as.numeric(sapply(gene.coords, function(g) GetStart(g)))
gene.ends <- as.numeric(sapply(gene.coords, function(g) GetEnd(g)))
gene.lengths <- gene.ends - gene.starts

count.mat.norm.tes <- sweep(count.mat.tes, MARGIN = 2, STATS = colSums(count.mat.tes), FUN = "/")
count.mat.norm.tes <- sweep(count.mat.norm.tes, MARGIN = 1, STATS = gene.lengths, FUN = "/")



jgene <- "Tead1" 

jgene <- "Tal1" 

jgenes.vec <- subset(genes.dat, jset == "HSPCs")$gsymbol
jgene <- paste(jgenes.vec, collapse = "|")

jgene <- "Meis1" 
jgene <- "Tal1" 
jgene <- "Ppp1r18" 


jgene <- "Ebf1" 
jgene <- "Hbb-y"
jgene <- "Pax5"
jgene <- "Hoxa9"
jgene <- "Tal1"
jgene <- "Tead1"
jgene <- "S100a8"
jgene <- "Sema6d"
if (! grepl("\\|", jgene)){
  jrow.tes <- grep(pattern = jgene, rownames(count.mat.norm.tes), value = TRUE)[[1]]
  jrow.tss <- grep(pattern = jgene, rownames(count.mat.norm.tss), value = TRUE)[[1]]
  jcuts.tes <- data.frame(cuts.tes = count.mat.norm.tes[jrow.tes, ], cell = colnames(count.mat.norm.tes), stringsAsFactors = FALSE)
  jcuts.tss <- data.frame(cuts.tss = count.mat.norm.tss[jrow.tss, ], cell = colnames(count.mat.norm.tss), stringsAsFactors = FALSE)
} else {
  jrow.tes <- grep(pattern = jgene, rownames(count.mat.norm.tes), value = TRUE)
  jrow.tss <- grep(pattern = jgene, rownames(count.mat.norm.tss), value = TRUE)
  jcuts.tes <- data.frame(cuts.tes = colSums(count.mat.norm.tes[jrow.tes, ]), cell = colnames(count.mat.norm.tes), stringsAsFactors = FALSE)
  jcuts.tss <- data.frame(cuts.tss = colSums(count.mat.norm.tss[jrow.tss, ]), cell = colnames(count.mat.norm.tss), stringsAsFactors = FALSE)
}

jcuts.merge <- left_join(jcuts.tes, jcuts.tss)

# jmerge <- left_join(left_join(dat.meta, jcuts), jcuts.tss)

jmerge  <- left_join(dat.meta, jcuts.merge)


# ggplot(jmerge, aes(x = umap1, y = umap2, color = log2(cuts.tss + 1))) + 
ggplot(jmerge, aes(x = umap1, y = umap2, color = ifelse(cuts.tss > 0, 1, 0), size = cuts.tss)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(paste("H3K27me3", jgene), jrow.tss) + 
  scale_color_viridis_c(direction = 1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
ggplot(jmerge, aes(x = umap1, y = umap2, color = cluster.fewer)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(paste("H3K27me3")) + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(jmerge, aes(x = cluster, y = log2(cuts.tss))) + 
  geom_boxplot()  + 
  geom_jitter(width = 0.2) + 
  theme_bw() + 
  ggtitle(paste0("H3K27me3", jgene), jrow.tss) + 
  scale_color_viridis_c(direction = -1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# ggplot(jmerge, aes(x = umap1, y = umap2, color = cuts.tes)) + 
#   geom_point() + 
#   theme_bw() + 
#   ggtitle(jgene, jrow.tes) + 
#   scale_color_viridis_c(direction = -1) + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  


# Find highly variables genes  --------------------------------------------

inf.lda.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_TSS/lda_outputs.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TSS.varfilt.K-30.binarize.FALSE/ldaOut.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TSS.varfilt.K-30.Robj"
load(inf.lda.tss, v=T)


library(topicmodels)
library(scchicFuncs)
library(hash)
library(igraph)
library(umap)
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

tm.result <- posterior(out.lda)
tm.result <- AddTopicToTmResult(tm.result)
test <-  sapply(colnames(tm.result$terms), function(x) AssignHash(x = x, jhash = coord2gene, null.fill = NA))
assertthat::assert_that(all(!is.na(test)))
colnames(tm.result$terms) <-  sapply(colnames(tm.result$terms), function(x) AssignHash(x = x, jhash = coord2gene, null.fill = NA))

dat.umap.tss <- DoUmapAndLouvain(tm.result$topics, jsettings)

topics.dat <- OrderTopicsByEntropy(tm.result = tm.result)

outpdf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/topics_loadings_H3K27me3_techrepmerged/TSS_topics_on_peaks_umap.pdf"
pdf(outpdf, useDingbats = FALSE)
for (jtopic in topics.dat$topic){
  print(jtopic)
  loadings <- data.frame(cell = rownames(tm.result$topics), cell.loadings = tm.result$topics[, jtopic], stringsAsFactors = FALSE)
  dat.tmp <- left_join(dat.meta, loadings)
  m <- ggplot(dat.tmp, aes(x = umap1, y = umap2, color = cell.loadings)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jtopic) + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}
dev.off()

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
ggplot(dat.meta, aes(x = umap1, y = umap2, color = cluster.fewer)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Get marker genes --------------------------------------------------------



keeptop <- 400
jannot <- list()
jannot[["topic8"]] <- "Eryths"
jannot[["topic18"]] <- "NKs"
jannot[["topic10"]] <- "ErythsProg"
jannot[["topic1"]] <- "DCs"
jannot[["topic28"]] <- "pDCs"
jannot[["topic15"]] <- "LymphProg"
jannot[["topic9"]] <- "HSPCs"
jannot[["topic22"]] <- "Bcells"
jannot[["topic16"]] <- "HSPCs2"
jannot[["topic7"]] <- "LymphProg2"
jannot[["topic12"]] <- "Granulocytes"
jannot[["topic20"]] <- "Basophils"


outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/topics_loadings_H3K27me3_techrepmerged"

jnames <- names(jannot)
names(jnames) <- jnames

gvec.dat.merged <- lapply(jnames, function(jtopic){
  jclst <- jannot[[jtopic]]
  gvec <- sort(tm.result$terms[jtopic, ], decreasing = TRUE)[1:keeptop]
  ranking <- seq(length(gvec))
  gvec.dat <- data.frame(gene = names(gvec), gene.loading = gvec, topic = jtopic, cluster = jclst, rnk = ranking, stringsAsFactors = FALSE)
}) %>%
  bind_rows() %>%
  arrange(cluster)

outtxt <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/topics_loadings_H3K27me3_techrepmerged/topics_gene_loadings_H3K27me3_techrepmerged.txt"
fwrite(x = gvec.dat.merged, file = outtxt, sep = "\t")



# Load H3K4me1 LDA TSS  ---------------------------------------------------



jmark.check <- "H3K4me3"
inf.annot.k4me1 <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.", jmark.check, ".txt"))
dat.annot.k4me1 <- fread(inf.annot.k4me1)

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
ggplot(dat.annot.k4me1, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point()  + 
  theme_bw() + 
  ggtitle(jmark.check) + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

inf.k4me1 <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_TSS.dist_10000.RemoveDupRows/lda_outputs.count_mat_from_TSS.", jmark.check, ".dist_10000.K-30.binarize.FALSE/ldaOut.count_mat_from_TSS.", jmark.check, ".dist_10000.K-30.Robj"))
load(inf.k4me1, v=T)

count.mat.k4me1 <- count.mat
out.lda.k4me1 <- out.lda

count.mat.k4me1 <- sweep(count.mat.k4me1, MARGIN = 2, STATS = colSums(count.mat.k4me1), FUN = "/")

# take hits and plot UMAP 
# jgene <- "Ppp1r18"
jgene <- "Ppp1r18"

jgene <- "Lrrtm1"
jgene <- "Sox6"
jgene <- "S100a8"
jgene <- "Lrrtm1"
jgene <- "Hbb-y"
jgene <- "Hoxa9"
(jrow <- grep(jgene, rownames(count.mat.k4me1), value = TRUE)[[1]])

jrow <- "17:35861127-35871127;NM_001146710.1..Ppp1r18"

cuts.k4me1 <- data.frame(cell = colnames(count.mat.k4me1), cuts.tss = count.mat.k4me1[jrow, ], stringsAsFactors = FALSE)

jmerged.k4me1 <- left_join(dat.annot.k4me1, cuts.k4me1)

ggplot(jmerged.k4me1, aes(x = umap1, y = umap2, color = ifelse(cuts.tss > 0, 1, 0))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jrow, jmark.check) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jmerged.k4me1, aes(x = cluster, y = log2(cuts.tss))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2) + 
  theme_bw() + 
  ggtitle(jrow) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Sum across sets ---------------------------------------------------------


library(DescTools)
print(unique(gvec.dat.merged$cluster))
jset <- "Granulocytes"
jset <- "Eryths"
jset <- "DCs"
jset <- "pDCs"
jset <- "HSPCs"
jset <- "LymphProg"
jset <- "HSPCs2"
jset <- "Eryths"
jset <- "NKs"
jset <- "Bcells"
jset <- "Granulocytes"
jgenes <- subset(gvec.dat.merged, cluster == jset)$gene

cuts.k4me1.multi <- data.frame(cell = colnames(count.mat.k4me1), cuts.tss = colMeans(count.mat.k4me1[jgenes, ]), stringsAsFactors = FALSE)

jmerged.k4me1.multi <- left_join(dat.annot.k4me1, cuts.k4me1.multi)

ggplot(jmerged.k4me1.multi, aes(x = umap1, y = umap2, color = Winsorize(log2(cuts.tss + 1), probs = c(0, 0.99)))) +
# ggplot(jmerged.k4me1.multi, aes(x = umap1, y = umap2, color = log2(cuts.tss + 1))) +
  geom_point() + 
  theme_bw() + 
  ggtitle("H3K4me1 cuts at TSS defined from H3K27me3 genes:", paste("H3K27me3 gene set:", jset)) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# summarize as a heatmap


clstrs.keep <- sort(unique(dat.annot.k4me1$cluster))
glst <- split(x = gvec.dat.merged$gene, f = gvec.dat.merged$cluster)

dat.ordered <- dat.annot.k4me1 %>% arrange(cluster, jrep)
cells.order <- dat.ordered$cell
clstrs.uniq <- unique(dat.ordered$cluster)
clstrs.new <- dat.ordered$cluster
jcol <- cbPalette[1:length(colsvec.uniq)]
colsvec.hash <- hash::hash(clstrs.uniq, jcol)
colsvec.new <- sapply(clstrs.new, function(x) AssignHash(x, jhash = colsvec.hash, null.fill = NULL))

colsvec.row <- sapply(clstrs.uniq, function(x) AssignHash(x, jhash = colsvec.hash, null.fill = NULL))

pseudogene.mat.k4me1 <- SumAcrossClusters(t(count.mat.k4me1[, cells.order]), cnames.keep.lst = glst[clstrs.keep])
pseudogene.mat.k4me1 <- do.call(rbind, pseudogene.mat.k4me1)

# heatmap3::heatmap3(log2(as.matrix(pseudogene.mat.k4me1) + 1), Rowv = NA, Colv = NA, ColSideColors = colsvec.new, scale = "row", RowSideColors = colsvec.row, revC = TRUE, main = "H3K27me3 new")




outpdf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/topics_loadings_H3K27me3_techrepmerged/check_K27me3_gene_sets.pdf"

pdf(outpdf, useDingbats = FALSE)

ggplot(dat.meta, aes(x = umap1, y = umap2, color = cluster.fewer)) + 
  geom_point() + 
  ggtitle("H3K27me3 UMAP from peaks") + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


heatmap3::heatmap3(log2(as.matrix(pseudogene.mat.k4me1) + 1), Rowv = NA, Colv = NA, ColSideColors = colsvec.new, scale = "row", RowSideColors = colsvec.row, revC = TRUE, main = paste0(jmark.check, " on K27me3-defined genes"))
heatmap3::heatmap3(Winsorize(log2(as.matrix(pseudogene.mat.k4me1) + 1), probs = c(0, 0.94)), Rowv = NA, Colv = NA, ColSideColors = colsvec.new, scale = "row", RowSideColors = colsvec.row, revC = TRUE, main = paste0(jmark.check, " on K27me3-defined genes"))




# Add in Giladi data and check which gene sets are more celltype s --------

m2c <- MarkerToCelltype()
inf.giladi <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData")
load(inf.giladi, v=T)

dat.sum.long <- data.table::melt(dat.sum.norm.quantnorm)
colnames(dat.sum.long) <- c("gene", "celltype", "exprs")

dat.sum.long <- dat.sum.long %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
  filter(!is.nan(zscore))


# take granu
names(glst)
jset <- "HSPCs"
jset <- "pDCs"
jset <- "Eryths"
jset <- "Granulocytes"

for (jset in names(glst)){
  
  genes.keep <- sapply(glst[[jset]], function(x) strsplit(x, "\\.")[[1]][[4]])[1:150]
  
  dat.giladi.filt <- subset(dat.sum.long, gene %in% genes.keep)  %>%
    rowwise() %>%
    mutate(celltype = m2c[[celltype]])
  
  m <- ggplot(dat.giladi.filt, aes(x = forcats::fct_reorder(.f = celltype, .x = zscore, .fun = median, .desc = TRUE), y = zscore)) + 
    geom_boxplot() + 
    geom_jitter(width = 0.2) + 
    ggtitle(paste("Gene set from H3K27me3:", jset)) + 
    theme_bw() +
    xlab("") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
  
  
}



# Compare with a H3K4me3-defined geenset ----------------------------------

inf.gset.k4me3 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-09.from_LDA_topics.condensed.heatmap.famousgenes.keepn_400.refmark_H3K4me3.2020-12-09.txt")
dat.gset.k4me3 <- fread(inf.gset.k4me3)

print(unique(dat.gset.k4me3$jset))

jset2 <- "Eryths"
jset2 <- "HSPCs"
jset2 <- "Basophils"

for (jset2 in unique(dat.gset.k4me3$jset)){
  
  rnames.keep.k4me3 <- subset(dat.gset.k4me3, jset == jset2)$gene[1:150]
  genes.keep.k4me3 <- sapply(rnames.keep.k4me3, function(x) strsplit(x, "\\.")[[1]][[4]])
  
  dat.giladi.filt.k4me3 <- subset(dat.sum.long, gene %in% genes.keep.k4me3) %>%
    rowwise() %>%
    mutate(celltype = m2c[[celltype]])
  
  m2 <- ggplot(dat.giladi.filt.k4me3, aes(x = forcats::fct_reorder(.f = celltype, .x = zscore, .fun = median, .desc = TRUE), y = zscore)) + 
    geom_boxplot() + 
    geom_jitter(width = 0.2) + 
    ggtitle(paste("Gene set from H3K4me3:", jset2)) + 
    theme_bw() +
    xlab("") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m2)
  
  
  
}
dev.off()


