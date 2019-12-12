# Jake Yeung
# Date of Creation: 2019-12-10
# File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/make_tx_dataset_zebrafish_WKM_check_thrombos_naivenorm.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)

library(xlsx)
library(DESeq2)

library(preprocessCore)

# Load  -------------------------------------------------------------------


inf <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/For_Jake/the_massive_complete_zf_dataset.csv.gz"
inf.meta <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/For_Jake/tsne_clusterID_zebrafish_GateID_dataset.csv"
inf.meta2 <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/For_Jake/cluster_info_JY_edited.xlsx"

dat <- as.data.frame(fread(inf, stringsAsFactors = FALSE))
colnames(dat)[[1]] <- "gene"
rownames(dat) <- dat$gene
dat$gene <- NULL

meta <- fread(inf.meta, stringsAsFactors = TRUE)
meta2 <- read.xlsx2(inf.meta2, sheetIndex = 1, header = FALSE, colClasses = c("integer", "character")); colnames(meta2) <- c("ClusterID", "celltype")

colnames(meta) <- c("rname", "V1", "V2", "ClusterID", "experi")

meta <- left_join(meta, meta2)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(meta, aes(x = V1, y = V2, color = celltype)) + geom_point() + scale_color_manual(values = cbPalette) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


meta.summary <- meta %>%
  group_by(celltype) %>%
  summarise(ncells = length(celltype))

# Filter out unwanted celltypes -------------------------------------------

# ctype.exclude <- "thrombocytes"
ctype.exclude <- c()
meta <- subset(meta, !celltype %in% ctype.exclude)
meta2 <- subset(meta2, !celltype %in% ctype.exclude)



# Get sum across celltypes  -----------------------------------------------

ctypes <- as.character(unique(meta$celltype))
names(ctypes) <- ctypes

print(paste("Keeping these celltypes:", unique(ctypes)))

ctype.sum <- lapply(ctypes, function(ctype){
  cells.tmp <- subset(meta, celltype == ctype)$rname
  cols.i <- which(colnames(dat) %in% cells.tmp)
  dat.tmp <- dat[, cols.i]
  genesums <- Matrix::rowSums(dat.tmp) 
  return(genesums)
})
genes <- names(ctype.sum[[1]])
assertthat::assert_that(all(sapply(ctype.sum, function(x) identical(names(x), genes))))

ctype.sum <- ctype.sum %>%
  bind_rows() %>%
  as.data.frame()

rownames(ctype.sum) <- genes

# boxplots
boxplot(log2(ctype.sum))

# normalize using DESeq2
ctype.metadata <- data.frame(celltype = colnames(ctype.sum))
rownames(ctype.metadata) <- colnames(ctype.sum)

ctype.sum.vst <- log2(sweep(ctype.sum, MARGIN = 2, STATS = colSums(ctype.sum), FUN = "/") * 10^6 + 1)

# dds <- DESeqDataSetFromMatrix(countData = round(ctype.sum), colData = ctype.metadata, design = ~1)
# vsd <- vst(dds)
# ctype.sum.vst <- assay(vsd)

# remove bad genes?
bad.genes <- apply(ctype.sum.vst, 1, var) == 0

ctype.sum.vst <- ctype.sum.vst[!bad.genes, ]

pca.out <- prcomp(t(ctype.sum.vst), center = TRUE, scale. = TRUE, retx = TRUE)

plot(pca.out$x[, 1], pca.out$x[, 2], pch = 20)
text(pca.out$x[, 1], pca.out$x[, 2], labels = colnames(ctype.sum.vst))

boxplot(ctype.sum.vst)

 
# Remove lowly expressed genes? -------------------------------------------

xfilt <- 3.5
plot(density(unlist(ctype.sum.vst)))
abline(v = xfilt)

genes.keep.i <- which(apply(ctype.sum.vst, 1, function(x) quantile(x, 0.8)) > xfilt)
# genes.keep.i <- which(apply(ctype.sum.vst, 1, median) > xfilt)

ctype.sum.vst.filt <- ctype.sum.vst[genes.keep.i, ]

# do quantile normalization anyways
dat.mat.filt <- normalize.quantiles(as.matrix(ctype.sum.vst.filt), copy=TRUE)
colnames(dat.mat.filt) <- colnames(ctype.sum.vst.filt); rownames(dat.mat.filt) <- rownames(ctype.sum.vst.filt)

# dont do quantile normalization??
# dat.mat.filt <- ctype.sum.vst.filt

boxplot(dat.mat.filt)

jgenes <- sapply(rownames(dat.mat.filt), function(x) strsplit(x, "_")[[1]][[2]])

jgenes.dup <- jgenes[duplicated(jgenes)]
names(jgenes.dup) <- jgenes.dup

gene.exprs.collapsed <- lapply(jgenes.dup, function(jg){
  # get mean of duplicated genes
  apply(dat.mat.filt[which(jgenes == jg), ], 2, mean)
}) %>%
  do.call(rbind, .)

rows.exclude <- lapply(jgenes.dup, function(jg){
  which(jgenes == jg)
}) %>%
  unlist()

rows.include <- !seq(nrow(dat.mat.filt)) %in% rows.exclude
dat.mat.filt.collapsed <- dat.mat.filt[rows.include, ]

jgenes.collapsed <- sapply(rownames(dat.mat.filt.collapsed), function(x) strsplit(x, "_")[[1]][[2]])

# Get zscores -------------------------------------------------------------

dat.mat.filt.long <- tidyr::gather(data.frame(gene = jgenes.collapsed, dat.mat.filt.collapsed), key = "celltype", value = "exprs", -gene) %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
  ungroup()

dat.mat.filt.zscore <- tidyr::spread(dat.mat.filt.long %>% dplyr::select(-exprs), key = "celltype", value = "zscore") %>%
  as.data.frame()

rownames(dat.mat.filt.zscore) <- dat.mat.filt.zscore$gene

dat.mat.filt.zscore$gene <- NULL

# plot example genes
# jgene <- c("scpp8", "adam8a", "odc1", "lta4h", "thy1", "timp2b")
jgene <- c("rpl38", "med25", "pcdh7b", "pros1", "med25")
jgene <- c("mafba", "dapk1", "lhfpl6", "mrc1b")

jgene <- c("runx1t1", "meis1b", "gata2a", "slc12a4", "znf423", "adarb2", "pax7b", "plxna4", "csf1b", "dnah7", "asphd2", "caskin1", "adgrd1")

ggplot(subset(dat.mat.filt.long, gene == "dnah7"), aes(x = celltype, y = zscore)) + geom_boxplot() + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(dat.mat.filt.long, gene %in% jgene) %>% mutate(celltype = forcats::fct_reorder(celltype, dplyr::desc(zscore))), 
       aes(x = celltype, y = zscore)) +
  geom_boxplot() + 
  geom_point() + 
  ggtitle(jgene) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Save output  ------------------------------------------------------------

outf <- paste0("/Users/yeung/data/scchic/public_data/Zebrafish_WKM/Baron_et_al_pseudobulk_Zebrafish_WKM_cpmnorm.", Sys.Date(), ".rds")
# save(dat.mat.filt.long, ctype.sum.vst, file = outRdata)
# saveRDS(dat.mat.filt.long, file = outf)


