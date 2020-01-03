# Jake Yeung
# Date of Creation: 2019-12-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/intestinal_all_merged/1-setup_pseudobulk_scrnaseq.R
# Pseudobulk analysis

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(DESeq2)

# Load data ---------------------------------------------------------------

inf <- "/home/jyeung/hpc/intestinal_scchic/public_data/GSE92332_Regional_UMIcounts.txt.gz"
assertthat::assert_that(file.exists(inf))

dat <- fread(inf)
# make matrix
dat.mat <- as.matrix(subset(dat, select = -V1))
rownames(dat.mat) <- dat$V1

# remove cells
cellsums <- colSums(dat.mat)
genesums <- rowSums(dat.mat)

plot(density(log10(cellsums)))
plot(density(log10(genesums)))

genes.filt <- which(genesums > 10)

dat.mat.filt <- dat.mat[genes.filt, ]

# remove bad genes
meta <- data.frame(samp = colnames(dat.mat), celltype = sapply(colnames(dat.mat), function(x) strsplit(x, "_")[[1]][[4]]))

# sum across celltypes
# Get sum across celltypes  -----------------------------------------------

ctypes <- as.character(unique(meta$celltype))
names(ctypes) <- ctypes

ctype.sum <- lapply(ctypes, function(ctype){
  cells.tmp <- subset(meta, celltype == ctype)$samp
  cols.i <- which(colnames(dat.mat.filt) %in% cells.tmp)
  dat.tmp <- dat.mat.filt[, cols.i]
  genesums <- Matrix::rowSums(dat.tmp)
  return(genesums)
})

genes <- names(ctype.sum[[1]])
assertthat::assert_that(all(sapply(ctype.sum, function(x) identical(names(x), genes))))

ctype.sum <- ctype.sum %>%
  bind_rows() %>%
  as.data.frame()
rownames(ctype.sum) <- genes

boxplot(log2(ctype.sum), ylab = "log2 UMIs")

# normalize using DESeq2
ctype.metadata <- data.frame(celltype = colnames(ctype.sum))
rownames(ctype.metadata) <- colnames(ctype.sum)

dds <- DESeqDataSetFromMatrix(countData = ctype.sum, colData = ctype.metadata, design = ~1)

vsd <- vst(dds)

ctype.sum.vst <- assay(vsd)

pca.out <- prcomp(t(ctype.sum.vst), center = TRUE, scale. = TRUE, retx = TRUE)

plot(pca.out$x[, 1], pca.out$x[, 2], pch = 20)
text(pca.out$x[, 1], pca.out$x[, 2], labels = colnames(ctype.sum.vst))

boxplot(ctype.sum.vst)

# Quant Norm not necessary? 


# Save to output?  --------------------------------------------------------

outf <- "/home/jyeung/hpc/intestinal_scchic/public_data/pseudobulk_data/pseudobulk_Haber_et_al_intestinal_celltypes_2019-12-30.VST_noQuantNorm.rds"
outpdf <- "/home/jyeung/hpc/intestinal_scchic/public_data/pseudobulk_data/pseudobulk_Haber_et_al_intestinal_celltypes_2019-12-30.pdf"


pdf(outpdf, useDingbats = FALSE)
  boxplot(log2(ctype.sum), ylab = "log2 UMIs")
  plot(pca.out$x[, 1], pca.out$x[, 2], pch = 20)
  text(pca.out$x[, 1], pca.out$x[, 2], labels = colnames(ctype.sum.vst))
  boxplot(ctype.sum.vst)
dev.off()

saveRDS(ctype.sum.vst, outf)





