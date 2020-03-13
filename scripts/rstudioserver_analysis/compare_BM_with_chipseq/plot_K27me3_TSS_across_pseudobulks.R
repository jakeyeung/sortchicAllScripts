# Jake Yeung
# Date of Creation: 2020-03-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/compare_BM_with_chipseq/plot_K27me3_TSS_across_pseudobulks.R
# Pseudoublks 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(reticulate)

library(DescTools)

library(purrr)
library(ggrastr)
library(hash)

reticulate::source_python("scripts/python_functions/parse_dictionary_text.py")


# Get gene annots ---------------------------------------------------------

inf.annots <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData"
load(inf.annots, v=T)

dat.sum.long <- data.table::melt(dat.sum.norm.quantnorm)
colnames(dat.sum.long) <- c("gene", "celltype", "exprs")

dat.sum.long <- dat.sum.long %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
  filter(!is.nan(zscore))



# Get gene from coord -----------------------------------------------------

inf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/giladi_filtered.allctypes/gene_tss_winsize.2.BM_giladi.keeptop_1000.bed"
dat.tss <- fread(inf.tss)
tss.hash <- hash(sapply(dat.tss$V4, function(x) strsplit(x, ";")[[1]][[1]]), 
                 paste(sapply(dat.tss$V4, function(x) strsplit(x, ";")[[1]][[2]])))

#x <- dat.tss$V4[[1]]
#strsplit(x, split = ";")[[1]][[2]]

# Load data  --------------------------------------------------------------

inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/bigwig_outputs/merged_bams.deeptools_outputs.tss.MAPQ_40.dist_5000.allctypes.H3K27me3.bsize_100/computeMatrix.MAPQ_40.H3K27me3.allctypes.genetss.sorted.tab.gz"
meta <- inf2dic(inf)
mat.all <- fread(inf, header = FALSE, skip = 1, sep = "\t", quote = "")


# name columns 
nbins <- (unique(meta$downstream) + unique(meta$upstream)) / unique(meta[["bin size"]])
nsamps <- length(meta$sample_labels)

# Stack matrix vertically  ------------------------------------------------

mat <- mat.all[, -seq(6)]
coords <- mat.all[, c(seq(6))]
colnames(coords) <- c("chromo", "start", "end", "coord", "meta1", "meta2")
assertthat::assert_that(ncol(mat) == nsamps * nbins)

# split into a list of matrices, then rbind
sampids <- ceiling(seq(ncol(mat)) / nbins)
sampids.uniq <- as.character(unique(sort(sampids)))
names(sampids.uniq) <- meta$sample_labels

colids <- ( (seq(ncol(mat)) - 1) %% nbins ) + 1  # - 1 and +1 so first element is 1, last element is 10

colnames(mat) <- as.character(sampids)

mats.lst <- lapply(sampids.uniq, function(sampid){
  cols.i <- which(colnames(mat) == sampid)
  mat.sub <- mat[, ..cols.i]
  colnames(mat.sub) <- paste("bin", seq(nbins), sep = "")
  return(mat.sub) 
})

# make long dataframe of bed locations
coords.lst <- lapply(sampids.uniq, function(sampid){
  jtmp <- coords[, c(1,2,3,4)]
  jtmp$sampid <- sampid
  jtmp$gene <- sapply(jtmp$coord, function(x) tss.hash[[x]])
  return(jtmp)
})  

mats.long <- do.call(rbind, mats.lst)
coords.long <- do.call(rbind, coords.lst)

mat.long.merge <- cbind(mats.long, coords.long)

quantile(unlist(mats.lst[[1]]), 0.99)

jmat <- as.matrix(mats.lst[[1]])
print(range(jmat))
jmat.win <- Winsorize(jmat, minval = 0.01, maxval = quantile(jmat, 0.98))
print(range(jmat.win))

# Plot heatmap for top neutrophil genes? ----------------------------------

x <- subset(dat.sum.long, grepl("S100a8", gene))
keeptop <- 100
jgenes <- subset(dat.sum.long, celltype == "Ltf") %>%
  ungroup() %>%
  top_n(n = 100, wt = zscore)


