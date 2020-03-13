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

jmark <- "H3K27me3"
do.zscore <- FALSE

# Get gene annots ---------------------------------------------------------

inf.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat.rds"
dat.annot <- readRDS(inf.annot)

# Get gene from coord -----------------------------------------------------

inf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/giladi_filtered.seurat/gene_tss_winsize.2.diff_exprs_Giladi_seurat.bed"
dat.tss <- fread(inf.tss)
tss.hash <- hash(sapply(dat.tss$V4, function(x) strsplit(x, ";")[[1]][[1]]), 
                 paste(sapply(dat.tss$V4, function(x) strsplit(x, ";")[[1]][[2]])))

#x <- dat.tss$V4[[1]]
#strsplit(x, split = ";")[[1]][[2]]

# Load data  --------------------------------------------------------------

inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/bigwig_outputs/merged_bams.deeptools_outputs.tss.MAPQ_40.dist_10000.allctypes_from_seurat.", jmark, ".bsize_100/computeMatrix.MAPQ_40.", jmark, ".gene_tss_winsize.2.diff_exprs_Giladi_seurat.tab.gz")
assertthat::assert_that(file.exists(inf))

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


# Plot gene levels  -------------------------------------------------------


mats.lst.clean <- lapply(mats.lst, function(jmat){
  jmat <- as.matrix(jmat)
  jmat.win <- Winsorize(jmat, minval = 0, maxval = quantile(jmat, 0.98))
  return(jmat.win)
}) 

gene.exprs <- lapply(mats.lst.clean, function(jmat){
  rowMeans(jmat)
})

cnames <- names(gene.exprs)[2:length(gene.exprs)]
names(cnames) <- cnames


# Make matrix -------------------------------------------------------------

exprs.mat <- do.call(cbind, gene.exprs)
rownames(exprs.mat) <- paste(coords.lst[[1]]$coord, coords.lst[[1]]$gene, sep = ";")


# Plot top neutrophil genes? ----------------------------------------------

jclusts <- unique(dat.annot$clust)

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_genes_analysis_seurat_TSS"
outpdf <- file.path(outdir, paste0(jmark, "_top_DE_genes_TSS_signal_across_pseudobulks.zscore_", do.zscore, ".pdf"))

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
pdf(outpdf, useDingbats = FALSE)
for (jclust in jclusts){
  # jclust <- "core"
  print(jclust)
  genes.keep <- subset(dat.annot, cluster == jclust & p_val_adj < 0.01 & avg_logFC > 0.5)$gene
  print(length(genes.keep))
  
  # match to coords
  jsub  <- subset(coords.lst[[1]], gene %in% genes.keep)
  coords.keep <- paste(jsub$coord, jsub$gene, sep = ";")
  
  # sort coords by exprs in neutrophil then plot? 
  exprs.sub <- exprs.mat[coords.keep, ]
  samps.keep <- which(colSums(exprs.sub) > 0)
  exprs.sub <- exprs.sub[, samps.keep]
  
  exprs.sub.long <- tidyr::gather(data.frame(coord = rownames(exprs.sub), exprs.sub, stringsAsFactors = FALSE), "pseudobulk", "cpm", -coord)
  
  if (do.zscore){
    exprs.sub.long <- exprs.sub.long %>%
      group_by(coord) %>%
      mutate(cpm = scale(cpm, center = TRUE, scale = TRUE))
  }
  m <- ggplot(exprs.sub.long, aes(x = cpm, group = pseudobulk, fill = pseudobulk)) + geom_density(alpha = 0.5) + 
    theme_bw() + theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_manual(values = cbPalette) + facet_wrap(~pseudobulk, ncol = 2) + ggtitle(paste("Top DE genes (n=", length(genes.keep), ") for ", jclust))
  print(m)
}
dev.off()


# ggplot(exprs.sub.long %>% filter(!grepl("HSC|Linneg", key)), aes(x = value, group = key, fill = key)) + geom_density(alpha = 0.5) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
#   scale_color_manual(values = cbPalette)  + facet_wrap(~key) + ggtitle(jclust)



