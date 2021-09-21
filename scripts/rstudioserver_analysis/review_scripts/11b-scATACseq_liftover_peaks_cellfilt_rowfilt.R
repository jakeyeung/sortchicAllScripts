# Jake Yeung
# Date of Creation: 2021-08-26
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/11-check_scATACseq_data.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(irlba)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"



# Load  -------------------------------------------------------------------

mapq <- 30
bnames <- c("BoneMarrow_62016", "BoneMarrow_62216")
names(bnames) <- bnames

infs.lst <- lapply(bnames, function(bname){
  # inf.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/public_data/Cusanovich_2018/count_tables/", bname, ".peaks_filt.mapq_", mapq, ".count_table.txt"))
  # inf.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/public_data/Cusanovich_2018/count_tables/", bname, ".10kb_TSS.mapq_", mapq, ".count_table.txt"))
  inf.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/public_data/Cusanovich_2018/count_tables/", bname, ".50kb_TSS.mapq_", mapq, ".count_table.txt"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

dat.mat <- lapply(infs.lst, function(jinf){
  mat.tmp <- fread(jinf)
  mat.sparse <- Matrix(as.matrix(mat.tmp[, -1]), sparse = TRUE)
  rownames(mat.sparse) <- mat.tmp$coordname
  return(mat.sparse)
}) 
dat.mat <- do.call(cbind, dat.mat)

dat.lsi <- scchicFuncs::RunLSI(as.matrix(dat.mat))

dat.umap <- DoUmapAndLouvain(dat.lsi$u, jsettings)

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Get annots --------------------------------------------------------------

inf.meta <- file.path(hubprefix, paste0("jyeung/data/scChiC/public_data/Cusanovich_2018/metadata/cell_metadata.bonemarrow_only.withheader.txt"))
dat.meta <- fread(inf.meta)
dat.umap.tmp <- left_join(dat.umap, subset(dat.meta, select = c(cell, cell_label)))

# Check S100a8  -----------------------------------------------------------

ggplot(dat.umap.tmp, aes(x = umap1, y = umap2, color = cell_label)) +
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jgene <- "Erg$"
jgene <- "Ebf1"
jgene <- "Hbb-"
jgene <- "S100a8$"
jrows.keep <- grep(pattern = jgene, rownames(dat.mat))
print(rownames(dat.mat)[jrows.keep])
if (length(jrows.keep) == 1){
  jcounts <- data.frame(gene_counts = dat.mat[jrows.keep, ], cell = colnames(dat.mat), total_counts = colSums(dat.mat), stringsAsFactors = FALSE)
} else {
  jcounts <- data.frame(gene_counts = colSums(dat.mat[jrows.keep, ]), cell = colnames(dat.mat), total_counts = colSums(dat.mat), stringsAsFactors = FALSE)
}

dat.umap.annot <- left_join(dat.umap.tmp, jcounts)

ggplot(dat.umap.annot %>% arrange(gene_counts), aes(x = umap1, y = umap2, color = log2(gene_counts / total_counts))) +
  geom_point() + 
  theme_bw() + 
  ggtitle(paste("Counts +/-5kb around", jgene)) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.annot %>% mutate(cell_label = cell_label == "Hematopoietic progenitors") %>% arrange(cell_label), aes(x = umap1, y = umap2, color = cell_label)) +
  geom_point() + 
  ggtitle("Highlighting hematopoietic progenitors as TRUE. 10kb bins around TSS") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# check S100a8



# Compare mm9 mat ---------------------------------------------------------

inf.mm9 <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Cusanovich_2018/processed_data/atac_matrix.binary.qc_filtered.bonemarrow_filt.rowsfilt.rds")
dat.mm9 <- readRDS(inf.mm9)

print(dim(dat.mm9))

# Check colsums -----------------------------------------------------------

ncuts.mm10 <- data.frame(cell = colnames(dat.mat), cuts_mm10 = colSums(dat.mat), stringsAsFactors = FALSE)
ncuts.mm9 <- data.frame(cell = colnames(dat.mm9), cuts_mm9 = colSums(dat.mm9), stringsAsFactors = FALSE)

ncuts.merged <- left_join(ncuts.mm9, ncuts.mm10)

ggplot(ncuts.merged, aes(x = cuts_mm9, y = cuts_mm10)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

