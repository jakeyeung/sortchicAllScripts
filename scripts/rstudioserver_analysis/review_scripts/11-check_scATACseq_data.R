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

hubprefix <- "/home/jyeung/hub_oudenaarden"


# functions --------------------------------------------------------------

FilterRows <- function(mat, site_frequency_threshold = 0.03){
  num_cells_ncounted = Matrix::rowSums(mat)
  threshold = ncol(mat) * site_frequency_threshold
  mat.filt = mat[num_cells_ncounted >= threshold,]
  return(mat.filt)
}



# Add meta  ---------------------------------------------------------------

inf.meta <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Cusanovich_2018/metadata/cell_metadata.bonemarrow_only.withheader.txt")
dat.meta <- fread(inf.meta)
dat.meta.filt <- subset(dat.meta, select = c(cell, cell_label, tissue.replicate))
dat.meta.filt$barcode <- paste(dat.meta.filt$tissue.replicate, dat.meta.filt$cell, sep = ".")

# ggplot(dat.meta, aes(x = tsne_1, y = tsne_2, color = cell_label)) + 
ggplot(dat.meta, aes(x = subset_tsne1, y = subset_tsne2, color = cell_label)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load  -------------------------------------------------------------------

indir <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Cusanovich_2018/count_tables")
bnames <- c("BoneMarrow_62016", "BoneMarrow_62216")
names(bnames) <- bnames

infs <- lapply(bnames, function(bname){
  # inf <- file.path(indir, paste0(bname, ".bam.remapped.sorted.count_table.txt"))
  inf <- file.path(indir, paste0(bname, ".bam.remapped.sorted.TSS_50kb.mapq_40.count_table.txt"))
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

counts.lst <- lapply(bnames, function(bname){
  inf <- infs[[bname]]
  mat <- as.data.frame(data.table::fread(inf))
  rownames(mat) <- mat$coordname
  mat$coordname <- NULL
  colnames(mat) <- paste(bname, colnames(mat), sep = ".")
  return(mat)
})

cuts.total.long <- lapply(bnames, function(bname){
  jmat <- counts.lst[[bname]]
  dat <- data.frame(ncuts = colSums(jmat), cell = colnames(jmat), stringsAsFactors = FALSE)
  return(dat)
}) %>%
  bind_rows()

# filter 
cuts.total.filt.long <- subset(cuts.total.long, cell %in% dat.meta.filt$barcode)

plot(density(log(cuts.total.filt.long$ncuts)))




# Filter cols -------------------------------------------------------------

counts.colsfilt <- lapply(counts.lst, function(jmat){
  cols.keep <- colnames(jmat) %in% dat.meta.filt$barcode
  jmat[, cols.keep]
}) 
counts.colsfilt <- do.call(cbind, counts.colsfilt)

# Filter rows  ------------------------------------------------------------

counts.colsfilt.rowsfilt <- FilterRows(mat = counts.colsfilt, site_frequency_threshold = 0.01)


# Run LSI  ----------------------------------------------------------------

library(irlba)
library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

lsi.out <- scchicFuncs::RunLSI(as.matrix(counts.colsfilt.rowsfilt))

dat.umap <- DoUmapAndLouvain(lsi.out$u, jsettings)

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) +  
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

dat.umap.merge <- left_join(dat.umap, subset(dat.meta.filt, select = c(cell_label, barcode)), by =c("cell" = "barcode"))

ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cell_label)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Get rds  ----------------------------------------------------------------


inf.rds.filt <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Cusanovich_2018/processed_data/atac_matrix.binary.qc_filtered.bonemarrow_filt.rds")
dat.rds.filt <- readRDS(inf.rds.filt)

# save rows of mm9
outdir <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Cusanovich_2018/processed_data")
outf.rows <- file.path(outdir, paste0("atac_matrix.bonemarrow_filt.rownames.bed"))
rnames <- rownames(dat.rds.filt)

dat.bed <- data.frame(Chromo = sapply(rnames, function(x) strsplit(x, "_")[[1]][[1]]), 
                      Start = sapply(rnames, function(x) strsplit(x, "_")[[1]][[2]]),
                      End = sapply(rnames, function(x) strsplit(x, "_")[[1]][[3]]),
                      stringsAsFactors = FALSE)
data.table::fwrite(dat.bed, file = outf.rows, sep = "\t", col.names = FALSE)


plot(density(as.numeric(dat.bed$End) - as.numeric(dat.bed$Start)))

# inf.rds <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Cusanovich_2018/processed_data/atac_matrix.binary.qc_filtered.rds")
# dat.rds.full <- readRDS(inf.rds)

# outf.rds.filt <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Cusanovich_2018/processed_data/atac_matrix.binary.qc_filtered.bonemarrow_filt.rds")
# dat.rds <- readRDS(inf.rds)
# 
# cells.keep <- dat.meta$cell
# cols.keep <- colnames(dat.rds) %in% cells.keep
# dat.rds.filt <- dat.rds[, cols.keep]
# 
# saveRDS(dat.rds.filt, file = outf.rds.filt)




# Run  --------------------------------------------------------------------


# filter peaks 
site_frequency_threshold <- 0.03
num_cells_ncounted = Matrix::rowSums(dat.rds.filt)
threshold = ncol(dat.rds.filt) * site_frequency_threshold
dat.rds.filt.filt = dat.rds.filt[num_cells_ncounted >= threshold,]

print(dim(dat.rds.filt.filt))

lsi.out.peaks <- scchicFuncs::RunLSI(as.matrix(dat.rds.filt.filt))

dat.umap.peaks <- DoUmapAndLouvain(lsi.out.peaks$u, jsettings = jsettings)

ggplot(dat.umap.peaks %>% left_join(., dat.meta.filt) %>% filter(cell_label != "Unknown"), aes(x = umap1, y = umap2, color = cell_label)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("scATAC-seq mouse bone marrow of Cusanovich 2018") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# write filt filt
outf.rds.filt.filt <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Cusanovich_2018/processed_data/atac_matrix.binary.qc_filtered.bonemarrow_filt.rowsfilt.rds")
saveRDS(dat.rds.filt.filt, file = outf.rds.filt.filt)

outf.rows.filt <- file.path(outdir, paste0("atac_matrix.bonemarrow_filt.rownames_filt.mm9.bed"))
rnames.filt <- rownames(dat.rds.filt.filt)

dat.bed.filt <- data.frame(Chromo = sapply(rnames.filt, function(x) strsplit(x, "_")[[1]][[1]]), 
                      Start = sapply(rnames.filt, function(x) strsplit(x, "_")[[1]][[2]]),
                      End = sapply(rnames.filt, function(x) strsplit(x, "_")[[1]][[3]]),
                      stringsAsFactors = FALSE)
data.table::fwrite(dat.bed.filt, file = outf.rows.filt, sep = "\t", col.names = FALSE)

