# Jake Yeung
# Date of Creation: 2020-06-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/find_bins_H3K9me3_for_background_cuts.R
# Use H3K9me3 bins for background cuts for Poisson model 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

species <- "ZebrafishWKM"

hubprefix <- "/home/jyeung/hub_oudenaarden"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/offsets_hetero_and_totalcuts"
outf <- file.path(outdir, paste0(species, "_HeteroTotalCounts_50kb_bins.", Sys.Date(), ".RData"))
pdfout <- file.path(outdir, paste0(species, "_HeteroTotalCounts_50kb_bins.", Sys.Date(), ".pdf"))

# assertthat::assert_that(!file.exists(outf))
# assertthat::assert_that(!file.exists(pdfout))

make.plots <- TRUE
if (make.plots){
  pdf(pdfout, useDingbats = FALSE)
}

# Load cell cluster annots ------------------------------------------------


wkm.rename <- hash(c("eryth1", "eryth2", "HSC1", "HSC2", "monocyte"), c("eryth", "eryth", "HSPCs", "HSPCs", "granu"))

dat.annot.lst.WKM <- lapply(jmarks, function(jmark){
  print(jmark)
  # filter by previously defined celltypes? 
  inf.annot.louv <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/LDA_downstream/LDA_downstream_ZF.2020-04-23.imputevarfilt.lessstringent/ZF_LDA_output.", jmark, ".keepn_150.final.ClusterTables.txt")
  assertthat::assert_that(file.exists(inf.annot.louv))
  
  inf.annot.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt")
  assertthat::assert_that(file.exists(inf.annot.glmpca))
  
  annot.louv <- fread(inf.annot.louv)
  annot.louv$clusterplate <- paste(annot.louv$cluster, annot.louv$plate, "_")
  
  annot.glmpca <- fread(inf.annot.glmpca)
  annot.glmpca.filt <- subset(annot.glmpca, cell %in% annot.louv$cell) %>%
    rowwise() %>%
    mutate(clusterplate = paste(cluster, plate, sep = "_")) %>%
    mutate(cluster = ifelse(cluster %in% c("lymph1", "lymph2"), "lymph", cluster)) %>%   # rename lymph1 and lymph2 into lymph
    ungroup() %>%
    filter(cluster != "Unknown") %>%
    left_join(., subset(annot.louv, select = c(cell, var.imputed)))
  # rename clusters
  annot.glmpca.filt$cluster <- sapply(annot.glmpca.filt$cluster, function(jclst) AssignHash(jclst, wkm.rename, null.fill = jclst))
  annot.glmpca.filt$cond <- sapply(annot.glmpca.filt$cell, ClipLast, jsep = "_")
  print("annot glmpca filt")
  print(annot.glmpca.filt)
  return(annot.glmpca.filt)
})


# Cells keep --------------------------------------------------------------

cells.keep <- lapply(dat.annot.lst.WKM, function(jdat){
  return(jdat$cell)
})


# Load the new 50kb tables ------------------------------------------------

jdist <- "50000"
# indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data/debug/ZellerRawData_B6_All_MergedByMarks_final.count_tables_SlidingWindow.slurm.noR2")
indir <- file.path(hubprefix, paste0("jyeung/data/zebrafish_scchic/count_tables.winsize_", jdist, ".imputevarfilt.lessstringent.mapq_40"))

mats.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0(jmark, ".imputevarfilt.lessstringent.mapq_40.remerged.countTable.csv")
  inf <- file.path(indir, fname)
  mat <- ReadMatSlideWinFormat(inf = inf, as.sparse = TRUE, sort.rnames = TRUE, add.chromo = FALSE)
  mat <- mat[, cells.keep[[jmark]]]
  return(mat)
})

lapply(mats.lst, dim)

# Use K27me3 and K9me3 to find bins that separate between two clus --------

rnames.lst <- lapply(mats.lst, rownames)
rows.common <- Reduce(f = intersect, x = rnames.lst)
k27me3.exprs <- data.frame(bin = rows.common, K27me3.cuts = rowMeans(mats.lst$H3K27me3[rows.common, ]), stringsAsFactors = FALSE)
k9me3.exprs <- data.frame(bin = rows.common, K9me3.cuts = rowMeans(mats.lst$H3K9me3[rows.common, ]), stringsAsFactors = FALSE)

merged.exprs <- left_join(k27me3.exprs, k9me3.exprs)

max.counts <- 1.25
merged.exprs.highcut <- merged.exprs %>% filter(K27me3.cuts < max.counts & K9me3.cuts < max.counts)
ggplot(merged.exprs.highcut,
       aes(x = K27me3.cuts, y = K9me3.cuts)) + geom_point(alpha = 0.1) +  
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  ggtitle(paste("K9 vs K27 relationship K9 reads less than", max.counts, "removed"))

  
# take top 90% K9me3, bottom 10%
jlow <- 0.25
jhigh <- 0.75

merged.exprs.filt <- merged.exprs.highcut %>%
  ungroup() %>%
  filter(K27me3.cuts <= quantile(K27me3.cuts, jlow),
         K9me3.cuts >= quantile(K9me3.cuts, jhigh))

range(merged.exprs.filt$K27me3.cuts)
range(merged.exprs.filt$K9me3.cuts)

ggplot(merged.exprs.filt, 
       aes(x = K27me3.cuts, y = K9me3.cuts)) + geom_point(alpha = 0.1) +  
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste("Bins kept after filtering K27me3<", jlow, "K9me3", jhigh))

k9.bins <- merged.exprs.filt$bin

ggplot(merged.exprs.highcut %>% mutate(is.heterochromatin = bin %in% k9.bins), 
       aes(x = K27me3.cuts, y = K9me3.cuts, color = is.heterochromatin)) + geom_point(alpha = 0.1) +  
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  geom_density_2d(color = 'grey45') + 
  ggtitle("K9 vs K27 relationship")

# Count cuts in background bins  ------------------------------------------

dat.cuts.hetero <- lapply(jmarks, function(jmark){
  jmat <- mats.lst[[jmark]][k9.bins, ]
  jcounts <- colSums(jmat)
  jdat <- data.frame(cell = names(jcounts), ncuts.hetero = jcounts, stringsAsFactors = FALSE)
  return(jdat)
})


# Explore background over total  ------------------------------------------

dat.ncuts.hetero.total <- lapply(jmarks, function(jmark){
  print(jmark)
  jmat <- mats.lst[[jmark]]
  data.frame(cell = colnames(jmat), ncuts.total = colSums(jmat), stringsAsFactors = FALSE) %>%
    left_join(dat.cuts.hetero[[jmark]])  %>%
    left_join(subset(dat.annot.lst.WKM[[jmark]], select = c(cell, cond)))
}) 

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

mlst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.ncuts.hetero.total[[jmark]], aes(x = ncuts.hetero / ncuts.total, fill = cond)) + geom_density(alpha = 0.25) + ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values = cbPalette)
})
print(mlst)
mlst.log <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.ncuts.hetero.total[[jmark]], aes(x = ncuts.hetero / ncuts.total, fill = cond)) + geom_density(alpha = 0.25) + ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values = cbPalette) + scale_x_log10() 
})
print(mlst.log)

mlst.total <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.ncuts.hetero.total[[jmark]], aes(x = ncuts.total, fill = cond)) + geom_density(alpha = 0.25) + ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values = cbPalette)
})
print(mlst.total)
mlst.total.log <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.ncuts.hetero.total[[jmark]], aes(x = ncuts.total, fill = cond)) + geom_density(alpha = 0.25) + ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values = cbPalette) + scale_x_log10() 
})
print(mlst.total.log)

mlst.hetero <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.ncuts.hetero.total[[jmark]], aes(x = ncuts.hetero, fill = cond)) + geom_density(alpha = 0.25) + ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values = cbPalette)
})
print(mlst.hetero)
mlst.total.hetero.log <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.ncuts.hetero.total[[jmark]], aes(x = ncuts.hetero, fill = cond)) + geom_density(alpha = 0.25) + ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values = cbPalette) + scale_x_log10() 
})
print(mlst.total.hetero.log)


# Write output ------------------------------------------------------------

if (make.plots){
  dev.off()
}

save(dat.ncuts.hetero.total, k9.bins, rows.common, cells.keep, file = outf)

