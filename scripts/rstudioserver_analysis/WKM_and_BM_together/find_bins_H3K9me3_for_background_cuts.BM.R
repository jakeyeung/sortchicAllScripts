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

hubprefix <- "/home/jyeung/hub_oudenaarden"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/offsets_hetero_and_totalcuts"
outf <- file.path(outdir, paste0("MouseBM_HeteroTotalCounts_50kb_bins.", Sys.Date(), ".chrfilt.RData"))
pdfout <- file.path(outdir, paste0("MouseBM_HeteroTotalCounts_50kb_bins.", Sys.Date(), ".chrfilt.pdf"))
bedout <- file.path(outdir, paste0("MouseBM_HeteroTotalCounts_50kb_bins.", Sys.Date(), ".chrfilt.bed"))

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")


pdf(pdfout, useDingbats = FALSE)


# Load cell cluster annots ------------------------------------------------

indir.annots <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping")
dat.annots.all <- lapply(jmarks, function(jmark){
  inf.annots <- file.path(indir.annots, paste0("GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData"))
  load(inf.annots, v=T)
  return(dat.umap.glm.fillNAs)
})


# Cells keep --------------------------------------------------------------

cells.keep <- lapply(dat.annots.all, function(jdat){
  return(jdat$cell)
})


# Load the new 50kb tables ------------------------------------------------

jdist <- "50000"
indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data/debug/ZellerRawData_B6_All_MergedByMarks_final.count_tables_SlidingWindow.slurm.noR2")

mats.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0(jmark, ".mapq_40.SlidingWindow_", jdist, ".blfiltered.csv")
  inf <- file.path(indir, fname)
  mat <- ReadMatSlideWinFormat(inf = inf, as.sparse = TRUE, sort.rnames = TRUE, add.chromo = FALSE)
  mat <- mat[, cells.keep[[jmark]]]
  return(mat)
})

lapply(mats.lst, dim)

# Use K27me3 and K9me3 to find bins that separate between two clus --------

rnames.lst <- lapply(mats.lst, rownames)
rows.common <- Reduce(f = intersect, x = rnames.lst)
k27me3.exprs <- data.frame(bin = rows.common, K27me3.cuts = rowSums(mats.lst$H3K27me3[rows.common, ]), stringsAsFactors = FALSE)
k9me3.exprs <- data.frame(bin = rows.common, K9me3.cuts = rowSums(mats.lst$H3K9me3[rows.common, ]), stringsAsFactors = FALSE)

merged.exprs <- left_join(k27me3.exprs, k9me3.exprs)

max.counts <- 2000
merged.exprs.highcut <- merged.exprs %>% filter(K27me3.cuts < max.counts & K9me3.cuts < max.counts)
ggplot(merged.exprs.highcut,
       aes(x = K27me3.cuts, y = K9me3.cuts)) + geom_point(alpha = 0.1) +  
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  geom_density_2d() + 
  geom_hline(yintercept = 1000, color = 'red') + 
  geom_vline(xintercept = 155, color = 'red') + 
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

# filter out chromosomes
k9.chrs <- sapply(k9.bins, function(x) GetChromo(x))

k9.bins.keep <- which(k9.chrs %in% jchromos)

print(paste("Nbins before:", length(k9.bins)))
k9.bins <- k9.bins[k9.bins.keep]
print(paste("Nbins after:", length(k9.bins)))


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
    rowwise() %>%
    mutate(cond = GetCondFromSamp(cell, mark = jmark))
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

dev.off()

save(dat.ncuts.hetero.total, k9.bins, rows.common, cells.keep, file = outf)


# Wriite bins to output ---------------------------------------------------

# as bed

datbed <- data.frame(chromo = sapply(k9.bins, GetChromo), 
                     Start = sapply(k9.bins, GetStart), 
                     End = sapply(k9.bins, GetEnd),
                     stringsAsFactors = FALSE)

fwrite(datbed, file = bedout, sep = "\t", col.names = FALSE, na="NA", quote = FALSE)




