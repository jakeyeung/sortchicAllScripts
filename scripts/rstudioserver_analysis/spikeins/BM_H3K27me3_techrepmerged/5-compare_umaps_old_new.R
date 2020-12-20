# Jake Yeung
# Date of Creation: 2020-12-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_techrepmerged/5-compare_umaps_old_new.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(hash)

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

# Load old umap  ----------------------------------------------------------

inf.old <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.H3K27me3.txt"
inf.new <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/metadata_umap_celltype_cuts.H3K27me3.txt"

dat.old <- fread(inf.old)
dat.new <- fread(inf.new)

clstr2hash <- hash(dat.new$cluster, dat.new$colorcode)

dat.old$colorcode <- sapply(dat.old$cluster, function(x) AssignHash(x = x, jhash = clstr2hash, null.fill = x))

ggplot(dat.old, aes(x = umap1, y = umap2, color = colorcode)) + 
  geom_point() + 
  facet_wrap(~jrep) + 
  theme_bw() + 
  scale_color_identity() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.new, aes(x = umap1, y = umap2, color = colorcode)) + 
  geom_point() + 
  facet_wrap(~jrep) + 
  theme_bw() + 
  geom_vline(xintercept = 2.2, linetype = "dotted") + 
  geom_hline(yintercept = 1, linetype = "dotted") + 
  scale_color_identity() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# # why some yellow cells? 
# jsub <- subset(dat.new, umap1 > 2.2 & umap2 < 1)
# table(jsub$cluster)


# Load old/new raw counts TSS -------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

infnewtss <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_TSS/lda_outputs.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TSS.varfilt.K-30.binarize.FALSE/ldaOut.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TSS.varfilt.K-30.Robj")
load(infnewtss, v=T)
count.mat.newtss <- count.mat

infoldtss <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_TSS.dist_10000.RemoveDupRows/lda_outputs.count_mat_from_TSS.H3K27me3.dist_10000.K-30.binarize.FALSE/ldaOut.count_mat_from_TSS.H3K27me3.dist_10000.K-30.Robj")
load(infoldtss, v=T)
count.mat.oldtss <- count.mat

# add old batch into matrix
cells.to.add <- subset(dat.old, jrep == "rep1old")$cell

mat.to.add <- count.mat.oldtss[, cells.to.add]

rnames <- rownames(count.mat.oldtss)
coords <- paste("chr", sapply(rownames(count.mat.oldtss), function(x) strsplit(x, ";")[[1]][[1]]), sep = "")

coords2rnames <- hash::hash(coords, rnames)
count.mat.newtss.rename <- count.mat.newtss
rownames(count.mat.newtss.rename) <- sapply(rownames(count.mat.newtss.rename), function(x) AssignHash(x = x, jhash = coords2rnames, null.fill = x))

bins.to.keep <- intersect(rownames(count.mat.newtss.rename), rownames(count.mat.oldtss))

mat.tss.cbind <- cbind(count.mat.newtss.rename[bins.to.keep, ], mat.to.add[bins.to.keep, ])

print(dim(count.mat.newtss))
print(dim(count.mat.oldtss))


# Load old/new raw counts peaks -------------------------------------------------


infnewpeak <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt_peaks/lda_outputs.PZ-BM-rep2rep3reseq-H3K27me3.peaks.varfilt.K-30.binarize.FALSE/ldaOut.PZ-BM-rep2rep3reseq-H3K27me3.peaks.varfilt.K-30.Robj")
load(infnewpeak, v=T)
count.mat.newpeak <- count.mat

infoldpeak <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins.from_sitecount_mat.from_same_annot_file/lda_outputs.count_mat_from_sitecount_mat.H3K27me3.filtNAcells_allbins.from_same_annot_file.K-30.binarize.FALSE/ldaOut.count_mat_from_sitecount_mat.H3K27me3.filtNAcells_allbins.from_same_annot_file.K-30.Robj")
load(infoldpeak, v=T)
count.mat.oldpeak <- count.mat


mat.to.add.peak <- count.mat.oldpeak[, cells.to.add]

rnames.peak <- rownames(count.mat.oldpeak)
coords.peak <- paste("chr", sapply(rownames(count.mat.oldpeak), function(x) strsplit(x, ";")[[1]][[1]]), sep = "")

coords2rnames.peak <- hash::hash(coords.peak, rnames.peak)
count.mat.newpeak.rename <- count.mat.newpeak
rownames(count.mat.newpeak.rename) <- sapply(rownames(count.mat.newpeak.rename), function(x) AssignHash(x = x, jhash = coords2rnames.peak, null.fill = x))

peak.to.keep <- intersect(rownames(count.mat.newpeak.rename), rownames(count.mat.oldpeak))

mat.peak.cbind <- cbind(count.mat.newpeak.rename[peak.to.keep, ], mat.to.add.peak[peak.to.keep, ])


# Write new metadata ------------------------------------------------------


dat.old.filt <- subset(dat.old, cell %in% cells.to.add)

colnames(dat.old.filt)
colnames(dat.new)

# combine metas old and new

# # why is there NAs and empty marks? 
# inf.metanew <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq/varfilt_2/BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.metadata.txt"
# dat.metanew <- fread(inf.metanew)

dat.new.rename <- subset(dat.new, select = c(cell, umap1, umap2, louvain, cluster, plate, batch, cuts_in_peak, cuts_total, spikein_cuts, rowcoord, colcoord, jrep, mark, colorcode))
# dat.new.rename <- subset(dat.new, select = c(cell, umap1, umap2, louvain, cluster, plate, batch, rowcoord, colcoord, jrep, colorcode))

# dat.new.rename <- left_join(dat.new.rename, subset(dat.metanew, select = c(cell, cuts_in_peak, cuts_total, spikein_cuts, mark)))

# update metadata cuts_in_peak, cuts_total, spikein_cuts, mark

dat.new.rename.withold <- rbind(dat.new.rename, dat.old.filt) %>%
  arrange(cluster, jrep)

# # combine with old 
# dat.new.withold <- rbind(dat.new.rename, dat.old.filt)
# dat.new.withold$mark <- "H3K27me3"


# Save to output  ---------------------------------------------------------

# save meta
outf.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq.with_old/mat_H3K27me3_rep1rep2rep3reseq.metadata.txt"
fwrite(x = dat.new.rename.withold, file = outf.meta)

outf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq.with_old/mat_H3K27me3_rep1rep2rep3reseq.TSS.rds"
outf.peaks <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq.with_old/mat_H3K27me3_rep1rep2rep3reseq.peaks.rds"

if (!file.exists(outf.tss)){
  saveRDS(object = mat.tss.cbind, file = outf.tss)
} 
if (!file.exists(outf.peaks)){
  saveRDS(object = mat.peak.cbind, file = outf.peaks)
}



