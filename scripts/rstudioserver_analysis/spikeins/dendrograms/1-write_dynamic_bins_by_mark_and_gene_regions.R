# Jake Yeung
# Date of Creation: 2021-02-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/1-write_dynamic_bins_by_mark_and_gene_regions.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

hubprefix <- "/home/jyeung/hub_oudenaarden"


indir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned")
outdir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.dynamic_bins_TSS_TES_regions")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
bsize <- 50000

# Load genomic regions  ---------------------------------------------------

de.bins.dat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0("DE_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.", jmark, ".2021-01-30.txt")
  jdat <- fread(file.path(indir, fname))
  jdat <- jdat %>%
    rowwise() %>%
    mutate(startExtend = start + 1 - bsize / 2,
           endExtend = end - 1 + bsize / 2)
  jdat$region_coordExtend <- paste(jdat$seqnames, paste(as.character(jdat$startExtend), as.character(jdat$endExtend), sep = "-"), sep = ":")
  return(jdat)
})


# make beds

de.bins.bed.lst <- lapply(jmarks, function(jmark){
  jcoords <- de.bins.dat.lst[[jmark]]$region_coordExtend
  GetBedFromCoords(jcoords, add.chr = FALSE, strip.chr = TRUE)
})





# Load TSS diff exprs -----------------------------------------------------

inf.mat <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/heatmap_pdfs_and_ordered_matrices/heatmap_ordered_with_labels.H3K4me1.2021-01-08.rearranged.RData")
load(inf.mat, v=T)

rnames <- rownames(mat.adj.tmp)
jgenes <- unique(sapply(rnames, function(x) strsplit(x, "\\.")[[1]][[4]], USE.NAMES = FALSE))
# jgenes <- unique(sapply(rnames, function(x) strsplit(x, ";")[[1]][[2]], USE.NAMES = FALSE))


# Get TES coordinates from gene  ------------------------------------------

inf.tes <- file.path(hubprefix, "jyeung/data/databases/refseq/MmRefseqTssToTes.chromorenamed.merged.rearranged.nochr.bed")
dat.tes <- fread(inf.tes, col.names = c("Chr", "Start", "End", "FullName", "Name")) %>%
  rowwise() %>%
  mutate(Gene = strsplit(Name, split = "\\.")[[1]][[4]],
         Width = End - Start) 

dat.tes.filt <- subset(dat.tes, Gene %in% jgenes)

genes.tes.check <- jgenes[which(!jgenes %in% dat.tes.filt$Gene)]


# Extend only up to 50kb  -------------------------------------------------

maxlength <- 50000

dat.tes.filt.clipped <- dat.tes.filt %>%
  group_by(Gene) %>%
  filter(Width == max(Width))  %>%
  mutate(EndClipped = ifelse(Width > maxlength, Start + maxlength, End),
         WidthClipped = EndClipped - Start) 


# Take +/- 5 kb from TSS  -------------------------------------------------

inf.tss <- file.path(hubprefix, "jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.10000.again.nochromo.sorted.bed")

dat.tss <- fread(inf.tss, col.names = c("Chr", "Start", "End", "FullName", "Name")) %>%
  rowwise() %>%
  mutate(Gene = strsplit(FullName, split = "\\.")[[1]][[4]],
         Width = End - Start)
  # rowwise() %>%
  # mutate(Chr = paste("chr", Chr, sep = ""))

dat.tss.filt <- subset(dat.tss, Gene %in% jgenes)
genes.tss.check <- jgenes[which(!jgenes %in% dat.tss.filt$Gene)]

# Write tables -------------------------------------------------------------

# write these binsjj


for (jmark in jmarks){
  outftmp <- file.path(outdir, paste0("dynamic_bins.50kb.", jmark, ".txt"))
  fwrite(de.bins.bed.lst[[jmark]], file = outftmp, sep = "\t", col.names = FALSE)
}


# write these TES regions
outftes <- file.path(outdir, paste0("celltype_specific_genes.TSS_TES.txt"))
fwrite(dat.tes.filt %>% dplyr::select(Chr, Start, End, FullName), file = outftes, sep = "\t", col.names = FALSE)


# write TSS regions 
outftss <- file.path(outdir, paste0("celltype_specific_genes.TSS_10kb.txt"))
fwrite(dat.tss.filt %>% dplyr::select(Chr, Start, End, FullName), file = outftss, sep = "\t", col.names = FALSE)

