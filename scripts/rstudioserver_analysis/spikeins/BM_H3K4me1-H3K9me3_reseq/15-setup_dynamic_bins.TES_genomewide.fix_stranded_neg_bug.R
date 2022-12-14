# Jake Yeung
# Date of Creation: 2021-06-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/15-setup_dynamic_bins.TES_genomewide.fix_stranded_neg_bug.R
# Negative strand needs to subtract from TES 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

# 
# # Load fits ---------------------------------------------------------------
# 
# jmark.ref <- "H3K9me3"
# 
# outs.lst <- lapply(jmarks, function(jmark){
#   inf.fits <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData")
#   load(inf.fits, v=T)
#   params.long <- SummarizeParamsPvalues(jfits.lst, jmark = jmark, paramname = "Cluster") %>%
#     mutate(log2fc = estimate / log(2))
#   params.long$padj <- p.adjust(params.long$pval.param)
#   jnames <- names(jfits.lst); names(jnames) <- jnames
#   pvals.long <- lapply(jnames, function(jname){
#     x <- jfits.lst[[jname]]
#     xvec <- x$pval
#     data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
#   }) %>%
#     bind_rows()
#   return(list(params.long = params.long, pvals.long = pvals.long))
# })
# 
# pvals.long.lst <- lapply(outs.lst, function(jout) jout$pvals.long)
# params.long.lst <- lapply(outs.lst, function(jout) jout$params.long)
# 
# 
# 
# # Get dyanmics bins -------------------------------------------------------
# 
# pval.cutoff1 <- 1e-50
# pval.cutoff2 <- 1e-10
# 
# jmark1 <- "H3K4me1"
# jmark2 <- "H3K9me3"
# 
# bins.filt1 <- subset(pvals.long.lst[[jmark1]], pval < pval.cutoff1)$bin
# bins.filt2 <- subset(pvals.long.lst[[jmark2]], pval < pval.cutoff2)$bin


# # Load TSS diff exprs -----------------------------------------------------
# 
# inf.mat <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/heatmap_pdfs_and_ordered_matrices/heatmap_ordered_with_labels.H3K4me1.2021-01-08.rearranged.RData")
# load(inf.mat, v=T)
# 
# rnames <- rownames(mat.adj.tmp)
# jgenes <- unique(sapply(rnames, function(x) strsplit(x, "\\.")[[1]][[4]], USE.NAMES = FALSE))


# Get TES coordinates from gene  ------------------------------------------

inf.tes <- file.path(hubprefix, "jyeung/data/databases/refseq/MmRefseqTssToTes.chromorenamed.merged.rearranged.nochr.bed")
dat.tes <- fread(inf.tes, col.names = c("Chr", "Start", "End", "FullName", "Name")) %>%
  rowwise() %>%
  mutate(Gene = strsplit(Name, split = "\\.")[[1]][[4]],
         Width = End - Start, 
         Strand = substr(FullName, nchar(FullName), nchar(FullName)))
# jgenes <- unique(dat.tes$Gene)
# dat.tes.filt <- subset(dat.tes, Gene %in% jgenes)
dat.tes.filt <- dat.tes

print(head(dat.tes.filt))

# genes.tes.check <- jgenes[which(!jgenes %in% dat.tes.filt$Gene)]


# Extend only up to 50kb  -------------------------------------------------

maxlength <- 50000

dat.tes.filt.clipped.pos <- dat.tes.filt %>%
  group_by(Gene) %>%
  filter(Width == max(Width) & Strand == "+")  %>%
  mutate(StartClipped = Start, 
         EndClipped = ifelse(Width > maxlength, Start + maxlength, End),
         WidthClipped = EndClipped - StartClipped)

dat.tes.filt.clipped.neg <- dat.tes.filt %>%
  group_by(Gene) %>%
  filter(Width == max(Width) & Strand == "-")  %>%
  mutate(StartClipped = ifelse(Width > maxlength, End - maxlength, Start), 
         EndClipped = End,
         WidthClipped = EndClipped - StartClipped)

dat.tes.filt.clipped <- rbind(dat.tes.filt.clipped.pos, dat.tes.filt.clipped.neg)


# Create a K9me3 - K4me1 region ilst --------------------------------------

GetBedFromCoords <- function(coords, add.chr = FALSE, strip.chr = FALSE){
  dat.bed <- data.frame(Chr = sapply(coords, JFuncs::GetChromo),
                        Start = sapply(coords, JFuncs::GetStart, returnAsInt = TRUE), 
                        End = sapply(coords, JFuncs::GetEnd, returnAsInt = TRUE),
                        Name = coords, 
                        stringsAsFactors = FALSE)
  if (strip.chr){
    dat.bed$Chr <- gsub("^chr", "", dat.bed$Chr)
  }
  return(dat.bed)
}

dat.tes.regions <- dat.tes.filt.clipped %>%
  ungroup() %>%
  dplyr::select(Chr, StartClipped, EndClipped, Name, Strand) %>%
  dplyr::rename(End = EndClipped,
                Start = StartClipped) %>%
  mutate(FullName = paste(Name, Strand, sep = "_")) %>%
  dplyr::select(Chr, Start, End, FullName)



# Write outputs -----------------------------------------------------------

dat.combined.regions <- dat.tes.regions
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/regions_H3K4me1_H3K9me3_dynamic_regions"
outfile <- file.path(outdir, paste0("TES_genomewide.neg_strand_bug_fixed.", Sys.Date(), ".txt"))
fwrite(dat.combined.regions, file = outfile, sep = "\t", col.names = FALSE)




