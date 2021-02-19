# Jake Yeung
# Date of Creation: 2021-01-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/15-setup_dynamic_bins.from_heatmap.R
# Take k9 dynamic bins from heatmap (fewer bins)


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"



# Load dynamic bins in heatmap (low and high) ---------------------------------------------------------------

jlow.in.k9 <- TRUE
jkeeptop <- 150

jfits.lst.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jinf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.RData"))
  load(jinf, v=T)
  return(jfits.lst)
})

params.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfits.lst <- jfits.lst.lst[[jmark]]
  jnames <- names(jfits.lst); names(jnames) <- jnames
  params.dat.all <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xkeep <- grepl("^Cluster.*.Estimate$", x = names(x))
    jparams <- x[xkeep]
    data.frame(bin = jname, param = names(jparams), estimate = unlist(jparams), mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  if (jmark == "H3K9me3"){
    params.dat.all <- params.dat.all %>%
      mutate(param = gsub("Eryth", "Eryths", param),
             param = gsub("Lymphoid", "Bcells", param))
  }
  # make params more readable
  params.dat.all$ctype <- params.dat.all$param
  params.dat.all$ctype <- gsub("Cluster", "", params.dat.all$ctype)
  params.dat.all$ctype <- gsub(".Estimate", "", params.dat.all$ctype)
  return(params.dat.all)
})

pvals.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfits.lst <- jfits.lst.lst[[jmark]]
  jnames <- names(jfits.lst); names(jnames) <- jnames
  lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
})

pvals.lst2 <- lapply(jfits.lst.lst$H3K9me3, function(x) x$pval)
k9.bins <- which(pvals.lst2 < 1e-10)


pval.k9.sub <- subset(pvals.lst$H3K9me3, pval < 1e-10) %>%
  arrange(desc(pval))

k9.bins.names <- unique(pval.k9.sub$bin)
ctypes.keep <- c("Eryths", "Bcells", "Granulocytes")
params.keep <- paste("Cluster", ctypes.keep, ".Estimate", sep = "")



params.dat.wide.lst <- lapply(jmarks, function(jmark){
  jsub <- subset(params.lst[[jmark]], bin %in% k9.bins.names & param %in% params.keep) %>%
    group_by(bin) %>% filter(max(abs(estimate)) < 5)
  jdat <- GetParamsWideFormat(jsub, jvalue.var = "estimate")
  # # keep only effect cnames ( do this later maybe ? )
  # cnames.keep.i <- grep("effect$", colnames(jdat))
  # cnames.new <- paste(colnames(jdat)[cnames.keep.i], jmark, sep = "_")
  # colnames(jdat)[cnames.keep.i] <- cnames.new
  # cnames.keep.bin.i <- grep("bin", colnames(jdat))
  # cnames.keep.merged.i <- c(cnames.keep.bin.i, cnames.keep.i)
  # jdat.filt <- jdat[, cnames.keep.merged.i]
  return(jdat)
})

bins.keep.low.lst <- GetK9CelltypeBins(params.dat.wide.lst$H3K9me3, low.in.k9 = TRUE, keeptop = jkeeptop)
bins.keep.high.lst <- GetK9CelltypeBins(params.dat.wide.lst$H3K9me3, low.in.k9 = FALSE, keeptop = jkeeptop)
bnames <- names(bins.keep.low.lst); names(bnames) <- bnames

jbins.eryth <- c(bins.keep.low.lst[["Eryths"]], bins.keep.high.lst[["Eryths"]])
jbins.bcell <- c(bins.keep.low.lst[["Bcells"]], bins.keep.high.lst[["Bcells"]])
jbins.granu <- c(bins.keep.low.lst[["Granulocytes"]], bins.keep.high.lst[["Granulocytes"]])
jbins.hspcs <- c(bins.keep.low.lst[["HSPCs"]], bins.keep.high.lst[["HSPCs"]])

# jbins.bcell <- bins.keep.lst[["Bcells"]]
# jbins.granu <- bins.keep.lst[["Granulocytes"]]
# jbins.hspcs <- bins.keep.lst[["HSPCs"]]

bins.keep <- c(jbins.eryth, jbins.bcell, jbins.granu, jbins.hspcs)




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




# Create a K9me3 - K4me1 region ilst --------------------------------------

dat.k9.regions <- GetBedFromCoords(bins.keep, add.chr = FALSE, strip.chr = TRUE)
dat.k9.regions <- dat.k9.regions[!duplicated(dat.k9.regions), ]

dat.k4.regions <- dat.tes.filt.clipped %>%
  ungroup() %>%
  dplyr::select(Chr, Start, EndClipped, Name) %>%
  dplyr::rename(End = EndClipped)



# Write outputs -----------------------------------------------------------

dat.combined.regions <- rbind(dat.k4.regions, dat.k9.regions)
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/regions_H3K4me1_H3K9me3_dynamic_regions"
outfile <- file.path(outdir, paste0("H3K4me1_H3K9me3_celltype_specific_genes_and_bins_from_heatmap.", Sys.Date(), ".txt"))
fwrite(dat.combined.regions, file = outfile, sep = "\t", col.names = FALSE)




