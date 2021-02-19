# Jake Yeung
# Date of Creation: 2021-02-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/4-annotate_6085_signif_bins.R
# See if there's anything interesting in the 6085 bins

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(hash)
library(topicmodels)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

binsize <- 50000

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"

# merge some celltypes
ctypes <- list("Eryths" = "Erythroid", 
               "Bcells" = "Lymphoid", 
               "NKs" = "Lymphoid", 
               "Granulocytes" = "Myeloid",
               "Basophils" = "Myeloid", 
               "pDCs" = "Lymphoid",
               "DCs" = "Myeloid", 
               "HSPCs" = "HSPCs",
               "Erythroid" = "Erythroid",
               "Lymphoid" = "Lymphoid",
               "Myeloid" = "Myeloid")


merge.ctypes.by.lineage <- FALSE

# Load metadata -----------------------------------------------------------

indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage")

dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  fread(file.path(indir.meta, fname)) %>%
    rowwise() %>%
    mutate(lineage = ctypes[[cluster]])
}) 

# cluster to col
cluster2col <- hash(dat.metas$H3K4me1$cluster, dat.metas$H3K4me1$clustercol)
cluster2col[["Erythroid"]] <- "#0072B2"
cluster2col[["Lymphoid"]] <- "#56B4E9"
cluster2col[["Myeloid"]] <- "#D55E00"



# Load LDA  ---------------------------------------------------------------


jmark <- "H3K4me3"
# jsuffixmain <- paste0("dynamic_bins.50kb")
jsuffixmain <- "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table"
# jsuffix <- "celltype_specific_genes.TSS_10kb.txt"
# jsuffix <- "celltype_specific_genes.TSS_TES.txt"
# dat.impute.pbulk.lst <- lapply(jmarks, function(jmark){

# jsuffix.lst <- list("dynamic_bins.50kb", "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table")
# names(jsuffix.lst) <- jsuffix.lst

# jsuffixmain <- "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table.H3K9me3.2021-02-15.txt"


out.lst <- lapply(jmarks, function(jmark){
  # jsuffixmain <- "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table.H3K9me3.2021-02-15.txt"
  
  if (jsuffixmain == "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table"){
    jsuffix <- paste0(jsuffixmain, ".", jmark, ".2021-02-15.txt")
  } else {
    jsuffix <- paste0("dynamic_bins.50kb.", jmark, ".txt")
  }
  
  # jsuffix <- paste0(jsuffixmain, ".", jmark, ".txt")
  print(jmark)
  
  if (jsuffixmain == "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table"){
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM.varfilt.dynamic_bins_genes.top_6085/lda_outputs.count_tables_merged.", jmark, ".", jsuffix, ".K-30.binarize.FALSE/ldaOut.count_tables_merged.", jmark, ".", jsuffix, ".K-30.Robj"))
  } else {
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM.varfilt.dynamic_bins_genes/lda_outputs.count_tables_merged.", jmark, ".", jsuffix, ".K-30.binarize.FALSE/ldaOut.count_tables_merged.", jmark, ".", jsuffix, ".K-30.Robj"))
  }
  assertthat::assert_that(file.exists(inf.lda))
  
  load(inf.lda, v=T)
  
  tm.result <- posterior(out.lda)
  
  dat.impute <- t(tm.result$topics %*% tm.result$terms)
  
  if (merge.ctypes.by.lineage){
    cnames.keep.lst <- split(x = dat.metas[[jmark]]$cell, f = dat.metas[[jmark]]$lineage)
  } else {
    cnames.keep.lst <- split(x = dat.metas[[jmark]]$cell, f = dat.metas[[jmark]]$cluster)
  }
  dat.impute.pbulk <- do.call(cbind, SumAcrossClusters(dat.impute, cnames.keep.lst))
  dat.impute.pbulk <- sweep(dat.impute.pbulk, MARGIN = 2, STATS = colSums(dat.impute.pbulk), FUN = "/", check.margin = TRUE)
  
  dat.raw.pbulk <- SumAcrossClusters(count.mat, cnames.keep.lst)
  dat.raw.pbulk <- do.call(cbind, dat.raw.pbulk)
  # normalize
  dat.raw.pbulk <- sweep(dat.raw.pbulk, MARGIN = 2, STATS = colSums(dat.raw.pbulk), FUN = "/", check.margin = TRUE) * 1000000 + 1
  return(list(dat.impute.pbulk = dat.impute.pbulk, dat.raw.pbulk = dat.raw.pbulk))
})

# jmark.ref <- "H3K9me3"
# dat.raw.pbulk <- out.lst[[jmark.ref]]$dat.raw.pbulk
# plot(density(dat.raw.pbulk))


bins.lst <- lapply(out.lst, function(jout){
  sapply(rownames(jout$dat.raw.pbulk), function(x) strsplit(x, ";")[[1]][[2]], USE.NAMES = FALSE)
})

bed.mid.lst <- lapply(bins.lst, function(jbins){
  jbed <- GetBedFromCoords(coords = jbins, add.chr = FALSE, strip.chr = FALSE) %>%
    rowwise() %>%
    mutate(Midpt = Start + (End - Start) / 2,
           StartNew = Midpt - 1,
           EndNew = Midpt + 1,
           CoordNew = paste(Chr, paste(StartNew, EndNew, sep = "-"), sep = ":"))
})


# Annotate bins  ----------------------------------------------------------

# jinf.tss <- file.path(hubprefix, "jyeung/data/databases/gene_tss/IGtypes/CodingGenes_IGtypes_gene_tss_winsize.50000.bed")
# jinf.tss <- file.path(hubprefix, paste0("jyeung/data/databases/gene_tss/IGtypes/CodingGenes_IGtypes_gene_tss_winsize.100000.bed"))
jinf.tss <- file.path(hubprefix, paste0("jyeung/data/databases/gene_tss/gene_tss_winsize.100000.bed"))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

bins.annot.genewise.lst <- lapply(jmarks, function(jmark){
  bins.annot <- AnnotateCoordsFromList.GeneWise(coords.vec = bed.mid.lst[[jmark]]$CoordNew, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)
  return(bins.annot)
})

# jinf.tss <- file.path(hubprefix, "jyeung/data/databases/gene_tss/IGtypes/CodingGenes_IGtypes_gene_tss_winsize.50000.bed")
bins.annot.binwise.lst <- lapply(jmarks, function(jmark){
  bins.annot <- AnnotateCoordsFromList(coords.vec = bed.mid.lst[[jmark]]$CoordNew, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)
  return(bins.annot)
})


# bins.dist.long <- 

# Check distances ---------------------------------------------------------

bins.dist.long.lst <- lapply(jmarks, function(jmark){
  jdat <- bins.annot.binwise.lst[[jmark]]$regions.annotated
  jdat$mark <- jmark
  
  jdat <- jdat %>%
    rowwise() %>%
    mutate(startExtend = start + 1 - binsize / 2,
           endExtend = end - 1 + binsize / 2)
  jdat$region_coordExtend <- paste(jdat$seqnames, paste(jdat$startExtend, jdat$endExtend, sep = "-"), sep = ":")
  return(jdat)
}) 


m <- ggplot(bins.dist.long.lst %>% bind_rows(), aes(x = log10(abs(distanceToTSS + 1)), fill = mark)) + 
  geom_density(alpha = 0.25) +
  # coord_cartesian(xlim = c(0, 5)) + 
  facet_wrap(~mark, ncol = 1) + 
  geom_vline(xintercept = log10(25000), linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)



# Check GCs  --------------------------------------------------------------

# load GCs 
load("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_gc_analysis/gcs_genomewide.RData", v=T)

gr.gc.dat.filt.lst <- lapply(jmarks, function(jmark){
  jbins <- bins.lst[[jmark]]
  jout <- subset(gr.gc.dat, bname %in% jbins) %>%
    mutate(mark = jmark)
  jout <- jout[!duplicated(jout$bname), ]
}) 

ggplot(gr.gc.dat.filt.lst %>% bind_rows(), aes(x = mark, y = gc)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point() + 
  theme_bw()  + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Check GOs  --------------------------------------------------------------
# do it later


# Write dynamic bins 6085 to output ---------------------------------------

outdir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/tables_top_6085_four_marks_dynamic_bins")
dir.create(outdir)
jmerge.lst <- lapply(jmarks, function(jmark){
  outf <- file.path(outdir, paste0("top_6085_bins_nearest_gene_gc.", jmark, ".", Sys.Date(), ".txt"))
  print(dim(bins.dist.long.lst[[jmark]]))
  jmerge <- left_join(bins.dist.long.lst[[jmark]], subset(gr.gc.dat.filt.lst[[jmark]], select = c(bname, gc)), by = c("region_coordExtend" = "bname"))
  # save to output
  print(dim(jmerge))
  fwrite(x = jmerge, file = outf, quote = FALSE, sep = "\t")
  return(jmerge)
})

 
# write bed file for motevo -----------------------------------------------

bed.lst <- lapply(jmarks, function(jmark){
  outbed <- file.path(outdir, paste0("top_6085.", jmark, ".", Sys.Date(), ".bed"))
  bed.tmp <- GetBedFromCoords(coords = bins.lst[[jmark]], add.chr = FALSE, strip.chr = FALSE)
  fwrite(bed.tmp, outbed, quote = FALSE, sep = "\t", col.names = FALSE)
})



# Write pbulk  ------------------------------------------------------------

dat.impute.pbulk.lst <- lapply(jmarks, function(jmark){
  outimpute <- file.path(outdir, paste0("dat_imputed_pbulk.top_6085.", jmark, ".", Sys.Date(), ".txt"))
  write.table(x = out.lst[[jmark]]$dat.impute.pbulk, file = outimpute, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
  return(out.lst[[jmark]]$dat.impute.pbulk)
})

dat.raw.pbulk.lst <- lapply(jmarks, function(jmark){
  outraw <- file.path(outdir, paste0("dat_raw_pbulk.top_6085.", jmark, ".", Sys.Date(), ".txt"))
  write.table(x = out.lst[[jmark]]$dat.raw.pbulk, file = outraw, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
  return(out.lst[[jmark]]$dat.raw.pbulk)
})



# Find lost bins ----------------------------------------------------------

jthres <- 1
pbulk.diff.lst <- lapply(jmarks, function(jmark){
  jtest <- log2(dat.impute.pbulk.lst[[jmark]]) %>%
    melt() 
  colnames(jtest) <- c("bname", "ctype", "log2signal")
  jtest <- jtest %>%
    rowwise() %>%
    mutate(is.ctype = ifelse(ctype == "HSPCs", "HSPCs", "notHSPCs")) %>%
    group_by(bname, is.ctype) %>%
    summarise(log2signal = mean(log2signal))
  jtest.sum <- jtest %>%
    group_by(bname) %>%
    summarise(log2signal.diff = log2signal[1] - log2signal[2])
  plot(density(jtest.sum$log2signal.diff), main = jmark)
  return(jtest.sum)
  # bins.lost <- sapply(as.character(subset(jtest.sum , log2signal.diff < 0)$bname), function(x) strsplit(x, ";")[[1]][[2]], USE.NAMES = FALSE)
})

bins.lost.lst <- lapply(pbulk.diff.lst, function(jtest.sum){
  bins.lost <- sapply(as.character(subset(jtest.sum , log2signal.diff > jthres)$bname), function(x) strsplit(x, ";")[[1]][[2]], USE.NAMES = FALSE)
})

# write lost bins

lapply(bins.lost.lst, length)

bed.lost.lst <- lapply(jmarks, function(jmark){
  outbed <- file.path(outdir, paste0("lost_bins.", jmark, ".minlog2fc_", jthres, ".", Sys.Date(), ".bed"))
  bed.tmp <- GetBedFromCoords(coords = bins.lost.lst[[jmark]], add.chr = FALSE, strip.chr = FALSE)
  fwrite(bed.tmp, outbed, quote = FALSE, sep = "\t", col.names = FALSE)
  return(bed.tmp)
})



# dat.raw.pbulk.lst <- lapply(jmarks, function(jmark){
#   outraw <- file.path(outdir, paste0("dat_raw_pbulk.top_6085.", jmark, ".txt"))
#   write.table(x = out.lst[[jmark]]$dat.raw.pbulk, file = outraw, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
#   return(out.lst[[jmark]]$dat.raw.pbulk)
# })


# Check Hox cluster -------------------------------------------------------

# dim(subset(bins.dist.long.lst$H3K27me3, grepl("Hox", SYMBOL)))

# subset(bins.dist.long.lst$H3K4me3, grepl("Hox", SYMBOL))



# inf.k9.checks2 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/tables_H3K9me3_dynamic_bins/coords_H3K9me3_dynamic_bins.bed")
# dat.k9.checks2 <- fread(inf.k9.checks2)
# common <- intersect(dat.k9.checks2$V4, bins.lst$H3K9me3)

# 
# # check old hits
# inf.k9.checks <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.H3K9me3.2021-01-30.txt"
# inf.k9.high <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/High_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.H3K9me3.2021-01-30.txt"
# dat.k9.checks <- fread(inf.k9.checks)
# dat.k9.ref <- bins.annot.binwise.lst[["H3K9me3"]]$regions.annotated
# dat.k9.high <- fread(inf.k9.high)
# 
# dat.k9.merge <- rbind(dat.k9.checks %>% mutate(type = "checks"), dat.k9.ref %>% mutate(mark = "H3K9me3", type = "refs")) 
# 
# dat.k9.high <- rbind(dat.k9.merge, dat.k9.high %>% mutate(type = "high"))
# 
# commons <- intersect(dat.k9.ref$region_coord, dat.k9.checks$region_coord)
# 
# ggplot(dat.k9.checks, aes(x = log10(abs(distanceToTSS + 1)))) + 
#   geom_density()  + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(dat.k9.high, aes(x = log10(abs(distanceToTSS + 1)), fill = type)) + 
#   geom_density(alpha = 0.25)  + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# bins.annot.checks <- AnnotateCoordsFromList(coords.vec = dat.k9.checks$region_coord, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)
# 
# ggplot(bins.annot.checks$regions.annotated, aes(x = abs(distanceToTSS + 1))) + 
#   geom_density() + 
#   theme_bw() + 
#   scale_x_log10() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# # Check lost for H3K27me3 and H3K9me3  ------------------------------------
# 
# jmark.ref <- "H3K9me3"
# dat.impute.pbulk <- out.lst[[jmark.ref]]$dat.impute.pbulk
# dat.vars <- sweep(log2(dat.impute.pbulk), MARGIN = 1, STATS = rowMeans(log2(dat.impute.pbulk)), FUN = "-")
# 
# # get HSPC high, and others low
# 
# 
# 
# # Do some GO term analysis ------------------------------------------------
# 
# 
# 

