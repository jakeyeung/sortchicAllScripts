# Jake Yeung
# Date of Creation: 2021-03-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/8-break_variance_by_celltype_vs_other_regions.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



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

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"


outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps/variance_breakdown_by_genespec_bins_hspcs.", Sys.Date(), ".pdf")
pdf(outpdf, useDingbats = FALSE)

# Load metadata -----------------------------------------------------------

indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage")


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

cname2color <- hash::hash(dat.metas$H3K4me1$cluster, dat.metas$H3K4me1$clustercol)
colpalette <- colorRampPalette(c("grey1", "grey35", "grey99"))(1024)



# Funmctions --------------------------------------------------------------

inf.pbulks <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pseudobulk_high_bins_matrix/dat_raw_four_marks_high_bins_pseudobulk.2021-03-04.rds"
pbulks.lst <- readRDS(inf.pbulks)


# Assign each bin to either celltype region or not ------------------------

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
jinf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed"

coords.lst <- lapply(pbulks.lst, function(x){
  rownames(x)
  # sapply(rownames(x$dat.raw.pbulk), function(x) strsplit(x, ";")[[1]][[2]], USE.NAMES = FALSE)
})

coords.annot.lst <- lapply(coords.lst, function(coords){
  AnnotateCoordsFromList.GeneWise(coords.vec = coords, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)
})



# celtype specific genes  -------------------------------------------------------------

jmarks.sub <- c("H3K4me1", "H3K9me3")
names(jmarks.sub) <- jmarks.sub

indir.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM_dbl_reseq.varfilt.k4_k9_dynamic_bins"
outs.lst <- lapply(jmarks.sub, function(jmark){
  print(jmark)
  inf.tmp <- file.path(indir.lda, paste0("lda_outputs.count_name.", jmark, ".k4_k9_dynamic_bins.2021-01-30.K-30.binarize.FALSE/ldaOut.count_name.", jmark, ".k4_k9_dynamic_bins.2021-01-30.K-30.Robj"))
  load(inf.tmp, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})

dat.raw.lst <- lapply(outs.lst, function(x){
  x$count.mat
})

rnames.both <- rownames(dat.raw.lst$H3K4me1)
rnames.ctypespec <- rnames.both[grepl(";NM", rnames.both)]

gene.names <- sapply(rnames.ctypespec, function(x) strsplit(x, "\\.")[[1]][[4]])

length(gene.names)

bins.filt.lst <- lapply(coords.annot.lst, function(coords.annot){
  subset(coords.annot$out2.df, gene %in% gene.names)$region_coord
})




# Calculate total variance for each bin  ----------------------------------


log2fc.mat.lst <- lapply(pbulks.lst, function(pbulk){
  jdat.log2 <- log2(pbulk)
  jcheck <- sweep(jdat.log2, MARGIN = 1, STATS = rowMeans(jdat.log2), FUN = "-") 
  return(jcheck)
})

var.mat.lst <- lapply(log2fc.mat.lst, function(jcheck){
  jcheck <- jcheck ^ 2
  return(jcheck)
})

var.total.lst <- lapply(var.mat.lst, function(varmat){
  jsums <- rowSums(varmat)
  # names(jsums) <- sapply(names(jsums), function(x) strsplit(x, ";")[[1]][[2]])
  jsums.dat <- data.frame(bin = names(jsums), genevar = jsums, stringsAsFactors = FALSE)
  return(jsums.dat)
})

var.total.annot.lst <- lapply(jmarks, function(jmark){
  bins.in.ctype <- bins.filt.lst[[jmark]]
  var.sub <- var.total.lst[[jmark]]
  var.sub$is.ctypespec <- sapply(var.sub$bin, function(x) x %in% bins.in.ctype)
  return(var.sub)
})

var.sum.annot.long <- lapply(jmarks, function(jmark){
  x <- var.total.annot.lst[[jmark]]
  xsum <- x %>%
    group_by(is.ctypespec) %>%
    summarise(genevar.sum = sum(genevar))
  xsum$mark <- jmark
  return(xsum)
}) %>%
  bind_rows() %>%
  group_by(mark) %>%
  mutate(genevar.norm = genevar.sum / sum(genevar.sum))

ggplot(var.sum.annot.long, aes(x = mark, y = genevar.norm, fill = is.ctypespec)) + 
  geom_col() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jfilt <- var.sum.annot.long %>% filter(is.ctypespec)

ggplot(jfilt, aes(x = forcats::fct_reorder(.f = mark, .x = genevar.norm, .desc = TRUE), y = genevar.norm)) + 
  geom_col() + 
  theme_bw() + 
  xlab("") + 
  ylab("Fraction of Total Variance at Celltype Specific Markers") + 
  ggtitle("Fraction of total variance at celltype specific markers \n (i.e. genes in Fig 2)") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())




# Go by TSS ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
inf.tss <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/heatmap_pdfs_and_ordered_matrices/heatmap_ordered_with_labels.H3K4me1.2021-01-08.rearranged.RData")
load(inf.tss, v=T)

rnames.tss <- sapply(paste("chr", rownames(mat.adj.tmp), sep = ""), function(x) strsplit(x, ";")[[1]][[1]], USE.NAMES = FALSE)
rnames.tss.orig <- paste("chr", rownames(mat.adj.tmp), sep = "")

# sapply(rnames.tss, function(x) strsplit(x, ";")[[1]][[1]])

# find overlaps for each set of bins

regions.tss <- data.frame(seqnames = sapply(rnames.tss, GetChromo),
                      start = sapply(rnames.tss, GetStart),
                      end = sapply(rnames.tss, GetEnd),
                      stringsAsFactors = FALSE)
rownames(regions.tss) <- rnames.tss.orig

regions.bins.lst <- lapply(coords.lst, function(coords.vec){
  regions <- data.frame(seqnames = sapply(coords.vec, GetChromo),
                        start = sapply(coords.vec, GetStart),
                        end = sapply(coords.vec, GetEnd),
                        stringsAsFactors = FALSE)
})

overlaps.dat.lst <- lapply(jmarks, function(jmark){
  annots.gr <- makeGRangesFromDataFrame(regions.bins.lst[[jmark]], keep.extra.columns = TRUE)
  annots.tss.gr <- makeGRangesFromDataFrame(regions.tss, keep.extra.columns = TRUE)
  out2 <- findOverlaps(query = annots.tss.gr, subject = annots.gr, type = "any")
  bins.overlap <- names(annots.gr[subjectHits(out2),])
  return(bins.overlap)
})

query.dat.lst <- lapply(jmarks, function(jmark){
  annots.gr <- makeGRangesFromDataFrame(regions.bins.lst[[jmark]], keep.extra.columns = TRUE)
  annots.tss.gr <- makeGRangesFromDataFrame(regions.tss, keep.extra.columns = TRUE)
  out2 <- findOverlaps(query = annots.tss.gr, subject = annots.gr, type = "any")
  out2.df = data.frame(bin = names(annots.gr[subjectHits(out2),]), annots.gr[subjectHits(out2),], tss = names(annots.tss.gr[queryHits(out2),]), annots.tss.gr[queryHits(out2),])
  # rownames(out2.df) <- out2.df$tss
  return(out2.df)
})


var.total.annot.overlap.lst <- lapply(jmarks, function(jmark){
  bins.in.ctype <- overlaps.dat.lst[[jmark]]
  var.sub <- var.total.lst[[jmark]]
  var.sub$is.ctypespec <- sapply(var.sub$bin, function(x) x %in% bins.in.ctype)
  return(var.sub)
})


var.sum.annot.overlap.long <- lapply(jmarks, function(jmark){
  x <- var.total.annot.overlap.lst[[jmark]]
  xsum <- x %>%
    group_by(is.ctypespec) %>%
    summarise(genevar.sum = sum(genevar))
  xsum$mark <- jmark
  return(xsum)
}) %>%
  bind_rows() %>%
  group_by(mark) %>%
  mutate(genevar.norm = genevar.sum / sum(genevar.sum))

jfilt.overlap <- var.sum.annot.overlap.long %>% filter(is.ctypespec)
ggplot(jfilt.overlap, aes(x = forcats::fct_reorder(.f = mark, .x = genevar.norm, .desc = TRUE), y = genevar.norm)) + 
  geom_col() + 
  theme_bw() + 
  xlab("") + 
  ylab("Fraction of Total Variance at Celltype Specific Markers") + 
  ggtitle("Fraction of total variance at celltype specific markers \n (i.e. genes in Fig 2)") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Get bins ssociated with each TSS  ---------------------------------------

# jmark.test <- "H3K9me3"
jmark.test <- "H3K4me1"
jmark.test <- "H3K27me3"
bins.filt <- overlaps.dat.lst[[jmark.test]]
jmethod <- "complete"

# jmat2 <- log2(pbulks.lst[[jmark.test]][bins.filt, c(-6, -8)])
jmat2 <- log2(pbulks.lst[[jmark.test]][bins.filt, ])
print(dim(jmat2))
Ngenes <- nrow(jmat2)
cnames.color <- sapply(colnames(jmat2), function(x) AssignHash(x, cname2color, null.fill = x))

heatmap3::heatmap3(jmat2, 
                   Rowv = NA, Colv = TRUE, scale = "row",  revC = TRUE, 
                   main = paste("peaks", jmark.test, jmethod, Ngenes), margins = c(5, 8), 
                   cexRow = 0.5, method = jmethod, ColSideColors = cnames.color, col = colpalette)




# Estimate log2fc at celltype specific regions all marks  -----------------

log2fc.long.lst <- lapply(jmarks, function(jmark){
  jmat <- log2fc.mat.lst[[jmark]]
  jlong <- jmat %>%
    melt()
  colnames(jlong) <- c("bin", "celltype", "log2fc")
  bins.in.ctype <- overlaps.dat.lst[[jmark]]
  jlong$is.ctype.spec <- sapply(jlong$bin, function(x) x %in% bins.in.ctype)
  jlong$mark <- jmark
  return(jlong)
})

ggplot(log2fc.long.lst$H3K27me3, aes(x = log2fc, fill = is.ctype.spec)) + 
  geom_density(alpha = 0.25) + 
  facet_wrap(~celltype) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(log2fc.long.lst %>% bind_rows() %>% filter(celltype == "HSPCs" & mark != "H3K9me3"), aes(x = log2fc, fill = is.ctype.spec)) + 
  geom_density(alpha = 0.25) + 
  facet_wrap(~mark) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(log2fc.long.lst %>% bind_rows() %>% filter(celltype == "HSPCs" & mark != "H3K9me3" & is.ctype.spec), aes(x = log2fc, fill = mark)) + 
  geom_density(alpha = 0.25) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw() +  
  xlab("log2fc from global mean") + 
  ggtitle("is.ctype.spec, HSPC") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(log2fc.long.lst %>% bind_rows() %>% filter(celltype == "HSPCs" & mark != "H3K9me3" & !is.ctype.spec), aes(x = log2fc, fill = mark)) + 
  geom_density(alpha = 0.25) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  xlab("log2fc from global mean") + 
  ggtitle("is.not.ctype.spec, HSPC") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dev.off()