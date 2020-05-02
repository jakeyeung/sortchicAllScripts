# Jake Yeung
# Date of Creation: 2020-05-02
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/8-finalize_downsampled_TSS_counts_integrated_all_3_marks.R
# 


rm(list=ls())

jstart <- Sys.time()

set.seed(0)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(purrr)

library(hash)
library(igraph)
library(umap)

library(JFuncs)
library(scchicFuncs)

library(topicmodels)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)
library(GenomicRanges)


library(ggrastr)

# Constants ---------------------------------------------------------------

jwinsize <- "10000"


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
jvarcutoffs <- c(0.75, 2, 1, 0.5)
names(jvarcutoffs) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

jchromos <- paste("chr", seq(25), sep = "")

jsize <- 4
cbPalette <- c("grey85", "#32CD32", "#56B4E9", "#CC79A7", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/pseudobulk_analysis.integrated_analysis.TSS"
dir.create(outdir)
outname <- paste0("integrated_analysis.", Sys.Date(), ".pdf")

outpdf <- file.path(outdir, outname)

pdf(outpdf, useDingbats = FALSE)


# Read TSS Signal to figure out which transcript to keep  -----------------

indir.tss <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.TSS.CodingOnly.winsize_", jwinsize)
assertthat::assert_that(dir.exists(indir.tss))

tss.out <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.tss <- file.path(indir.tss, paste0("PZ-ChIC-ZF_", jmark, "_2020-04-07.countTable.TSS.csv"))
  mat.tss <- ReadMatTSSFormat(inf.tss)
  mat.tss <- CollapseRowsByGene(mat.tss, as.long=FALSE, track.kept.gene = TRUE)
  return(list(mat.tss = mat.tss, tss.exprs = rowSums(mat.tss)))
})

tss.exprs.lst <- lapply(tss.out, function(x) x$tss.exprs)
tss.mats.singlecell <- lapply(tss.out, function(x) x$mat.tss)



# exprs.vec <- tss.exprs.lst$H3K4me1

lapply(jmarks, function(jmark){
  plot(density(tss.exprs.lst[[jmark]]), main = jmark)
})


# go with the K4me3 definition...

ref.mark <- "H3K4me3"
jthres <- 50  # maybe not exactly at hump? what about tissuespecific stuff? rare celltypes? complicated from the bulk 
plot(density(tss.exprs.lst[[ref.mark]]))
abline(v = jthres)

# # pick TSS's that are active, the rest are thrown out 
# tss.keep.i <- which(tss.exprs.lst[[ref.mark]] >= jthres)
# tss.keep <- tss.exprs.lst[[ref.mark]][tss.keep.i]



# Get common rows ---------------------------------------------------------

lapply(tss.exprs.lst, length)

tss.all <- lapply(tss.exprs.lst, function(exprs.lst){
  names(exprs.lst)
}) %>%
  unlist() %>%
  unique()

tss.common <- lapply(tss.exprs.lst, function(exprs.lst){
  names(exprs.lst) 
}) %>%
  Reduce(f = intersect, .)

# get ensembl names ? 
genes.common <- sapply(tss.common, function(x) strsplit(x, ";")[[1]][[2]])
ens.common <- Gene2Ensembl.ZF(genes.common, return.original = TRUE, species = "drerio")

# create tss, genes, ens dat
genes.annot <- data.frame(bin = tss.common, gene = genes.common, ens = ens.common, stringsAsFactors = FALSE)

# Load UMAP and GLM annots ------------------------------------------------

annot.glmpca.filt.lst <- lapply(jmarks, function(jmark){
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
    filter(cluster != "Unknown")  # remove the small cluster Unknown
  annot.glmpca.filt <- left_join(annot.glmpca.filt, subset(annot.louv, select = c(cell, var.imputed)))
  return(annot.glmpca.filt)
})



# Make pseudobulks --------------------------------------------------------


cells.keep.lst <- lapply(annot.glmpca.filt.lst, function(x){
  return(x$cell)
})

# filter for good cells and good bins
tss.mats.singlecell.filt <- lapply(jmarks, function(jmark){
  jmat <- tss.mats.singlecell[[jmark]]
  cells.keep <- cells.keep.lst[[jmark]]
  jmat.filt <- jmat[tss.common, cells.keep]
  return(jmat.filt)
})



# Pseudobulks  ------------------------------------------------------------

cnames.keep.lst.lst <- lapply(jmarks, function(jmark){
  lapply(split(annot.glmpca.filt.lst[[jmark]], annot.glmpca.filt.lst[[jmark]]$cluster), function(x){
    x$cell
  })
})

# try to downsample to lowest number of cells? 
lapply(cnames.keep.lst.lst, function(cnames.by.mark) lapply(cnames.by.mark, function(cnames.by.cluster) length(cnames.by.cluster)))


# down sample pbulk to 100 cells
min.cells.keep <- 100
cnames.keep.lst.lst <- lapply(cnames.keep.lst.lst, function(cnames.by.mark){
  lapply(cnames.by.mark, function(cnames.by.cluster){
    return(sample(cnames.by.cluster, size = min.cells.keep, replace = FALSE))
  })
})

lapply(cnames.keep.lst.lst, function(cnames.by.mark) lapply(cnames.by.mark, function(cnames.by.cluster) length(cnames.by.cluster)))

tss.mats.pbulk <- lapply(jmarks, function(jmark){
  pbulk <- SumAcrossClusters(tss.mats.singlecell.filt[[jmark]], cnames.keep.lst.lst[[jmark]])
  return(do.call(cbind, pbulk))
})


# Make long  --------------------------------------------------------------

jlong.lst <- lapply(jmarks, function(jmark){
  jmat.tmp <- tss.mats.pbulk[[jmark]]
  dat.tss <- data.frame(bin = rownames(jmat.tmp), jmat.tmp, stringsAsFactors = FALSE) %>%
    reshape2::melt(., id.var = "bin", variable.name = "cluster", value.name = "counts")
})


# Plot K4me3 vs K27me3  ---------------------------------------------------

jname.hash <- hash::hash(c("eryth1", "HSC1"), c("eryth", "HSC"))
jnames.keep <- c("eryth", "HSC", "lymph", "monocyte")


jmerged <- lapply(jmarks, function(jmark){
  print(jmark)
  jtmp <- jlong.lst[[jmark]] %>%
    rowwise() %>%
    mutate(cluster = as.character(cluster), 
           cluster = ifelse(is.null(jname.hash[[cluster]]), cluster, jname.hash[[cluster]])) %>%
    ungroup() %>%
    filter(cluster %in% jnames.keep)
  colnames(jtmp) <- c("bin", "cluster", jmark)
  return(jtmp)
}) %>%
  purrr::reduce(., left_join)

print(head(jmerged))



# color genes by genesets

ggplot(jmerged, aes(x = H3K4me3, y = H3K27me3, color = cluster)) + 
  geom_point_rast(alpha = 0.25, size = jsize) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster) + 
  geom_density_2d(color = "black")

ggplot(jmerged, aes(x = H3K4me1, y = H3K27me3, color = cluster)) + 
  geom_point_rast(alpha = 0.25, size = jsize) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster) + 
  geom_density_2d(color = "black")

ggplot(jmerged, aes(x = H3K4me1, y = H3K4me3, color = cluster)) + 
  geom_point_rast(alpha = 0.25, size = jsize) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster) + 
  geom_density_2d(color = "black")


# Define sets of genes from literature ------------------------------------

g2e <- hash::hash(genes.common, ens.common)

genesets <- list()

inf.kob.erythmyeloid <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Kobayashi_2019_zf_genelists/erythmyeloid_genes_list.txt"
kobayashi.genes.erythmyeloid <- fread(inf.kob.erythmyeloid, header = FALSE)$V1
jens.erythmyeloid <- kobayashi.genes.erythmyeloid

inf.kob <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Kobayashi_2019_zf_genelists/hspc_genes_list.txt"
kobayashi.genes <- fread(inf.kob, header = FALSE)$V1
jens.hspc <- kobayashi.genes
# HSPCs
jgenes.choose <- c("meis1b", "pmp22b", "ahnak", "krt8", "anxa5b", "mrc1b", "pa2g4a", "npm1a", "ahcy", "adh5", "fabp3", "myb")
# add kobayashi?
jgenes.choose <- c(jgenes.choose, kobayashi.genes)
jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)

jens.choose <- jens.choose[which(!is.na(jens.choose))]

genesets[["HSC"]] <- jens.choose

# lymphs
jgenes.choose <- c("pax5", "cd79a", "bhlhe40", "cd83", "cxcr4a", "cd74b", "cd74a", "cd37", "zfp36l1a")  # lymphs
jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)

jens.choose <- jens.choose[which(!is.na(jens.choose))]
genesets[["lymph"]] <- jens.choose

# monocytes
jgenes.choose1 <- c("adam8a", "odc1", "lta4h", "thy1", "scpp8", "illr4", "timp2b", "mmp9", "mmp13a", "scinlb")  # monocytes
jgenes.choose2 <- c("cpa5", "lyz", "lect2l", "npsn", "sms", "abcb9", "ch25hl2", "papss2b", "hsd3b7", "cfd")  # neutros
jgenes.choose <- c(jgenes.choose1, jgenes.choose2)
jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)


jens.choose <- jens.choose[which(!is.na(jens.choose))]
genesets[["granu"]] <- jens.choose

# eryths
jgenes.choose <- c("rhag", "prdx2", "epor", "gata1a", "tspo", "slc4a1a", "sptb", "cahz", "hbba1", "hbba2", "alas2", "epb41b", "nt5c2l1")
jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)

# ba1 is hbb-like, maybe another name? chr3:55,098,010-55,099,718 use hbba1 and hbba2 instead
# jsubset(annot.out$out2.df, seqnames == "chr3" & start > 55000000 & end < 55500000)

jens.choose <- jens.choose[which(!is.na(jens.choose))]
genesets[["eryth"]] <- jens.choose

jnames <- names(genesets)
names(jnames) <- jnames

genesets.dat <- lapply(jnames, function(jname){
  jgenes <- genesets[[jname]]
  data.frame(geneset = jname, ens = genesets[[jname]], stringsAsFactors = FALSE)
}) %>%
  bind_rows()


# Add annotations ---------------------------------------------------------

jmerged.annot <- left_join(jmerged, genes.annot) %>%
  left_join(., genesets.dat)

jmerged.annot$geneset[is.na(jmerged.annot$geneset)] <- "aOther"

jmerged.annot <- jmerged.annot %>%
  arrange(geneset)


ggplot(jmerged.annot, aes(x = H3K4me3, y = H3K27me3, color = geneset)) + 
  geom_point_rast(size = jsize) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster) + 
  scale_color_manual(values = cbPalette) + 
  geom_density_2d(color = "black", alpha = 0.5)

ggplot(jmerged.annot, aes(x = H3K4me1, y = H3K27me3, color = geneset)) + 
  geom_point_rast(size = jsize) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster) + 
  scale_color_manual(values = cbPalette) + 
  geom_density_2d(color = "black", alpha = 0.5)

ggplot(jmerged.annot, aes(x = H3K4me3, y = H3K4me1, color = geneset)) + 
  geom_point_rast(size = jsize) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster) + 
  scale_color_manual(values = cbPalette) + 
  geom_density_2d(color = "black", alpha = 0.5)


dev.off()