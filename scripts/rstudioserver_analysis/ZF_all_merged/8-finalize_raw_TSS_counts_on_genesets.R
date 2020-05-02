# Jake Yeung
# Date of Creation: 2020-05-02
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/8-finalize_raw_TSS_counts_on_genesets.R
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


# Constants ---------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
jvarcutoffs <- c(0.75, 2, 1, 0.5)
names(jvarcutoffs) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

jchromos <- paste("chr", seq(25), sep = "")

# Read count matrix -------------------------------------------------------

# jbinsize <- 5000
# indir.cmat <- file.path(hubprefix, "jyeung/data/zebrafish_scchic/count_tables.winsize_5000")
indir.cmat <- file.path(hubprefix, "jyeung/data/zebrafish_scchic/count_tables.winsize_5000")
# indir.cmat <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.TSS.winsize_", jbinsize)
assertthat::assert_that(dir.exists(indir.cmat))


count.mats <- lapply(jmarks, function(jmark){
  inf.cmat <- file.path(indir.cmat, paste0("PZ-ChIC-ZF_", jmark, "_2020-04-07.countTable.csv"))
  # inf.cmat <- file.path(indir.cmat, paste0("PZ-ChIC-ZF_", jmark, "_2020-04-07.countTable.TSS.csv"))
  cmat <- ReadMatSlideWinFormat(inf.cmat)
  # keep bins starting with chromos
  return(cmat)
})

lapply(count.mats, dim)

bins.all <- lapply(count.mats, function(count.mat){
  bins <- rownames(count.mat)
}) %>%
  Reduce(f = intersect, .)
  

# keep bins in chromosomes
chromos.all <- sapply(bins.all, GetChromo)
chromos.keep <- chromos.all %in% jchromos
bins.keep <- bins.all[chromos.keep]

print(length(bins.all))
print(length(bins.keep))


# # Read TSS Signal to figure out which transcript to keep  -----------------
# 
# indir.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.TSS.CodingOnly.winsize_10000"
# assertthat::assert_that(dir.exists(indir.tss))
# 
# tss.exprs.lst <- lapply(jmarks, function(jmark){
#   print(jmark)
#   inf.tss <- file.path(indir.tss, paste0("PZ-ChIC-ZF_", jmark, "_2020-04-07.countTable.TSS.csv"))
#   mat.tss <- ReadMatTSSFormat(inf.tss)
#   mat.tss <- CollapseRowsByGene(mat.tss, as.long=FALSE, track.kept.gene = TRUE)
#   return(rowSums(mat.tss))
# })
# 
# # exprs.vec <- tss.exprs.lst$H3K4me1
# 
# lapply(jmarks, function(jmark){
#   plot(density(tss.exprs.lst[[jmark]]), main = jmark)
# })
# 
# # go with the K4me3 definition...
# 
# ref.mark <- "H3K4me3"
# jthres <- 50  # maybe not exactly at hump? what about tissuespecific stuff? rare celltypes? complicated from the bulk 
# plot(density(tss.exprs.lst[[ref.mark]]))
# abline(v = jthres)
# 
# # pick TSS's that are active, the rest are thrown out 
# tss.keep.i <- which(tss.exprs.lst[[ref.mark]] >= jthres)
# tss.keep <- tss.exprs.lst[[ref.mark]][tss.keep.i]


# Load pseudobullks  ------------------------------------------------------

ref.mark <- "H3K4me3"
jthres.log <- 0.7

pbulk.chic.mat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.pbulk.chic <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/rdata_robjs/WKM_scChICseq_pseudobulk_downsampled.ens.", jmark, ".RData")
  load(inf.pbulk.chic, v=T)
  pbulk.filt.ds <- CollapseRowsByGene(pbulk.filt.ds, as.long = FALSE, track.kept.gene = TRUE)
  return(pbulk.filt.ds)
})

jnames <- rownames(pbulk.chic.mat.lst[[ref.mark]])
exprs.lst <- lapply(as.list(as.data.frame(pbulk.chic.mat.lst[[ref.mark]])), function(x){
  names(x) <- jnames
  return(x)
})


lapply(exprs.lst, function(jexprs){
  plot(density(log10(jexprs + 1)))
  abline(v = jthres.log)
})

# kep TSS's that have greater than threshold
tss.keep.by.pbulk <- lapply(exprs.lst, function(x){
  xkeep <- x >= jthres.log
  tss.keep.tmp <- names(x)[xkeep]
  return(tss.keep.tmp)
})

tss.keep.final <- unlist(tss.keep.by.pbulk) %>%
  unique()

gene.txi.final <- sapply(tss.keep.final, function(x) paste(strsplit(x, ";")[[1]][2:3], collapse = ";"))


# Annotate bins  ----------------------------------------------------------

jwin <- 500
# take any mark
inf.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/zebrafish/gene_tss.CodingOnly.winsize_", jwin, ".species_drerio.bed")
assertthat::assert_that(file.exists(inf.annot))
jchromos <- paste("chr", seq(25), sep = "")

annot.out <- AnnotateCoordsFromList.GeneWise(bins.keep, inf.tss = inf.annot, txdb = TxDb.Drerio.UCSC.danRer11.refGene, annodb = "org.Dr.eg.db", chromos.keep = jchromos)

# take just the first TSS? 

annot.regions <- subset(annot.out$out2.df, select = c(dist.to.tss, region_coord, gene, tssname)) %>%
  filter(gene != "")

# add ensembl
g2e.dat <- data.frame(region_coord = as.character(annot.out$regions.annotated$region_coord), ens = as.character(annot.out$regions.annotated$ENSEMBL), stringsAsFactors = FALSE)
# remove duplicates
g2e.dat <- g2e.dat[!duplicated(g2e.dat), ]

annot.regions <- left_join(annot.regions, g2e.dat)

print(dim(annot.regions))
print(length(bins.keep))

print(paste("Fraction of bins with TSS annotations:", signif(nrow(annot.regions) / length(bins.keep), digits = 2)))



# 2 more rows now... because there is a gene matched to multiple ensembl IDs? 
  
# # add ensembl properly
# jgenes <- annot.regions$gene
# jens <- Gene2Ensembl.ZF(jgenes)


# # annot.regions.sub <- subset(annot.regions, distanceToTSS < 10000)

# why NAs?
genes.na <- unique(subset(annot.regions, is.na(ens))$gene)
print(length(genes.na))
ens.na <- Gene2Ensembl.ZF(genes.na)
print(length(ens.na))
jhash.na <- hash::hash(genes.na, ens.na)

annot.regions <- annot.regions %>%
  rowwise() %>%
  mutate(ens = ifelse(is.na(ens), AssignHash(gene, jhash.na), ens))


# genes.split <- sapply(annot.regions$tssname, function(x) strsplit(x, ";")[[1]][[2]], USE.NAMES = FALSE)
# jgenes.unsplit <- annot.regions$gene
# assertthat::assert_that(identical(jgenes.split, jgenes.unsplit))

annot.regions$gene.txi <- sapply(annot.regions$tssname, function(x) paste(strsplit(x, ";")[[1]][2:3], collapse = ";"))

# filter out regions that are "active" tss
annot.regions.filt <- subset(annot.regions, gene.txi %in% gene.txi.final)

# # keep one region per gene?
# annot.regions.filt2 <- annots.regions.filt %>%
#   group_by(gene) %>%
#   top_n(., 1, -abs(dist.to.tss))



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


# Annotate bins -----------------------------------------------------------



cells.keep.lst <- lapply(annot.glmpca.filt.lst, function(x){
  return(x$cell)
})

# filter for good cells and good bins
count.mats.filt <- lapply(jmarks, function(jmark){
  jmat <- count.mats[[jmark]]
  cells.keep <- cells.keep.lst[[jmark]]
  jmat.filt <- jmat[bins.keep, cells.keep]
  return(jmat.filt)
})

# assign bins as TSS or not 
bins.tss <- unique(annot.regions.filt$region_coord)  # a region can have multiple genes
bins.outside <- bins.keep[which(!bins.keep %in% bins.tss)]
assertthat::assert_that(length(bins.tss) + length(bins.outside) == length(bins.keep))

bins.tss.dat <- data.frame(region_coord = bins.tss, is.tss = TRUE, stringsAsFactors = FALSE)
bins.outside.dat <- data.frame(region_coord = bins.outside, is.tss = FALSE, stringsAsFactors = FALSE)

bins.dat <- rbind(bins.tss.dat, bins.outside.dat)


# Naively plot TSS signal over all bins  ----------------------------------


tss.mats <- lapply(jmarks, function(jmark){
  jmat <- count.mats.filt[[jmark]]
  jmat[bins.tss, ]
})

outside.mats <- lapply(jmarks, function(jmark){
  jmat <- count.mats.filt[[jmark]]
  jmat[bins.outside, ]
})

lapply(tss.mats, dim)
lapply(outside.mats, dim)


total.counts <- lapply(count.mats.filt, function(count.mat) data.frame(cell = colnames(count.mat), total.cuts = colSums(count.mat), stringsAsFactors = FALSE))
tss.counts <- lapply(tss.mats, function(tss.mat) data.frame(cell = colnames(tss.mat), tss.cuts = colSums(tss.mat), stringsAsFactors = FALSE))
outside.counts <- lapply(outside.mats, function(outside.mat) data.frame(cell = colnames(outside.mat), outside.cuts = colSums(outside.mat), stringsAsFactors = FALSE))

dat.counts <- lapply(jmarks, function(jmark){
  dat.count <- left_join(total.counts[[jmark]], tss.counts[[jmark]]) %>%
    left_join(., outside.counts[[jmark]]) %>%
    left_join(., annot.glmpca.filt.lst[[jmark]])
})


# check for one mark
lapply(jmarks, function(jmark){
  m1 <- ggplot(dat.counts[[jmark]], aes(x = umap1, y = umap2, color = tss.cuts / total.cuts)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_viridis_c(direction = 1) + ggtitle(jmark)
  m2 <- ggplot(dat.counts[[jmark]], aes(x = umap1, y = umap2, color = var.imputed)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_viridis_c(direction = -1) + ggtitle(jmark)
  multiplot(m1, m2, cols = 2)
})

lapply(jmarks, function(jmark){
  m1 <- ggplot(dat.counts[[jmark]], aes(x = umap1, y = umap2, color = tss.cuts / total.cuts / var.imputed)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_viridis_c(direction = 1) + ggtitle(jmark)
  m2 <- ggplot(dat.counts[[jmark]], aes(x = umap1, y = umap2, color = var.imputed)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_viridis_c(direction = -1) + ggtitle(jmark)
  multiplot(m1, m2, cols = 2)
})


lapply(jmarks, function(jmark){
  m1 <- ggplot(dat.counts[[jmark]], aes(x = umap1, y = umap2, color = log10(outside.cuts))) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_viridis_c(direction = -1) + ggtitle(jmark)
  print(m1)
})

# Check TSS signal for outliers -------------------------------------------


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
  pbulk <- SumAcrossClusters(tss.mats[[jmark]], cnames.keep.lst.lst[[jmark]])
  return(do.call(cbind, pbulk))
})

outside.mats.pbulk <- lapply(jmarks, function(jmark){
  pbulk <- SumAcrossClusters(outside.mats[[jmark]], cnames.keep.lst.lst[[jmark]])
  return(do.call(cbind, pbulk))
})

total.mats.pbulk <- lapply(jmarks, function(jmark){
  pbulk <- SumAcrossClusters(count.mats.filt[[jmark]], cnames.keep.lst.lst[[jmark]])
  return(do.call(cbind, pbulk))
})

jlong.lst <- lapply(jmarks, function (jmark){
  jmat.lst <- list(tss.cuts = tss.mats.pbulk, outside.cuts = outside.mats.pbulk, total.cuts = total.mats.pbulk)
  jmat.names <- names(jmat.lst)
  names(jmat.names) <- jmat.names
  dat.lst <- lapply(jmat.names, function(jmat.name){
    jmat.tmp <- jmat.lst[[jmat.name]][[jmark]]
    dat.tss <- data.frame(bin = rownames(jmat.tmp), jmat.tmp, stringsAsFactors = FALSE) %>%
      reshape2::melt(., id.var = "bin", variable.name = "cluster", value.name = "counts")
  })
})

jmark.test <- "H3K27me3"
jfeature <- "tss.cuts"
ggplot(jlong.lst[[jmark.test]][[jfeature]], aes(x = log10(counts + 1), fill = cluster)) + geom_density() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster, ncol = 1) + ggtitle(jmark.test, jfeature)


# normalize by number of cells?
# ncells.by.pbulk.lst <- lapply(jmarks, function(jmark){
#   annot.sum <- annot.glmpca.filt.lst[[jmark]] %>%
#     group_by(cluster) %>%
#     summarise(ncells = length(cell))
# })
# counts.by.pbulk.lst <- lapply(jmarks, function(jmark){
#   dat.out <- data.frame(cluster = colnames(total.mats.pbulk[[jmark]]), N = colSums(total.mats.pbulk[[jmark]]), stringsAsFactors = FALSE)
# })
# jlong.lst.annot <- lapply(jmarks, function(jmark){
#   lapply(jlong.lst[[jmark]], function(jtmp){
#     left_join(jtmp, ncells.by.pbulk.lst[[jmark]]) %>%
#       left_join(., counts.by.pbulk.lst[[jmark]])
#   })
# })

jmark.test <- "H3K4me1"
jfeature <- "tss.cuts"

ggplot(jlong.lst[[jmark.test]][[jfeature]], aes(x = counts + 1, fill = cluster)) + geom_density() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster, ncol = 1) + ggtitle(jmark.test, jfeature) + scale_x_log10()

# barplot?
jlong.lst[[jmark.test]][[jfeature]] %>%
  group_by(cluster) %>%
  summarise(counts = sum(counts)) %>% 
  ggplot(., aes(x = cluster, y = counts)) + geom_col()

# ncells greater than threshold?
jlong.lst[[jmark.test]][[jfeature]] %>%
  group_by(cluster) %>%
  summarise(ngenes.great.than.thres = length(which(counts > 3)))


# ggplot(jlong.lst.annot[[jmark.test]][[jfeature]], aes(x = counts/N * 10^6 + 1, fill = cluster)) + geom_density() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~cluster, ncol = 1) + ggtitle(jmark.test, jfeature) + scale_x_log10()


# There are some regions with TONS of cuts in eryth?  ---------------------

# check single cells

# count mat by pseudobulks??
jmark.test <- "H3K27me3"
jmat.by.pbulk <- lapply(cnames.keep.lst.lst[[jmark.test]], function(cells.keep){
  jmat.tmp <- count.mats.filt[[jmark.test]][, cells.keep]
})
lapply(jmat.by.pbulk, sum)

lapply(jmat.by.pbulk, dim)

# lapply(jmat.by.pbulk, function(x) table(as.matrix(x)))

# show the UMAP of the 100 cells?


# check for one mark (100 cells each, but we still see this effect)
lapply(jmarks, function(jmark){
  m1 <- ggplot(dat.counts[[jmark]] %>% filter(cell %in% unlist(cnames.keep.lst.lst[[jmark]])), aes(x = umap1, y = umap2, color = tss.cuts / total.cuts)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_viridis_c(direction = 1) + ggtitle(jmark)
  print(m1)
})


# Plot K4me3 vs K27me3  ---------------------------------------------------

jfeature <- "tss.cuts"

jname.hash <- hash::hash(c("eryth1", "HSC1"), c("eryth", "HSC"))

jnames.keep <- c("eryth", "HSC", "lymph", "monocyte")

jmerged <- lapply(jmarks, function(jmark){
  print(jmark)
  jtmp <- jlong.lst[[jmark]]$tss.cuts %>%
    rowwise() %>%
    mutate(cluster = as.character(cluster), 
           cluster = ifelse(is.null(jname.hash[[cluster]]), cluster, jname.hash[[cluster]])) %>%
    ungroup() %>%
    filter(cluster %in% jnames.keep)
  colnames(jtmp) <- c("bin", "cluster", jmark)
  return(jtmp)
}) %>%
  purrr::reduce(., left_join)

ggplot(jmerged, aes(x = H3K4me3, y = H3K27me3, color = cluster)) + geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster)

ggplot(jmerged, aes(x = H3K4me1, y = H3K27me3, color = cluster)) + geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster)

ggplot(jmerged, aes(x = H3K4me3, fill = cluster)) + geom_density(alpha = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster, ncol = 1)

ggplot(jmerged, aes(x = H3K27me3, fill = cluster)) + geom_density(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster, ncol = 1)
