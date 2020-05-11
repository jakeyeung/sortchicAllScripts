t# Jake Yeung
# Date of Creation: 2020-05-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/8-finalize_downsampled_TSS_counts_integrated_all_3_marks_readsDownsampled.DE_genes_like_BM.R
# 


rm(list=ls())

jstart <- Sys.time()

set.seed(0)

library(DropletUtils)

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

jwinsize <- "20000"


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
jvarcutoffs <- c(0.75, 2, 1, 0.5)
names(jvarcutoffs) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

jchromos <- paste("chr", seq(25), sep = "")

jsize <- 4
cbPalette <- c("grey85", "#32CD32", "#56B4E9", "#CC79A7", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

downsample.reads <- TRUE
downsample.cells <- FALSE

if (downsample.cells){
  min.cells.keep <- 100
}
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/pseudobulk_analysis.integrated_analysis.TSS.readsDS.likeBM.OtherWinsizes"
dir.create(outdir)

outprefix <- paste0("integrated_analysis.", Sys.Date(), ".UseTSSfromH3K4me3.DSreads_", downsample.reads, ".DScells_", downsample.cells, ".likeBM.winsize_", jwinsize)
outname <- paste0(outprefix, ".pdf")
outname.rdata <- paste0(outprefix, ".RData")

outpdf <- file.path(outdir, outname)
outrdata <- file.path(outdir, outname.rdata)

pdf(outpdf, useDingbats = FALSE)


# Read TSS Signal to figure out which transcript to keep  -----------------

indir.tss <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.TSS.CodingOnly.imputevarfilt.lessstringent.mapq_40.winsize_", jwinsize)
# indir.tss <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.TSS.CodingOnly.winsize_", jwinsize)
assertthat::assert_that(dir.exists(indir.tss))

tss.out <- lapply(jmarks, function(jmark){
  print(jmark)
  # inf.tss <- file.path(indir.tss, paste0("PZ-ChIC-ZF_", jmark, "_2020-04-07.countTable.TSS.csv"))
  inf.tss <- file.path(indir.tss, paste0(jmark, ".imputevarfilt.lessstringent.mapq_40.remerged.countTable.TSS.csv"))
  mat.tss <- ReadMatTSSFormat(inf.tss)
  # mat.tss <- CollapseRowsByGene(mat.tss, as.long=FALSE, track.kept.gene = TRUE)  # do this for just one mark
  return(list(mat.tss = mat.tss, tss.exprs = rowSums(mat.tss)))
})

tss.exprs.lst.unfilt <- lapply(tss.out, function(x) x$tss.exprs)
tss.mats.singlecell.unfilt <- lapply(tss.out, function(x) x$mat.tss)



# exprs.vec <- tss.exprs.lst$H3K4me1

lapply(jmarks, function(jmark){
  plot(density(tss.exprs.lst.unfilt[[jmark]]), main = jmark)
})


# go with the K4me3 definition...

ref.mark <- "H3K4me3"
jthres <- 50  # maybe not exactly at hump? what about tissuespecific stuff? rare celltypes? complicated from the bulk 
plot(density(tss.exprs.lst.unfilt[[ref.mark]]))
abline(v = jthres)

tss.mat.ref <- CollapseRowsByGene(count.mat = tss.mats.singlecell.unfilt[[ref.mark]], as.long = FALSE, track.kept.gene = TRUE)

tss.keep <- rownames(tss.mat.ref)

# # pick TSS's that are active, the rest are thrown out 
# tss.keep.i <- which(tss.exprs.lst[[ref.mark]] >= jthres)
# tss.keep <- tss.exprs.lst[[ref.mark]][tss.keep.i]

tss.exprs.lst <- lapply(tss.exprs.lst.unfilt, function(exprs.vec){
  jkeep <- names(exprs.vec) %in% tss.keep
  return(exprs.vec[jkeep])
})

print("Dimensions of TSS raw keeping all TSS")
lapply(tss.mats.singlecell.unfilt, length)
tss.mats.singlecell <- lapply(tss.mats.singlecell.unfilt, function(tss.mat){
  jkeep <- rownames(tss.mat) %in% tss.keep
  return(tss.mat[jkeep, ])
})

print("Dimensions of TSS after keeping one TSS for each gene, defined by highest expression in H3K4me3")
lapply(tss.mats.singlecell, length)

# Get common rows ---------------------------------------------------------

lapply(tss.exprs.lst.unfilt, length)

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
if (downsample.cells){
  cnames.keep.lst.lst <- lapply(cnames.keep.lst.lst, function(cnames.by.mark){
    lapply(cnames.by.mark, function(cnames.by.cluster){
      return(sample(cnames.by.cluster, size = min.cells.keep, replace = FALSE))
    })
  })
}

lapply(cnames.keep.lst.lst, function(cnames.by.mark) lapply(cnames.by.mark, function(cnames.by.cluster) length(cnames.by.cluster)))

tss.mats.pbulk <- lapply(jmarks, function(jmark){
  pbulk <- SumAcrossClusters(tss.mats.singlecell.filt[[jmark]], cnames.keep.lst.lst[[jmark]])
  return(do.call(cbind, pbulk))
})

for (jmark in jmarks){
  boxplot(tss.mats.pbulk[[jmark]], main = jmark)
  # boxplot(log10(tss.mats.pbulk + 1), main = jmark)
}

# Further downsample? -----------------------------------------------------


if (downsample.reads){
  # down sample

  # do further down sampling?
  ncuts.by.pbulk <- lapply(tss.mats.pbulk, function(jmat) colSums(jmat))
  
  # min across all? yes
  min.reads <- min(unlist(ncuts.by.pbulk))


  tss.mats.pbulk.ds <- lapply(tss.mats.pbulk, function(jmat){
    jprop <- min.reads / colSums(jmat)
    jmat.ds <- DropletUtils::downsampleMatrix(jmat, jprop, bycol = TRUE)
    return(jmat.ds)
  })
  (ncuts.by.pbulk.ds <- lapply(tss.mats.pbulk.ds, function(jmat) colSums(jmat)))
  tss.mats.pbulk <- tss.mats.pbulk.ds
  
  for (jmark in jmarks){
    boxplot(tss.mats.pbulk[[jmark]], main = paste("Reads downsampled:", jmark))
    # boxplot(log10(tss.mats.pbulk + 1), main = paste("Reads downsampled:", jmark))
  }
}



# Make long  --------------------------------------------------------------

jname.hash <- hash::hash(c("eryth1", "HSC1"), c("eryth", "HSC"))
jnames.keep <- c("eryth", "HSC", "lymph", "monocyte")

jlong.lst <- lapply(jmarks, function(jmark){
  jmat.tmp <- tss.mats.pbulk[[jmark]]
  dat.tss <- data.frame(bin = rownames(jmat.tmp), jmat.tmp, stringsAsFactors = FALSE) %>%
    reshape2::melt(., id.var = "bin", variable.name = "cluster", value.name = "counts") %>%
    rowwise() %>%
    mutate(cluster = as.character(cluster), 
           cluster = ifelse(is.null(jname.hash[[cluster]]), cluster, jname.hash[[cluster]])) %>%
    ungroup() %>%
    filter(cluster %in% jnames.keep)
  return(dat.tss)
})


# Plot K4me3 vs K27me3  ---------------------------------------------------



jmerged <- lapply(jmarks, function(jmark){
  print(jmark)
  jtmp <- jlong.lst[[jmark]]
    # rowwise() %>%
    # mutate(cluster = as.character(cluster), 
    #        cluster = ifelse(is.null(jname.hash[[cluster]]), cluster, jname.hash[[cluster]])) %>%
    # ungroup() %>%
    # filter(cluster %in% jnames.keep)
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

inf.zf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/pseudobulk_analysis.get_DE_genes_from_pbulk_scrnaseq/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.exprsmax_5.logfcmin_2.logfcmax_1.RData"
load(inf.zf.de, v=T)

genesets <- de.ens.zf.stringent


# inf.kob.erythmyeloid <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Kobayashi_2019_zf_genelists/erythmyeloid_genes_list.txt"
# kobayashi.genes.erythmyeloid <- fread(inf.kob.erythmyeloid, header = FALSE)$V1
# jens.erythmyeloid <- kobayashi.genes.erythmyeloid
# 
# inf.kob <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Kobayashi_2019_zf_genelists/hspc_genes_list.txt"
# kobayashi.genes <- fread(inf.kob, header = FALSE)$V1
# jens.hspc <- kobayashi.genes
# # HSPCs
# jgenes.choose <- c("meis1b", "pmp22b", "ahnak", "krt8", "anxa5b", "mrc1b", "pa2g4a", "npm1a", "ahcy", "adh5", "fabp3", "myb")
# # add kobayashi?
# jgenes.choose <- c(jgenes.choose, kobayashi.genes)
# jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)
# 
# jens.choose <- jens.choose[which(!is.na(jens.choose))]
# 
# genesets[["HSC"]] <- jens.choose
# 
# # lymphs
# jgenes.choose <- c("pax5", "cd79a", "bhlhe40", "cd83", "cxcr4a", "cd74b", "cd74a", "cd37", "zfp36l1a")  # lymphs
# jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)
# 
# jens.choose <- jens.choose[which(!is.na(jens.choose))]
# genesets[["lymph"]] <- jens.choose
# 
# # monocytes
# jgenes.choose1 <- c("adam8a", "odc1", "lta4h", "thy1", "scpp8", "illr4", "timp2b", "mmp9", "mmp13a", "scinlb")  # monocytes
# jgenes.choose2 <- c("cpa5", "lyz", "lect2l", "npsn", "sms", "abcb9", "ch25hl2", "papss2b", "hsd3b7", "cfd")  # neutros
# jgenes.choose <- c(jgenes.choose1, jgenes.choose2)
# jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)
# 
# 
# jens.choose <- jens.choose[which(!is.na(jens.choose))]
# genesets[["granu"]] <- jens.choose
# 
# # eryths
# jgenes.choose <- c("rhag", "prdx2", "epor", "gata1a", "tspo", "slc4a1a", "sptb", "cahz", "hbba1", "hbba2", "alas2", "epb41b", "nt5c2l1")
# jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)
# 
# # ba1 is hbb-like, maybe another name? chr3:55,098,010-55,099,718 use hbba1 and hbba2 instead
# # jsubset(annot.out$out2.df, seqnames == "chr3" & start > 55000000 & end < 55500000)
# 
# jens.choose <- jens.choose[which(!is.na(jens.choose))]
# genesets[["eryth"]] <- jens.choose

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

jmerged.annot$geneset[is.na(jmerged.annot$geneset)] <- "zOther"

jmerged.annot <- jmerged.annot %>%
  arrange(desc(geneset))


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

# plot marginals
m.dens <- lapply(jmarks, function(jmark){
  m <- ggplot(jmerged.annot, aes_string(x = jmark, fill = "cluster")) + 
    geom_density(alpha = 0.3) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~cluster) + 
    scale_color_manual(values = cbPalette)
  return(m)
})
print(m.dens)

m.hist <- lapply(jmarks, function(jmark){
  m <- ggplot(jmerged.annot, aes_string(x = jmark, fill = "cluster")) + 
    geom_histogram(alpha = 0.3) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~cluster) + 
    scale_color_manual(values = cbPalette)
  return(m)
})
print(m.hist)

m.dens.log <- lapply(jmarks, function(jmark){
  m <- ggplot(jmerged.annot, aes_string(x = jmark, fill = "cluster")) + 
    geom_histogram(alpha = 0.3) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~cluster) + 
    scale_color_manual(values = cbPalette) + 
    scale_x_log10() 
  return(m)
})
print(m.dens.log)


m.dens.sqrt <- lapply(jmarks, function(jmark){
  m <- ggplot(jmerged, aes_string(x = paste0("sqrt(", jmark, ")"), fill = "cluster")) + 
    geom_histogram() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~cluster, ncol = 1)
  return(m)
})
print(m.dens.sqrt)


# Create differential analysis --------------------------------------------

jlong.diff <- lapply(jmarks, function(jmark){
  jlong <- jlong.lst[[jmark]]
  jtmp <- jlong %>%
    group_by(bin) %>%
    mutate(log2p1counts = log2(counts + 1),
           log2fc = log2p1counts - mean(log2p1counts),
           zscore = log2fc / sd(log2p1counts))
  jtmp$gene <- sapply(jtmp$bin, function(b) strsplit(b, ";")[[1]][[2]])
  jtmp$ens <- sapply(jtmp$gene, function(g) AssignHash(g, g2e, null.fill = g))
  jtmp$mark <- jmark
  return(jtmp)
}) %>%
  bind_rows()

jlong.diff.genesets <- left_join(jlong.diff, genesets.dat)
jlong.diff.genesets$geneset[is.na(jlong.diff.genesets$geneset)] <- "zOther"

# jlong.diff.genesets$cluster <- jlong.diff.genesets$cluster.new


ggplot(jlong.diff.genesets, aes(x = cluster, y = zscore, fill = mark)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset) + 
  geom_hline(yintercept = 0, linetype = "dotted")

ggplot(jlong.diff.genesets, aes(x = cluster, y = log2p1counts, fill = mark)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset) 

ggplot(jlong.diff.genesets, aes(x = cluster, y = sqrt(counts), fill = mark)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset) 

jmerged.log2fc <- reshape2::dcast(jlong.diff, formula = "bin + cluster ~ mark", value.var = "log2fc")
jmerged.zscore <- reshape2::dcast(jlong.diff, formula = "bin + cluster ~ mark", value.var = "zscore")

ggplot(jmerged, aes(x = H3K4me3, y = H3K27me3, color = cluster)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~cluster) 

ggplot(jmerged, aes(x = H3K4me1, y = H3K27me3, color = cluster)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~cluster) 

ggplot(jmerged, aes(x = H3K4me3, y = H3K4me1, color = cluster)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~cluster) 

ggplot(jmerged.log2fc, aes(x = H3K4me3, y = H3K27me3, color = cluster)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~cluster) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

ggplot(jmerged.log2fc, aes(x = H3K4me1, y = H3K27me3, color = cluster)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~cluster) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

ggplot(jmerged.log2fc, aes(x = H3K4me3, y = H3K4me1, color = cluster)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~cluster) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

ggplot(jmerged.zscore, aes(x = H3K4me3, y = H3K27me3, color = cluster)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~cluster) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

ggplot(jmerged.zscore, aes(x = H3K4me1, y = H3K27me3, color = cluster)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~cluster) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

ggplot(jmerged.zscore, aes(x = H3K4me3, y = H3K4me1, color = cluster)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~cluster) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)




dev.off()

print(Sys.time() - jstart)

# save objects
save(jmerged.annot, genesets, jlong.lst, tss.mats.pbulk, cnames.keep.lst.lst, annot.glmpca.filt.lst, jlong.diff.genesets, file = outrdata)