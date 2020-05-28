# Jake Yeung
# Date of Creation: 2020-05-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/redo_multiomics_integration_TSS_pseudobulk_stringent_DE.R
# 

# After analyzing ZF, redo the multiomics analysis  
# Things to consider: downsampling (100 cells), assigning bin to gene in a gene-wise manner, duplicating DE genes in different gene sets

rm(list=ls())

jstart <- Sys.time()

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(ggrastr)

  library(DropletUtils)


make.plots <- TRUE


# Load DE genes -----------------------------------------------------------

# load this first because it loads a lot of objects, might disuprt things

inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData"
load(inf.de, v=T)

# Load TSS  ---------------------------------------------------------------

fewer.k27me3 <- TRUE
downsample.reads <- TRUE
downsample.cells <- FALSE

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.stringentDE"
dir.create(outdir, showWarnings = TRUE)

outprefix <- file.path(outdir, paste0("integrated_analysis_3_marks.stringentDE.WithHighLowExprs.pseudobulk.fewerk27me3_", fewer.k27me3, ".DSreads_", downsample.reads, "DScells.", downsample.cells, ".", Sys.Date()))

outpdf <- paste0(outprefix, ".pdf")
outrdata <- paste0(outprefix, ".RData")

set.seed(0)
cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks


# we did 10kb for zf
jwinsize <- 10000L
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS"

tss.mats <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0(jmark, ".countTableTSS.mapq_40.TSS_", jwinsize, ".blfiltered.csv")
  inf <- file.path(indir, fname)
  mat <- ReadMatTSSFormat(inf, as.sparse = TRUE, add.coord = FALSE, sort.rnames = TRUE)
  return(mat)
})

ref.mark <- "H3K4me3"
tss.mat.ref.filt <- CollapseRowsByGene(tss.mats[[ref.mark]], as.long = FALSE, track.kept.gene = TRUE)

tss.keep.fromref <- rownames(tss.mat.ref.filt)

# filter TSS for H3K4me3 as reference, then all the other marks will take the same TSS (TODO maybe do this for ZF also)
tss.mats.filt.fromref <- lapply(tss.mats, function(jmat){
  jkeep <- rownames(jmat) %in% tss.keep.fromref
  return(jmat[jkeep, ])
})


# Keep common rows --------------------------------------------------------


lapply(tss.mats.filt.fromref, length)

tss.all <- lapply(tss.mats.filt.fromref, function(jmat){
  rownames(jmat)
}) %>%
  unlist() %>%
  unique()

tss.common <- lapply(tss.mats.filt.fromref, function(jmat){
  rownames(jmat) 
}) %>%
  Reduce(f = intersect, .)

# filter TSS for H3K4me3 as reference, then all the other marks will take the same TSS (TODO maybe do this for ZF also)
tss.mats.filt <- lapply(tss.mats.filt.fromref, function(jmat){
  jkeep <- rownames(jmat) %in% tss.keep.fromref
  return(jmat[jkeep, ])
})


# get ensembl names ? 
genes.common <- sapply(tss.common, function(x) strsplit(x, ";")[[1]][[2]])
ens.common <- Gene2Ensembl.ZF(genes.common, return.original = TRUE, species = "mmusculus")
g2e <- hash::hash(genes.common, ens.common)

# create tss, genes, ens dat
genes.annot <- data.frame(bin = tss.common, gene = genes.common, ens = ens.common, stringsAsFactors = FALSE)


# Load cell cluster annots ------------------------------------------------

indir.annots <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping")
dat.annots.all <- lapply(jmarks, function(jmark){
  inf.annots <- file.path(indir.annots, paste0("GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData"))
  load(inf.annots, v=T)
  return(dat.umap.glm.fillNAs)
})


# Show UMAPs --------------------------------------------------------------

if (make.plots){
  pdf(outpdf, useDingbats = FALSE)
}

m.umaps <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.annots.all[[jmark]], aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_manual(values = cbPalette)
  return(m)
})

JFuncs::multiplot(m.umaps[[1]], m.umaps[[2]], m.umaps[[3]], cols = 3)

# check all possible clusters
lapply(dat.annots.all, function(x) sort(unique(x$cluster)))


# let's keep some common celltypes 


# Bcells, Eryth, HSCs, Neutrophils (aligns with zebrafish)
clstrs.all <- lapply(dat.annots.all, function(x) sort(unique(x$cluster))) %>%
  unlist()
print(length(clstrs.all))

# clstrs.remove <- c("Eryth-Gfi1-_topic17", "HSCs-Lrp5_topic14", "HSCs-Msi2_topic20", "HSCs-Ephb2_topic5")
clstrs.remove <- c("HSCs-Lrp5_topic14", "HSCs-Msi2_topic20", "HSCs-Ephb2_topic5")  # K27me3 has a questionable eryth cluster, Gfi1, but let's keep it otherwise only 152 eryth

clstrs.all <- clstrs.all[!clstrs.all %in% clstrs.remove]
print(length(clstrs.all))

bcells.names <- grep("^Bcell", clstrs.all, ignore.case = TRUE, value = TRUE)
neutro.names <- grep("^Neut", clstrs.all, ignore.case = TRUE, value = TRUE)
eryth.names <- grep("^Eryth", clstrs.all, ignore.case = TRUE, value = TRUE)
hsc.names <- grep("^HSC", clstrs.all, ignore.case = TRUE, value = TRUE)

clstr.names <- c("Bcells", "Granulocytes", "Erythroblasts", "HSPCs")
names(clstr.names) <- clstr.names
clstrs.filt.lst <- list(bcells.names, neutro.names, eryth.names, hsc.names)
names(clstrs.filt.lst) <- clstr.names

# replot maps 
m.umaps.filt <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.annots.all[[jmark]] %>% filter(cluster %in% unlist(clstrs.filt.lst)), aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_manual(values = cbPalette)
  return(m)
})
JFuncs::multiplot(m.umaps.filt[[1]], m.umaps.filt[[2]], m.umaps.filt[[3]], cols = 3)

clstrs.filt.dat <- lapply(clstr.names, function(cname){
  jdat <- data.frame(cluster.new = cname, cluster = clstrs.filt.lst[[cname]], stringsAsFactors = FALSE)
  return(jdat)
}) %>%
  bind_rows()

# Reannotate so that filtered cells have same cluster names  --------------

dat.annots.filt <- lapply(jmarks, function(jmark){
  jdat.new <- left_join(dat.annots.all[[jmark]], clstrs.filt.dat) %>%
    filter(!is.na(cluster.new))
})

lapply(dat.annots.all, nrow)
lapply(dat.annots.filt, nrow)

# check number of cells kept
m.umaps.beforeafter <- lapply(jmarks, function(jmark){
  n.before <- paste("N:", nrow(dat.annots.all[[jmark]]))
  n.after <- paste("N:", nrow(dat.annots.filt[[jmark]]))
  multiplot(m.umaps[[jmark]] + ggtitle(n.before), m.umaps.filt[[jmark]] + ggtitle(n.after), cols = 2)
})


# go to smallest number of cells (across all marks and clusters) ---------

annots.split <- lapply(jmarks, function(jmark){
  cells.vec <- lapply(split(dat.annots.filt[[jmark]], f = dat.annots.filt[[jmark]]$cluster.new), function(jout2){
    print(paste(jmark, unique(jout2$cluster.new)))
    print(nrow(jout2))
    return(jout2$cell)
  })
})

ncells.by.clstr <- lapply(jmarks, function(jmark){
  lapply(annots.split[[jmark]], function(cells.vec){
    return(length(cells.vec))
  })
})

print(ncells.by.clstr)
hist(unlist(unlist(ncells.by.clstr)))
(min.ncells <- min(unlist(unlist(ncells.by.clstr))))  # Bcells in K27me3 smallest 


# Downsmaple by lowest cells --------------------------------------------------------------

if (downsample.cells){
  cnames.keep.lst.lst <- lapply(annots.split, function(cnames.by.mark){
    lapply(cnames.by.mark, function(cnames.by.cluster){
      return(sample(cnames.by.cluster, size = min.ncells, replace = FALSE))
    })
  })
} else {
  cnames.keep.lst.lst <- lapply(annots.split, function(cnames.by.mark){
    lapply(cnames.by.mark, function(cnames.by.cluster){
      return(cnames.by.cluster)  # do nothing, usually down sample reads will take care of this later 
    })
  })
}

ncells.by.clstr.after <- lapply(jmarks, function(jmark){
  lapply(cnames.keep.lst.lst[[jmark]], function(cells.vec){
    return(length(cells.vec))
  })
})

print(ncells.by.clstr)
print(ncells.by.clstr.after)


# Show cell fractions?  ---------------------------------------------------

cells.downsampled.keep <-unlist(unlist(cnames.keep.lst.lst)) 
  
dat.annots.filt.downsamp.long <- lapply(jmarks, function(jmark){
  jdat <- dat.annots.filt[[jmark]]
  jsub <- subset(jdat, cell %in% cells.downsampled.keep)
  jsub$mark <- jmark
  return(jsub)
}) %>%
  bind_rows()

ggplot(dat.annots.filt.downsamp.long, aes(x = cluster.new, group = cond, fill = cond)) + geom_bar() + 
  facet_wrap(~mark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


# Make pseudobulks  -------------------------------------------------------

tss.mats.pbulk <- lapply(jmarks, function(jmark){
  pbulk <- SumAcrossClusters(tss.mats.filt[[jmark]], cnames.keep.lst.lst[[jmark]])
  return(do.call(cbind, pbulk))
})

if (downsample.reads){
  
  
  # do further down sampling?
  ncuts.by.pbulk <- lapply(tss.mats.pbulk, function(jmat) colSums(jmat))
  
  # min across all? yes
  min.reads <- min(unlist(ncuts.by.pbulk))
  
  # down sample
  tss.mats.pbulk.ds <- lapply(tss.mats.pbulk, function(jmat){
    jprop <- min.reads / colSums(jmat)
    jmat.ds <- DropletUtils::downsampleMatrix(jmat, jprop, bycol = TRUE)
    return(jmat.ds)
  })
  
  (ncuts.by.pbulk.ds <- lapply(tss.mats.pbulk.ds, function(jmat) colSums(jmat)))
  
} else {
  tss.mats.pbulk.ds <- tss.mats.pbulk
}


for (jmark in jmarks){
  boxplot(tss.mats.pbulk[[jmark]], main = paste("before reads downsampling", jmark))
  boxplot(tss.mats.pbulk.ds[[jmark]], main = paste("after reads downsampling", jmark))
}


# Make long  --------------------------------------------------------------

jlong.lst <- lapply(jmarks, function(jmark){
  # jmat.tmp <- tss.mats.pbulk[[jmark]]
  jmat.tmp <- tss.mats.pbulk.ds[[jmark]]  # try downsample
  dat.tss <- data.frame(bin = rownames(jmat.tmp), jmat.tmp, stringsAsFactors = FALSE) %>%
    reshape2::melt(., id.var = "bin", variable.name = "cluster", value.name = "counts")
  dat.tss$cluster <- as.character(dat.tss$cluster)
  return(dat.tss)
})

# check very large counts
subset(jlong.lst$H3K27me3, counts > 150)

jmerged <- lapply(jmarks, function(jmark){
  print(jmark)
  jtmp <- jlong.lst[[jmark]]
  colnames(jtmp) <- c("bin", "cluster", jmark)
  return(jtmp)
}) %>%
  purrr::reduce(., left_join)

print(head(jmerged))

jsize <- 4


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

ggplot(jmerged, aes(x = sqrt(H3K4me3), y = sqrt(H3K27me3), color = cluster)) + 
  geom_point_rast(alpha = 0.25, size = jsize) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster) + 
  geom_density_2d(color = "black")

ggplot(jmerged, aes(x = sqrt(H3K4me1), y = sqrt(H3K27me3), color = cluster)) + 
  geom_point_rast(alpha = 0.25, size = jsize) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster) + 
  geom_density_2d(color = "black")



ggplot(jmerged, aes(x = H3K4me3, y = H3K4me1, color = cluster)) + 
  geom_point_rast(alpha = 0.25, size = jsize) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster) + 
  geom_density_2d(color = "black")

ggplot(jmerged, aes(x = log10(H3K27me3 + 1), fill = cluster)) + 
  # geom_density(alpha = 0.25) + 
  geom_histogram(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster, ncol = 1)

# plot marginals
m.dens <- lapply(jmarks, function(jmark){
  m <- ggplot(jmerged, aes_string(x = jmark, fill = "cluster")) + 
    geom_density() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~cluster, ncol = 1)
  return(m)
})
print(m.dens)

# m.dens.log <- lapply(jmarks, function(jmark){
#   m <- ggplot(jmerged, aes_string(x = jmark, fill = "cluster")) + 
#     geom_histogram() + 
#     scale_x_log10() + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     facet_wrap(~cluster, ncol = 1)
#   return(m)
# })
# print(m.dens.log)

m.dens.sqrt <- lapply(jmarks, function(jmark){
  m <- ggplot(jmerged, aes_string(x = paste0("sqrt(", jmark, ")"), fill = "cluster")) + 
    geom_histogram() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~cluster, ncol = 1)
  return(m)
})
print(m.dens.sqrt)


# Get gene sets -----------------------------------------------------------



de.ens.lst <- de.ens.sorted.stringent


# Add background genes ----------------------------------------------------

# add background genes 
de.ens.lst$zOther <- sample(dat.de.giladi.seurat$ens, size = 500)
jclsts.after <- names(de.ens.lst)
names(jclsts.after) <- jclsts.after



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

jlong.diff.genesets <- lapply(jclsts.after, function(jclst){
  jlong.sub <- subset(jlong.diff, ens %in% de.ens.lst[[jclst]]) %>%
        ungroup() %>%
        mutate(geneset = jclst)
  return(jlong.sub)
}) %>%
  bind_rows()


# # create different gene sets for boxplot, loop mark and cluster
# jlong.diff.genesets <- lapply(jmarks, function(jmark){
#   jlong.genesets <- lapply(jclsts.after, function(jclst){
#     jlong.sub <- subset(jlong.diff.lst[[jmark]], ens %in% de.ens.lst[[jclst]]) %>%
#       ungroup() %>%
#       mutate(geneset = jclst)
#     jlong.sub$mark <- jmark
#     return(jlong.sub)
#   }) %>%
#     bind_rows()
#   return(jlong.genesets)
# }) %>%
#   bind_rows()

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

ggplot(jmerged.log2fc, aes(x = H3K4me3, y = H3K27me3, color = cluster)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~cluster) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

ggplot(jmerged.zscore, aes(x = H3K4me3, y = H3K27me3, color = cluster)) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~cluster) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

# 
# # plot 2D of zscores and log2fc
# jmerged.counts <- subset(jmerged, select = c(-bin, -cluster))
# rownames(jmerged.counts) <- paste(jmerged$bin, jmerged$cluster, sep = "_")
# 
# jmerged.log2p1counts <- log2(jmerged.counts + 1)
# jmerged.log2fc <- sweep(jmerged.log2p1counts, MARGIN = 1, STATS = rowMeans(jmerged.log2p1counts), FUN = "-")
# jmerged.zscore <- sweep(jmerged.log2fc, MARGIN = 1, STATS = apply(jmerged.log2fc, 1, sd), FUN = "/")
# 
# # # check
# # jmerged.log2fc.check <- t(scale(t(jmerged.log2p1counts), center = TRUE, scale = FALSE))
# # jmerged.zscore.check <- t(scale(t(jmerged.log2p1counts), center = TRUE, scale = TRUE))
# 
# # add back rownames
# jmerged.log2fc <- data.frame(bin = jmerged$bin, cluster = jmerged$cluster, jmerged.log2fc, stringsAsFactors = FALSE)
# jmerged.zscore <- data.frame(bin = jmerged$bin, cluster = jmerged$cluster, jmerged.zscore, stringsAsFactors = FALSE)
# 
# # plot the 2D
# ggplot(jmerged.log2fc, aes(x = H3K4me3, y = H3K27me3, color = cluster)) + 
#   geom_point_rast(size = jsize, alpha = 0.3) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~cluster)
# 
# ggplot(jmerged.log2fc, aes(x = H3K4me3, color = cluster)) + 
#   geom_point_rast(size = jsize, alpha = 0.3) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~cluster)


# Wrap up and save objects ------------------------------------------------




if (make.plots){
  dev.off()
}

# write objects
save(clstrs.filt.dat, dat.annots.filt, tss.mats.filt, tss.mats.pbulk.ds, jlong.lst, jmerged, jlong.diff.genesets, de.ens.lst, file = outrdata)


print(Sys.time() - jstart)


