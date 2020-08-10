# Jake Yeung
# Date of Creation: 2020-06-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/redo_multiomics_integration_TSS_pseudobulk_stringent_DE.set_up_tables_for_poisson_regression.R
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
library(reshape2)

make.plots <- TRUE
save.long <- FALSE   # faster if false

# Load DE genes -----------------------------------------------------------

# load this first because it loads a lot of objects, might disuprt things

inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData"
load(inf.de, v=T)

# Load TSS  ---------------------------------------------------------------

fewer.k27me3 <- TRUE
downsample.reads <- TRUE
downsample.cells <- FALSE

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression"
dir.create(outdir, showWarnings = TRUE)

outprefix <- file.path(outdir, paste0("integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_", fewer.k27me3, ".forPoissonRegression.CountR1only.", Sys.Date()))

outpdf <- paste0(outprefix, ".pdf")
if (!save.long){
  outrdata <- paste0(outprefix, ".smaller.RData")
} else {
  outrdata <- paste0(outprefix, ".RData")
}

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

print(m.umaps.beforeafter)


# Filter TSS for cells of relevant cluster  ---------------------------------------------------

lapply(dat.annots.filt, function(jdat){
  sort(unique(jdat$cluster.new))
})

cells.keep.lst <- lapply(dat.annots.filt, function(jdat){
  jdat$cell
})

tss.mats.filt.fromref.cellfilt <- lapply(jmarks, function(jmark){
  jmat <- tss.mats.filt.fromref[[jmark]]
  jmat[, cells.keep.lst[[jmark]]]
})

# count mat bug: round up to nearest even, then divide by 2. Must be all integers
tss.mats.filt.fromref.cellfilt <- lapply(tss.mats.filt.fromref.cellfilt, function(jmat){
  jmat.fixed <- 2 * ceiling(jmat / 2)
  jmat.fixed <- jmat.fixed / 2
  return(jmat.fixed)
})


lapply(tss.mats.filt.fromref.cellfilt, dim)

ncuts.cells <- lapply(tss.mats.filt.fromref.cellfilt, function(jmat){
  ncuts.dat <- data.frame(cell = colnames(jmat), ncuts.total = colSums(jmat), stringsAsFactors = FALSE)
})

# Make long object --------------------------------------------------------


if (save.long){
  jlong.cellfilt.counts <- lapply(jmarks, function(jmark){
    jmat <- tss.mats.filt.fromref.cellfilt[[jmark]]
    # set up the data frame for melting
    jdat <- data.frame(bin = rownames(jmat), as.data.frame(as.matrix(jmat)), stringsAsFactors = FALSE)
    colnames(jdat) <- c("bin", colnames(jmat))
    # annotate bins
    jdat$mark <- jmark
    jdat$gene <- sapply(as.character(jdat$bin), function(b) strsplit(b, ";")[[1]][[2]])
    jdat$ens <- sapply(jdat$gene, function(g) AssignHash(g, g2e, null.fill = g))
    jlong <- reshape2::melt(jdat, id.vars = c("bin", "mark", "gene", "ens"), variable.name = "cell", value.name = "ncuts")
    jlong <- left_join(jlong, ncuts.cells[[jmark]])  # add total counts
    return(jlong)
  })
}

# Wrap up and save objects ------------------------------------------------


if (make.plots){
  dev.off()
}

# write objects
if (!save.long){
  save(clstrs.filt.dat, dat.annots.filt, tss.mats.filt.fromref.cellfilt, de.ens.sorted.stringent, g2e, ncuts.cells, file = outrdata)
} else {
  save(clstrs.filt.dat, dat.annots.filt, jlong.cellfilt.counts, de.ens.sorted.stringent, g2e, ncuts.cells, file = outrdata)
}

print(Sys.time() - jstart)


