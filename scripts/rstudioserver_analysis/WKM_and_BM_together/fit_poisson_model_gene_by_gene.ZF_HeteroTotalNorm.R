# Jake Yeung
# Date of Creation: 2020-06-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/fit_poisson_model_gene_by_gene.ZF_HeteroTotalNorm.R
# description


rm(list=ls())

jstart <- Sys.time()

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(hash)
library(JFuncs)
library(scchicFuncs)

ncores <- 16

hubprefix <- "/home/jyeung/hub_oudenaarden"

jdist <- "10000"
# jdist <- "50000"

jnorm <- "ByHetero"
# jnorm <- "ByTotalFromBins"

jspecies <- "ZebrafishWKM"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm"
outfits <- file.path(outdir, paste0("fit_poisson_model_on_TSS.", jspecies, ".NormMeth_", jnorm, ".RData"))

assertthat::assert_that(!file.exists(outfits))

inf.offsets <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/offsets_hetero_and_totalcuts/MouseBM_HeteroTotalCounts_50kb_bins.RData"

# Load annots -------------------------------------------------------------

wkm.rename <- hash(c("eryth1", "eryth2", "HSC1", "HSC2", "monocyte", "HSC"), c("eryth", "eryth", "HSPCs", "HSPCs", "granu", "HSPCs"))
dat.annot.lst.WKM <- lapply(jmarks, function(jmark){
  print(jmark)
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
    filter(cluster != "Unknown") %>%
    filter(cluster != "thrombo") %>%
    left_join(., subset(annot.louv, select = c(cell, var.imputed)))
  # rename clusters
  annot.glmpca.filt$cluster.new <- sapply(annot.glmpca.filt$cluster, function(jclst) AssignHash(jclst, wkm.rename, null.fill = jclst))
  annot.glmpca.filt$cluster.new <- sapply(annot.glmpca.filt$cluster.new, function(jclst) ifelse(jclst == "lymph", "Bcell", jclst))
  annot.glmpca.filt$cond <- sapply(annot.glmpca.filt$cell, ClipLast, jsep = "_")
  print("annot glmpca filt")
  print(annot.glmpca.filt)
  return(annot.glmpca.filt)
})

load(inf.offsets, v=T)

lapply(dat.annot.lst.WKM, function(jdat) print(unique(jdat$cluster.new)))

dat.annots.filt.forfit <- lapply(dat.annot.lst.WKM, function(jdat){
  jdat <- subset(jdat, select = c(cell, cluster.new)) %>%
    mutate(Cluster = ifelse(cluster.new == "HSPCs", "aHSPCs", cluster.new))  # set HSPC as intercept
  return(jdat)
})

if (jnorm == "ByHetero"){
  ncuts.for.fit <- lapply(dat.ncuts.hetero.total, function(jdat){
    # expects ncuts.total
    subset(jdat, select = c(cell, ncuts.hetero)) %>%
      dplyr::rename(ncuts.total = ncuts.hetero)
  })
} else if (jnorm == "ByTotalFromBins"){
  ncuts.for.fit <- lapply(dat.ncuts.hetero.total, function(jdat){
    # expects ncuts.total
    subset(jdat, select = c(cell, ncuts.total))
  })
} else {
  warning("jnorm must be ByHetero or ByTotalFromBins", jnorm)
}


# Load the new tables ------------------------------------------------

indir <- file.path(hubprefix, paste0("jyeung/data/zebrafish_scchic/count_tables.winsize_", jdist, ".imputevarfilt.lessstringent.mapq_40"))
tss.mats.filt.fromref.cellfilt <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0(jmark, ".imputevarfilt.lessstringent.mapq_40.remerged.countTable.csv")
  inf <- file.path(indir, fname)
  mat <- ReadMatSlideWinFormat(inf = inf, as.sparse = TRUE, sort.rnames = TRUE, add.chromo = FALSE)
  mat <- mat[rows.common, cells.keep[[jmark]]]
  return(mat)
})

lapply(mats.lst, dim)



print("Fitting... ")

system.time(
  jfits.lst.bymark <- lapply(jmarks, function(jmark){
  # jfits.lst.bymark <- parallel::mclapply(jmarks, function(jmark){
    print(jmark)
    jmat.mark <- tss.mats.filt.fromref.cellfilt[[jmark]]
    dat.annots.filt.mark <- dat.annots.filt.forfit[[jmark]]
    ncuts.for.fit.mark <- ncuts.for.fit[[jmark]]
    cnames <- colnames(jmat.mark)
    
    jrow.names <- rownames(jmat.mark)
    names(jrow.names) <- jrow.names
    jfits.lst <- parallel::mclapply(jrow.names, function(jrow.name){
      jrow <- jmat.mark[jrow.name, ]
      jout <- FitGlmRowClusters(jrow, cnames, dat.annots.filt.mark, ncuts.for.fit.mark, jbin = NULL, returnobj = FALSE)
      return(jout)
    }, mc.cores = ncores)
    
    return(jfits.lst)
  })
)

print("Saving objects")
save(tss.mats.filt.fromref.cellfilt, dat.annots.filt.forfit, ncuts.for.fit, jfits.lst.bymark, file = outfits)
print("Done Saving objects")

print(Sys.time() - jstart)



