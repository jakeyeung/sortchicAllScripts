# Jake Yeung
# Date of Creation: 2020-06-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/downstream_poisson_model_analysis.ZF.R
# 



rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(ggrastr)
library(scchicFuncs)
library(JFuncs)


# Functions ---------------------------------------------------------------



# Constants ---------------------------------------------------------------

ncores <- 16

make.plots <- TRUE
save.objs <- TRUE

jdate <- "2020-06-08"
jdate2 <- "2020-06-09"

# gset.names.new[which(gset.names.new == "lymphocytes")] <- "Bcell"
# gset.names.new[which(gset.names.new == "granulocytes")] <- "Neutrophil"
# gset.names.new[which(gset.names.new == "erythrocytes")] <- "Erythro"

# gsets.filt <- c("lymphocytes", "erythrocytes", "HSPCs", "granulocytes")
gsets.filt <- c("Bcell", "Erythro", "Neutrophil")
# jgenes <- c("meis1b", "tal1", "pax5", "s100z")
jgenes <- c("meis1b;", "myb;", "pmp22b;", "tal1;", "gata1a;", "pax5;", "cd79a;", "lyz;")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/setup_matrix_for_poisson_regression.likeBM.redo_count_tables"
assertthat::assert_that(dir.exists(indir))

infrdata <- file.path(indir, paste0("integrated_analysis.", jdate2, ".UseTSSfromH3K4me3.likeBM.RData"))

datdir <- "/home/jyeung/data/from_rstudioserver/poisson_fits.redo_count_tables"
dir.create(datdir)

outfits <- file.path(indir, paste0("fit_poisson_model_on_TSS_ZF.", jdate2, ".RData"))
outfits.ci <- file.path(datdir, paste0("fit_poisson_model_on_TSS_ZF.DownstreamWrangled.", Sys.Date(), ".ClusterRenamed.ConfidenceIntervals.smaller.RData"))
# outpdf <- file.path(indir, paste0("fit_poisson_model_on_TSS_ZF.DownstreamWrangled.", Sys.Date(), ".ClusterRenamed.pdf"))

# infrdata <- paste0(jprefix, ".smaller.RData")

assertthat::assert_that(file.exists(infrdata))

load(infrdata, v=T)
load(outfits, v=T)


# Wrangle -----------------------------------------------------------------

# create g2e object
genes.vec <- sapply(rownames(tss.mats.sc.filt.zf$H3K4me1), function(x) strsplit(x, ";")[[1]][[2]])
ens.vec <- Gene2Ensembl.ZF(genes.vec, return.original = TRUE, species = "drerio")

g2e <- hash::hash(genes.vec, ens.vec)
e2g <- hash::invert(g2e)


jfits.dat.lst <- lapply(jmarks, function(jmark){
  jfit.dat <- do.call(rbind, jfits.lst.bymark[[jmark]])
  jfit.dat$bin <- rownames(jfit.dat)
  jfit.dat$gene <- sapply(as.character(jfit.dat$bin), function(b) strsplit(b, ";")[[1]][[2]])
  jfit.dat$ens <- sapply(jfit.dat$gene, function(g) AssignHash(g, g2e, null.fill = g))
  return(jfit.dat)
})

jfits.long.lst <- lapply(jmarks, function(jmark){
  jlong <- jfits.dat.lst[[jmark]] %>%
    dplyr::select(-c(dev.diff, df.diff)) %>%
    reshape2::melt(., id.vars = c("bin", "gene", "ens", "pval"), variable.name = "cluster", value.name = "logLambda")
  # mutate(cluster = ifelse(cluster == "X.Intercept.", "ClusterHSPCs", as.character(cluster)))
  jlong$mark <- jmark
  # extract the X intercept as separate column 
  jlong.noint <- subset(jlong, cluster != "X.Intercept.")
  jlong.int <- subset(jlong, cluster == "X.Intercept.")  %>%
    dplyr::rename(logintercept = logLambda) %>%
    dplyr::select(bin, logintercept, mark)
  jlong.merge <- left_join(jlong.noint, jlong.int)
  jlong.merge$cluster <- as.character(jlong.merge$cluster)
  return(jlong.merge)
})



# Rename clusters to match BM ---------------------------------------------

unique(jfits.long.lst$H3K4me1$cluster)

cnames.orig <- c("Clustereryth", "Clustergranu", "Clusterlymph")
cnames.new <- c("ClusterErythroblasts", "ClusterGranulocytes", "ClusterBcells")
cnames.hash <- hash::hash(cnames.orig, cnames.new)

jfits.long.lst <- lapply(jfits.long.lst, function(jfits){
  jfits$cluster <- sapply(jfits$cluster, AssignHash, cnames.hash)
  return(jfits)
})

 
unique(jfits.long.lst$H3K4me1$cluster)
# unique(jfits.long.test.lst$H3K4me1$cluster) 


# Rename genesets ---------------------------------------------------------

gset.names.old <- names(genesets)
gset.names.new <- gset.names.old
gset.names.new[which(gset.names.new == "lymphocytes")] <- "Bcell"
gset.names.new[which(gset.names.new == "granulocytes")] <- "Neutrophil"
gset.names.new[which(gset.names.new == "erythrocytes")] <- "Erythro"
gset.names.hash <- hash::hash(gset.names.old, gset.names.new)

de.ens.sorted.stringent <- genesets
gset.names <- names(de.ens.sorted.stringent)
names(gset.names) <- gset.names.new
names(de.ens.sorted.stringent) <- gset.names.new
gset.names <- gset.names.new



# Run refits with CI ------------------------------------------------------

print(paste("Running fits CI, ncores:", ncores))

system.time(
  jfits.ci.lst <- lapply(jmarks, function(jmark){
    jrow.names <- rownames(tss.mats.sc.filt.zf[[jmark]])
    names(jrow.names) <- jrow.names
    cnames <- colnames(tss.mats.sc.filt.zf[[jmark]])
    dat.annots.filt.mark <- dat.annots.filt.forfit[[jmark]]
    ncuts.cells.mark <- ncuts.cells[[jmark]]
    refits.lst <- parallel::mclapply(jrow.names, function(jrow.name){
      refit <- RefitPoissonForPlot(jrow = tss.mats.sc.filt.zf[[jmark]][jrow.name, ], cnames = cnames, dat.annots.filt.mark = dat.annots.filt.mark, ncuts.cells.mark = ncuts.cells.mark, return.means = FALSE)
      return(refit)
    }, mc.cores = ncores)
  })
)


if (save.objs){
    save(jfits.ci.lst, file = outfits.ci)
}


