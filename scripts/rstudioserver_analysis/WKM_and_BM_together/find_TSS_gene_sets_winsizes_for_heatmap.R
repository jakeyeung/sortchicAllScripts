# Jake Yeung
# Date of Creation: 2020-06-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/find_TSS_gene_sets_winsizes_for_heatmap.R
# Get winsizes for heatmap (see spatial distribution)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(JFuncs)

outdir.bm <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations/MouseBM"
outdir.zf <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations/ZebrafishWKM"

# Constants ---------------------------------------------------------------

# bsizes <- c(6000, 10000, 20000, 50000)
bsizes <- c(2)

# Load genesets -----------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

jdate <- "2020-06-06"
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression"
infit.wrangled <- file.path(indir, paste0("fit_poisson_model_on_TSS.DownstreamWrangled.", jdate, ".RData"))
load(infit.wrangled, v=T)

fits.bygenesets.long.bm <- fits.bygenesets.long

gsets.bm <- unique(fits.bygenesets.long.bm$geneset)
names(gsets.bm) <- gsets.bm


# remove IGK genes from Bcells ---------------------------------------------------

fits.bygenesets.long.filt.bm <- subset(fits.bygenesets.long.bm, !grepl("Igk", gene, ignore.case = TRUE))

# Wrangle  ----------------------------------------------------------------

tss.bygeneset.bm <- lapply(gsets.bm, function(gset.bm){
  print(gset.bm)
  jbins <- unique(subset(fits.bygenesets.long.filt.bm, geneset == gset.bm)$bin)
  jcoords <- sapply(jbins, function(jbin) strsplit(jbin, ";")[[1]][[1]])
  jgenes <- sapply(jbins, function(jbin) strsplit(jbin, ";")[[1]][[2]])
  # get chromo and midpoint
  jmidpts <- sapply(jcoords, function(jcoord){
    jchromo <- JFuncs::GetChromo(jcoord)
    jstart <- as.numeric(JFuncs::GetStart(jcoord))
    jend <- as.numeric(JFuncs::GetEnd(jcoord))
    # get midpoint, nearest integer
    jmidpoint <- as.integer((jstart + jend) / 2)
    # jout <- paste(jchromo, jmidpoint, sep = ":")
    return(jmidpoint)
  })
  # jout2 <- paste(jchromo.midpt, jgene, sep = ';')
  jout.dat <- data.frame(chromo = sapply(jcoords, GetChromo), tss = jmidpts, gene = jgenes, stringsAsFactors = FALSE)
})


# Adjust different bin sizes  ---------------------------------------------

for (gset.bm in gsets){
  tss.tmp.bm <- tss.bygeneset.bm[[gset.bm]]
  for (bsize in bsizes){
    print(paste(gset.bm, bsize))
    outf.tmp.bm <- file.path(outdir.bm, paste0("MouseBM_TSS.", gset.bm, ".bsize_", bsize, ".bed"))
    hstep <- as.integer(bsize / 2)
    
    # write bedfile
    bed.tmp.bm <- tss.tmp.bm %>%
      rowwise() %>%
      mutate(Start = tss - hstep, 
             End = tss + hstep) %>%
      dplyr::select(chromo, Start, End, gene)
    # write outputs
    fwrite(bed.tmp.bm, file = outf.tmp.bm, sep = "\t", col.names = FALSE)
  }
}


# Do fr zebrafish  --------------------------------------------------------


jdate2 <- "2020-06-09"
indir.zf <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/setup_matrix_for_poisson_regression.likeBM.redo_count_tables"
assertthat::assert_that(dir.exists(indir.zf))
infits.wrangled.zf <- file.path(indir.zf, paste0("fit_poisson_model_on_TSS_ZF.DownstreamWrangled.",jdate2, ".ClusterRenamed.RData"))

load(infits.wrangled.zf, v=T)


fits.bygenesets.long.zf <- fits.bygenesets.long

gsets.zf <- unique(fits.bygenesets.long.zf$geneset)
names(gsets.zf) <- gsets.zf

fits.bygenesets.long.filt.zf <- subset(fits.bygenesets.long.zf, !grepl("Igk", gene, ignore.case = TRUE))


# Wrangle  ----------------------------------------------------------------

tss.bygeneset.zf <- lapply(gsets.zf, function(gset.zf){
  print(gset.zf)
  jbins <- unique(subset(fits.bygenesets.long.filt.zf, geneset == gset.zf)$bin)
  jcoords <- sapply(jbins, function(jbin) strsplit(jbin, ";")[[1]][[1]])
  jgenes <- sapply(jbins, function(jbin) strsplit(jbin, ";")[[1]][[2]])
  # get chromo and midpoint
  jmidpts <- sapply(jcoords, function(jcoord){
    jchromo <- JFuncs::GetChromo(jcoord)
    jstart <- as.numeric(JFuncs::GetStart(jcoord))
    jend <- as.numeric(JFuncs::GetEnd(jcoord))
    # get midpoint, nearest integer
    jmidpoint <- as.integer((jstart + jend) / 2)
    # jout <- paste(jchromo, jmidpoint, sep = ":")
    return(jmidpoint)
  })
  # jout2 <- paste(jchromo.midpt, jgene, sep = ';')
  jout.dat <- data.frame(chromo = sapply(jcoords, GetChromo), tss = jmidpts, gene = jgenes, stringsAsFactors = FALSE)
})




# Adjust different bin sizes  ---------------------------------------------

for (gset.zf in gsets.zf){
  tss.tmp.zf <- tss.bygeneset.zf[[gset.zf]]
  for (bsize in bsizes){
    print(paste(gset.zf, bsize))
    outf.tmp.zf <- file.path(outdir.zf, paste0("ZebrafishWKM_TSS.", gset.zf, ".bsize_", bsize, ".bed"))
    hstep <- as.integer(bsize / 2)
    
    # write bedfile
    bed.tmp.zf <- tss.tmp.zf %>%
      rowwise() %>%
      mutate(Start = tss - hstep, 
             End = tss + hstep) %>%
      dplyr::select(chromo, Start, End, gene)
    # write outputs
    fwrite(bed.tmp.zf, file = outf.tmp.zf, sep = "\t", col.names = FALSE)
  }
}

