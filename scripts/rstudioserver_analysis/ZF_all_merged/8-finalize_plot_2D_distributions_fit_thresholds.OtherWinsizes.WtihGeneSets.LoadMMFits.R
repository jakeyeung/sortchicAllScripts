# Jake Yeung
# Date of Creation: 2020-05-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/8-finalize_plot_2D_distributions_fit_thresholds.OtherWinsizes.WtihGeneSets.R
# description

rm(list=ls())

library(ggrastr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(preprocessCore)

library(mixtools)

library(JFuncs)

jalpha <- 0.15
jdotsize <- 5

# try the hardest one first, K27me3 sqrt for Bcells?
jlambda <- c(0.5, 0.5)
jmu <- c(1, 3)
jk <- 2
jsigma <- c(1, 1)
jthres <- 0.5
mincounts <- 4

jwinsize <- 10000L
jdate <- "2020-05-11"

jpos <- "bottom"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/pseudobulk_analysis.integrated_analysis.TSS.plot_2D_distributions.OtherWinsizes"
dir.create(outdir)
outname <- paste0("WKM_2D_active_vs_repressed_fit_thresholds.", Sys.Date(), ".winsize_", jwinsize, ".WithGeneSets.Outputs.pdf")
outname.rdata <- paste0("WKM_2D_active_vs_repressed_fit_thresholds.", jdate, ".winsize_", jwinsize, ".WithGeneSets.RData")

outf <- file.path(outdir, outname)
outf.rdata <- file.path(outdir, outname.rdata)

make.plots <- TRUE
if (make.plots){
  pdf(outf, useDingbats = FALSE)
}



# Load 2D pseudobulks -----------------------------------------------------


load(outf.rdata, v=T)

# zebrafish WKM

jwinsize <- 10000L
downsample.reads <- TRUE
downsample.cells <- FALSE
jdate <- "2020-05-09"  # downsample.cells FALSE now does not down sample 
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/pseudobulk_analysis.integrated_analysis.TSS.readsDS.likeBM.OtherWinsizes"

# inprefix <- "integrated_analysis.2020-05-09.UseTSSfromH3K4me3.DSreads_TRUE.DScells_FALSE.likeBM.winsize_10000.RData"
inprefix <- paste0("integrated_analysis.", jdate, ".UseTSSfromH3K4me3.DSreads_", downsample.reads, ".DScells_", downsample.cells, ".likeBM.winsize_", jwinsize)
inrdata <- paste0(inprefix, ".RData")
infrdata <- file.path(indir, inrdata)

assertthat::assert_that(file.exists(infrdata))
load(infrdata, v=T)
jlong.diff.genesets.WKM <- jlong.diff.genesets
jlong.diff.genesets.WKM$cluster <- sapply(jlong.diff.genesets.WKM$cluster, function(x) gsub("monocyte", "granu", x = x))
jlong.diff.genesets.WKM$cluster <- sapply(jlong.diff.genesets.WKM$cluster, function(x) gsub("HSC", "HSPCs", x = x))
jlong.diff.genesets.WKM$cluster <- factor(jlong.diff.genesets.WKM$cluster, levels = c("HSPCs", "granu", "lymph", "eryth"))
jlong.diff.genesets.WKM$mark <- factor(jlong.diff.genesets.WKM$mark, jmarks)
jlong.diff.genesets.WKM$geneset <- factor(jlong.diff.genesets.WKM$geneset, c("HSPCs", "granulocytes", "lymphocytes", "erythrocytes", "HighExprs", "LowExprs", "zOther"))
gsetfilt.WKM <- c("HSPCs", "granulocytes", "lymphocytes", "erythrocytes")


# Plot 2D with different gene sets ----------------------------------------

# make mat 

head(jlong.thres.lst$H3K4me1)

