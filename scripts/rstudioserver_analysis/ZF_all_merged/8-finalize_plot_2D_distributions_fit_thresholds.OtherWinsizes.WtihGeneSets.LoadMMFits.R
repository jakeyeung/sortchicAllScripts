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


# Functions ---------------------------------------------------------------


CalculateGeometricMedian <- function(x1, x2, x3, cname.out, cnames = c("H3K4me1", "H3K4me3", "H3K27me3")){
  # create matrix
  X <- as.matrix(data.frame(x1 = x1, x2 = x2, x3 = x3, stringsAsFactors = FALSE))
  colnames(X) <- cnames
  gmed.out <- Gmedian::Gmedian(X)
  colnames(gmed.out) <- cnames
  gmed.out <- as.data.frame(gmed.out)
  return(gmed.out[[cname.out]])
}


# Constants ---------------------------------------------------------------


make.plots <- TRUE

jalpha <- 0.2
jdotsize <- 10

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

gene.annot <- subset(jlong.diff.genesets.WKM, select = c(bin, gene, ens, geneset))
gene.annot <- gene.annot[!duplicated(gene.annot), ]


# Plot 2D with different gene sets ----------------------------------------

# make mat 
head(jlong.thres.lst$H3K4me1)
head(jlong.thres.lst$H3K4me3)
head(jlong.thres.lst$H3K27me3)

# add zscore
jlong.thres.lst <- lapply(jlong.thres.lst, function(jdat){
  jdat <- jdat %>%
    group_by(cluster, mark) %>%
    mutate(zscore.counts = scale(counts, center = TRUE, scale = TRUE))
})

jmerged.s2n.fc <- reshape2::dcast(as.data.frame(jlong.thres.lst %>% bind_rows()), formula = "bin + cluster ~ mark", value.var = "s2n.fc") %>%
  left_join(., gene.annot)
jmerged.s2n.sqrt <- reshape2::dcast(as.data.frame(jlong.thres.lst %>% bind_rows()), formula = "bin + cluster ~ mark", value.var = "s2n.diffsqrt") %>%
  left_join(., gene.annot)
jmerged.s2n.lin <- reshape2::dcast(as.data.frame(jlong.thres.lst %>% bind_rows()), formula = "bin + cluster ~ mark", value.var = "s2n.diff") %>%
  left_join(., gene.annot)
jmerged.s2n.zscore <- reshape2::dcast(as.data.frame(jlong.thres.lst %>% bind_rows()), formula = "bin + cluster ~ mark", value.var = "zscore.counts") %>%
  left_join(., gene.annot)

jmerged.lst <- list(jmerged.s2n.fc, jmerged.s2n.sqrt, jmerged.s2n.lin, jmerged.s2n.zscore)
names(jmerged.lst) <- c("Signal2Noise log2FoldChange", "Sqrt Diff", "Lin Diff", "Counts Zscore")
jnames <- names(jmerged.lst); names(jnames) <- jnames

m.2d.plots.all.h3k4me1 <- lapply(jnames, function(jname){
  m <- ggplot(jmerged.lst[[jname]], aes(x = H3K4me1, y = H3K27me3, color = cluster)) + 
    ggrastr::geom_point_rast(size = jdotsize, alpha = jalpha) + 
    facet_wrap(~cluster, ncol = 2) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
    geom_density2d(color = "black") + 
    ggtitle(jname)
  return(m)
})

m.2d.plots.all.h3k4me3 <- lapply(jnames, function(jname){
  m <- ggplot(jmerged.lst[[jname]], aes(x = H3K4me3, y = H3K27me3, color = cluster)) + 
    ggrastr::geom_point_rast(size = jdotsize, alpha = jalpha) + 
    facet_wrap(~cluster, ncol = 2) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
    geom_density2d(color = "black") + 
    ggtitle(jname)
  return(m)
})

print(m.2d.plots.all.h3k4me1)
print(m.2d.plots.all.h3k4me3)


# Mark gene sets  ---------------------------------------------------------

m.2d.plots.all.h3k4me1.geneset <- lapply(jnames, function(jname){
  m <- ggplot(jmerged.lst[[jname]], aes(x = H3K4me1, y = H3K27me3, color = cluster)) + 
    ggrastr::geom_point_rast(size = jdotsize, alpha = jalpha) + 
    # facet_wrap(~cluster, ncol = 2) + 
    facet_grid(geneset~cluster) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
    geom_density2d(color = "black") + 
    ggtitle(jname)
  return(m)
})
print(m.2d.plots.all.h3k4me1.geneset$`Signal2Noise log2FoldChange`)

m.2d.plots.all.h3k4me3.geneset <- lapply(jnames, function(jname){
  m <- ggplot(jmerged.lst[[jname]], aes(x = H3K4me3, y = H3K27me3, color = cluster)) + 
    ggrastr::geom_point_rast(size = jdotsize, alpha = jalpha) + 
    facet_grid(geneset~cluster) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
    geom_density2d(color = "black") + 
    ggtitle(jname)
  return(m)
})
print(m.2d.plots.all.h3k4me3.geneset$`Signal2Noise log2FoldChange`)


# Add arrows  -------------------------------------------------------------


jmeth <- "Mean"
jmeth <- "GeometricMedian"
jmerged.sum.lst <- lapply(jmerged.lst, function(jmerged){
  if (jmeth == "GeometricMedian"){
    jmerged.sum.tmp <- jmerged %>%
      group_by(geneset, cluster) %>%
      summarise(H3K4me3 = CalculateGeometricMedian(x1 = H3K4me1, x2 = H3K4me3, x3 = H3K27me3, cname.out = "H3K4me3", cnames = c("H3K4me1", "H3K4me3", "H3K27me3")),
                H3K27me3 = CalculateGeometricMedian(x1 = H3K4me1, x2 = H3K4me3, x3 = H3K27me3, cname.out = "H3K27me3", cnames = c("H3K4me1", "H3K4me3", "H3K27me3")),
                H3K4me1 = CalculateGeometricMedian(x1 = H3K4me1, x2 = H3K4me3, x3 = H3K27me3, cname.out = "H3K4me1", cnames = c("H3K4me1", "H3K4me3", "H3K27me3")),
                H3K4me1.start = 0,
                H3K4me3.start = 0,
                H3K27me3.start = 0)
  } else if (jmeth == "Mean") {
    jmerged.sum.tmp <- jmerged %>%
      group_by(geneset, cluster) %>%
      summarise(H3K4me3 = mean(H3K4me1),
                H3K27me3 = mean(H3K27me3),
                H3K4me1 = mean(H3K4me1),
                H3K4me1.start = 0,
                H3K4me3.start = 0,
                H3K27me3.start = 0)
  } else {
    stop(jmeth, "not yet coded")
  }
})

# make 2D plot with arrows
m.2d.witharrows.h3k4me3 <- lapply(jnames, function(jname){
  m <- ggplot(mapping = aes(x = H3K4me3, y = H3K27me3, color = geneset)) + 
    geom_point(data = jmerged.lst[[jname]], mapping = aes(x = H3K4me3, y = H3K27me3), alpha = 0.1) + 
    geom_segment(data = jmerged.sum.lst[[jname]], mapping = aes(xend = H3K4me3.start, yend = H3K27me3.start), arrow = arrow(length=unit(0.10,"cm"), ends = "first"), color = "black") + 
    geom_vline(xintercept = 0, linetype = "dotted") + geom_hline(yintercept = 0, linetype = "dotted") +  
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_grid(geneset ~ cluster) + ggtitle(paste(jname, "arrow points to:", jmeth))
  return(m)
})
# m.2d.witharrows.h3k4me3$`Signal2Noise log2FoldChange`
# m.2d.witharrows.h3k4me3$`Sqrt Diff`
# m.2d.witharrows.h3k4me3$`Lin Diff`
# m.2d.witharrows.h3k4me3$`Counts Zscore`

# jmerged.sum.lst$`Counts Zscore`

m.2d.witharrows.h3k4me1 <- lapply(jnames, function(jname){
  m <- ggplot(mapping = aes(x = H3K4me1, y = H3K27me3, color = geneset)) + 
    geom_point(data = jmerged.lst[[jname]], mapping = aes(x = H3K4me3, y = H3K27me3), alpha = 0.1) + 
    geom_segment(data = jmerged.sum.lst[[jname]], mapping = aes(xend = H3K4me3.start, yend = H3K27me3.start), arrow = arrow(length=unit(0.10,"cm"), ends = "first"), color = "black") + 
    geom_vline(xintercept = 0, linetype = "dotted") + geom_hline(yintercept = 0, linetype = "dotted") +  
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_grid(geneset ~ cluster) + ggtitle(paste(jname, jmeth))
  return(m)
})
# m.2d.witharrows.h3k4me1$`Signal2Noise log2FoldChange`
# m.2d.witharrows.h3k4me1$`Counts Zscore`

print(m.2d.witharrows.h3k4me3)
print(m.2d.witharrows.h3k4me1)


if (make.plots){
  dev.off()
}

# Write matrices to output ------------------------------------------------

lapply(jnames, function(jname){
  print(jname)
  jname.str <- gsub(" ", "_", jname)
  outname.jname <- paste0("WKM_2d_matrix_values", Sys.Date(), ".winsize_", jwinsize, ".value_", jname.str, ".txt")
  outf.jname <- file.path(outdir, outname.jname)
  fwrite(jmerged.lst[[jname]], file = outf.jname, sep = "\t")
})


