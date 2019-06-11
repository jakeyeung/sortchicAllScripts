# Jake Yeung
# Date of Creation: 2019-06-08
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/filter_stringent_H3K4me3/plot_UMAP_with_FACS.R
# Plot UMAP with FACS to show granulocyte trajectory 


library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(JFuncs)

source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")


# Load H3K4me3  -----------------------------------------------------------

inf.dat.stringent <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"
assertthat::assert_that(file.exists(inf.dat.stringent))
dat.umap.long.stringent <- LoadUmap(inf.dat.stringent)

dat.umap.long.stringent$mark <- "H3K4me3"



# Load umaps  -------------------------------------------------------------

jmarks.noh3k4me3 <- c("H3K4me1", "H3K27me3", "H3K9me3"); names(jmarks.noh3k4me3) <- jmarks.noh3k4me3

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks



inmain <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
inf.dats <- lapply(jmarks.noh3k4me3, function(jmark) file.path(inmain, paste0("dat_umap_long_with_louvain.", jmark, ".RData")))

dat.umap.long <- lapply(inf.dats, LoadUmap)

dat.umap.long$H3K4me3 <- dat.umap.long.stringent

dat.umap.long <- dat.umap.long %>%
  bind_rows()


# Plot umaps for sanity ---------------------------------------------------

PlotXYNoColor(dat.umap.long %>% filter(mark == "H3K4me1"), xvar = "umap1", yvar = "umap2")
PlotXYNoColor(dat.umap.long %>% filter(mark == "H3K4me3"), xvar = "umap1", yvar = "umap2")
PlotXYNoColor(dat.umap.long %>% filter(mark == "H3K27me3"), xvar = "umap1", yvar = "umap2")
PlotXYNoColor(dat.umap.long %>% filter(mark == "H3K9me3"), xvar = "umap1", yvar = "umap2")


# Check the FACS data -----------------------------------------------------

dat.facs.filt <- lapply(jmarks, LoadFACSGetLoadings) %>%
  bind_rows()

# Do the PCA plot  --------------------------------------------------------

dat.facs.filt.lst <- lapply(jmarks, LoadFACSGetLoadings)

# bind to umap long
dat.merge <- left_join(dat.umap.long, dat.facs.filt %>% dplyr::select(cell, loadings, mark)) %>%
  filter(!is.na(loadings))

pdf(paste0("/Users/yeung/data/scchic/pdfs/B6_figures/stringent_pdfs/FACS_umaps/FACS_umaps_all_marks_stringent.", Sys.Date(), ".pdf"), useDingbats = FALSE)
for (jmark in jmarks){
  jsub <- dat.facs.filt.lst[[jmark]]
  inmat <- jsub %>%
    dplyr::select(-cell, -loadings, -mark)
  # rownames(inmat) <- jsub$cell
  inmat <- scale(as.matrix(inmat), center = TRUE, scale = TRUE)
  pca.out <- prcomp(inmat, center = FALSE, scale. = FALSE)
  if (jmark == "H3K4me1"){
    print("Switching rotation space for H3K9me3")
    pca.out$rotation[, 1] <- pca.out$rotation[, 1] * -1
    pca.out$x[, 1] <- pca.out$x[, 1] * -1
  } 
  biplot(pca.out, xlabs=rep("o", nrow(inmat)), cex = 1, main = jmark, scale = 0, cex.main = 2.5, cex.axis = 1, cex.lab = 2.5)
  # plot just the PC1 and PC2 without the arrows 
  dat.proj <- scale(as.matrix(inmat) %*% pca.out$rotation, center = TRUE, scale = TRUE)
  dat.proj.long <- data.frame(PC1 = dat.proj[, 1] * pca.out$sdev[1], 
                              PC2 = dat.proj[, 2] * pca.out$sdev[2], 
                              PC1.Loadings = pca.out$x[, 1])
  m1 <- PlotXYWithColor(dat.proj.long, xvar = "PC1", yvar = "PC2", cname = "PC1.Loadings", remove.axis.info = FALSE, jtitle = jmark, jsize = 6)
  m1 <- m1 + geom_vline(xintercept = 0, color = "black", linetype = "dotted") + geom_hline(yintercept = 0, color = "black", linetype = "dotted") + xlab("PC1") + ylab("PC2")
  # m1 <- ggplot(dat.proj.long, aes(x = PC1, y = PC2, color = PC1.Loadings)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m1)
}

mlst <- lapply(jmarks, function(jmark){
  m <- PlotXYWithColor(dat.merge %>% filter(mark == jmark), xvar = "umap1", yvar = "umap2", cname = "loadings", jtitle = jmark, jsize = 4) 
  # print(m)
})
multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 4)
multiplot(mlst[[1]], mlst[[3]], mlst[[2]], mlst[[4]], cols = 2)
dev.off()