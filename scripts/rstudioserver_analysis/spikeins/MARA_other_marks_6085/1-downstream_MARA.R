# Jake Yeung
# Date of Creation: 2021-02-18
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/MARA_other_marks_6085/1-downstream_MARA.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(JFuncs)
library(scchicFuncs)
library(ggrepel)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"

outdir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/MARA_output_BM_k9dynamicbins")

zscores.cutoff <- 1.25
jsize <- 1.5
outpdf <- file.path(outdir, paste0("H3K27me3_MARA_outputs.zscorecutoff_", zscores.cutoff, ".dotsize_", jsize, "." , Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)



# Load UMAP  --------------------------------------------------------------

indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage")

dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  fread(file.path(indir.meta, fname))
}) 


# Load MARA  --------------------------------------------------------------

jmark <- "H3K27me3"
indir.mara <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_dynamic_bins_top_6085/", jmark, "/mara_output/ldaOut.count_tables_merged.", jmark, ".DE_bins_all_marks_top_6085_dists_to_TSS.annot_table.", jmark, ".2021-02-15.txt.K-30.keepNbins_0.withchr2-dynamic_bins_50kb__motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, "/ldaOut.count_tables_merged.", jmark, ".DE_bins_all_marks_top_6085_dists_to_TSS.annot_table.", jmark, ".2021-02-15.txt.K-30.keepNbins_0.withchr2"))
print(indir.mara)
assertthat::assert_that(dir.exists(indir.mara))

mara.out <- LoadMARA(indir.mara, make.cnames = FALSE)

act.mat.clean <- t(subset(mara.out$act.mat, select = -motif))
colnames(act.mat.clean) <- mara.out$act.mat$motif
act.mat.clean.dat <- data.frame(cell = rownames(act.mat.clean), act.mat.clean, stringsAsFactors = FALSE) %>% 
  ungroup() %>%
  mutate(cell = gsub("\\.", "-", cell))

motifs.keep <- subset(mara.out$zscores, zscore > zscores.cutoff)$motif

m.zscores <- mara.out$zscores %>%
  ungroup() %>% 
  mutate(rnk = seq(length(motif)),
         motiflab = ifelse(zscore > zscores.cutoff, toupper(motif), NA)) %>%
  ggplot(., aes(x = rnk, y = zscore, label = motiflab)) + 
  geom_point() + 
  geom_text_repel() + 
  theme_bw() + 
  geom_hline(yintercept = zscores.cutoff, linetype = "dotted") + 
  xlab("Rank") + 
  ylab("Motif zscore") + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.zscores)

m.zscores2 <- mara.out$zscores %>%
  ungroup() %>% 
  mutate(rnk = seq(length(motif))) %>%
  ggplot(., aes(x = zscore)) + 
  geom_density(fill = "red", alpha = 0.25) + 
  geom_vline(xintercept = zscores.cutoff, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.zscores2)


# add motifs to UMAP  -----------------------------------------------------

dat.merge.motifs <- left_join(dat.metas[[jmark]], act.mat.clean.dat, by = "cell")

print(head(motifs.keep))

jmotif <- "Zbtb16"
jmotif <- "Hic1"
jmotif <- "Gata3"
jmotif <- "Wt1"
jmotif <- "Ebf1"
jmotif <- "Tfdp1"
jmotif <- "Pax8"
jmotif <- "Yy2"
jmotif <- "Yy1"
jzscore <- signif(subset(mara.out$zscores, motif == jmotif)$zscore, digits = 2)

for (jmotif in motifs.keep){
# pdf("/home/jyeung/hub_oudenaarden/jyeung/tmp/motiftest.pdf", useDingbats = FALSE)
  (jtitle <- paste(jmotif, "Zscore:", jzscore))
  m <- ggplot(dat.merge.motifs, aes_string(x = "umap1", y = "umap2", color = jmotif)) + 
    geom_point(size = jsize) + 
    ggtitle(jtitle) + 
    theme_bw() + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  # m <- PlotXYWithColor(dat.merge.motifs, xvar = "umap1", yvar = "umap2", cname = jmotif, jtitle = jtitle, cont.color = TRUE, jsize = 0.5) + scale_color_viridis_c()
  print(m)
}


# Add heatmap -------------------------------------------------------------

jmat <- as.matrix(subset(act.mat.clean.dat, select = -cell))
rownames(jmat) <- act.mat.clean.dat$cell

# jmat[dat.metas$H3K9me3$cell, ][1:5, 1:5]

cells.ordered <- dat.metas[[jmark]]$cell
jsub <- jmat[cells.ordered, motifs.keep]


library(heatmap3)

colvec <- dat.metas[[jmark]]$clustercol

jmeth <- "ward.D2"

par(mfrow=c(1,1), mar=c(1,1,1,1), mgp=c(3, 1, 0), las=0)
hm.out <- heatmap3(jsub, margins = c(5, 8), cexCol = 0.35, Colv = TRUE, Rowv = NA, 
                   main = jmark, 
                   # ColSideColors = rep("blue", ncol(jsub)), 
                   # ColSideColors = FALSE,
                   RowSideColors = colvec, 
                   # RowSideColors = rep("red", nrow(jsub)), 
                   RowSideLabs = "celltype", 
                   labRow = FALSE, scale = "column", revC = TRUE,
                   distfun = dist, hclustfun = hclust, method = jmeth)

hm.out <- heatmap3(jsub, margins = c(5, 8), cexCol = 0.35, Colv = TRUE, Rowv = NA, 
                   main = jmark, 
                   # ColSideColors = rep("blue", ncol(jsub)), 
                   # ColSideColors = FALSE,
                   RowSideColors = colvec, 
                   # RowSideColors = rep("red", nrow(jsub)), 
                   RowSideLabs = "celltype", 
                   labRow = FALSE, scale = "column", revC = FALSE,
                   distfun = dist, hclustfun = hclust, method = jmeth)


hm.out.transpose <- heatmap3(t(jsub), margins = c(5, 8), cexCol = 0.35, Colv = NA, Rowv = TRUE, 
                             main = jmark, 
                             # ColSideColors = rep("blue", ncol(jsub)), 
                             # ColSideColors = FALSE,
                             ColSideColors = colvec, 
                             # RowSideColors = rep("red", nrow(jsub)), 
                             ColSideLabs = "celltype", 
                             labCol = FALSE, scale = "row", revC = FALSE, 
                             distfun = dist, hclustfun = hclust, method = jmeth)

hm.out.transpose <- heatmap3(t(jsub), margins = c(5, 8), cexCol = 0.35, Colv = NA, Rowv = TRUE, 
                             main = jmark, 
                             # ColSideColors = rep("blue", ncol(jsub)), 
                             # ColSideColors = FALSE,
                             ColSideColors = colvec, 
                             # RowSideColors = rep("red", nrow(jsub)), 
                             ColSideLabs = "celltype", 
                             labCol = FALSE, scale = "row", revC = TRUE,
                             distfun = dist, hclustfun = hclust, method = jmeth)




dev.off()



# Plot Zscores ------------------------------------------------------------




# Plot activities ---------------------------------------------------------





