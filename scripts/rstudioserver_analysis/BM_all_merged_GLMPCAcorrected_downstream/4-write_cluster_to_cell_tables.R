# Jake Yeung
# Date of Creation: 2020-02-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged_GLMPCAcorrected_downstream/4-write_cluster_to_cell_tables.R
# Make cluster to cell tables so we can create pseudobulk bams 
# output textfile should have column names "cell", "cluster" ... 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Load cluster labels  ----------------------------------------------------

outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/cell_cluster_tables"
assertthat::assert_that(dir.exists(outdir))


jmark1 <- "H3K4me1"
jmark2 <- "H3K4me3"
jmark3 <- "H3K27me3"
jmark4 <- "H3K9me3"
jmarks <- c(jmark1, jmark2, jmark3, jmark4)
names(jmarks) <- jmarks
inf1 <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark1, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
load(inf1, v=T)

dat.umap.glm.fillNAs1 <- dat.umap.glm.fillNAs
dat.umap.glm1 <- dat.umap.glm
dat.umap.lda1 <- dat.umap.lda
mm.celltype.lst1 <- mm.celltype.lst

inf2 <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark2, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
load(inf2, v=T)

dat.umap.glm.fillNAs2 <- dat.umap.glm.fillNAs
dat.umap.glm2 <- dat.umap.glm
dat.umap.lda2 <- dat.umap.lda
mm.celltype.lst2 <- mm.celltype.lst

inf3 <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark3, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
load(inf3, v=T)

dat.umap.glm.fillNAs3 <- dat.umap.glm.fillNAs
dat.umap.glm3 <- dat.umap.glm
dat.umap.lda3 <- dat.umap.lda
mm.celltype.lst3 <- mm.celltype.lst

inf4 <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark4, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
load(inf4, v=T)

dat.umap.glm.fillNAs4 <- dat.umap.glm.fillNAs
dat.umap.glm4 <- dat.umap.glm
dat.umap.lda4 <- dat.umap.lda
mm.celltype.lst4 <- mm.celltype.lst



# Write to output ---------------------------------------------------------

# for each mark separate
dat.lst <- list(dat.umap.glm.fillNAs1, dat.umap.glm.fillNAs2, dat.umap.glm.fillNAs3, dat.umap.glm.fillNAs4)
names(dat.lst) <- c(jmark1, jmark2, jmark3, jmark4)

# write each table with firs tcolumn cell, second column cluster

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

lapply(jmarks, function(jmark){
  print(jmark)
  outbase <- paste0("BM_AllMerged.", jmark, ".cell_cluster_table")
  outtxt <- file.path(outdir, paste0(outbase, ".txt"))
  outpdf <- file.path(outdir, paste0(outbase, ".pdf"))
  pdf(outpdf, useDingbats = FALSE)
    dat <- subset(dat.lst[[jmark]], select = c(cell, cluster, umap1, umap2, topic, topic.weight, cluster.orig, plate))
    data.table::fwrite(dat, file = outtxt, sep = "\t")
    m <- ggplot(dat, aes(x = umap1, y = umap2, color = cluster)) + geom_point() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      scale_color_manual(values = cbPalette, na.value = "grey85")
    print(m)
  dev.off()
})
dev.off()




