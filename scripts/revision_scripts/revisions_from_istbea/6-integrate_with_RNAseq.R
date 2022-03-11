# Jake Yeung
# Date of Creation: 2022-01-29
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/6-integrate_with_RNAseq.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(Seurat)

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/analysis_scrnaseq_atacseq"
outpdf <- file.path(outdir, paste0("plots_CCA_celltype_colored.", Sys.Date(), ".pdf"))
outmeta <- file.path(outdir, paste0("metadata_CCA_celltype_colored.", Sys.Date(), ".txt"))


# Load data ---------------------------------------------------------------

inf.cca <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/analysis_scrnaseq_atacseq/scRNAseqBaccin_sortChICZeller_CCA_comparison.RData"
load(inf.cca, v=T)


# Show UMAP  --------------------------------------------------------------

pdf(outpdf, useDingbats = FALSE)

ggplot(dat.umap.cca.rna.chic.annot, aes(x = umap1, y = umap2, color = experi, shape = experi)) + 
  geom_point(alpha = 0.5) + 
  scale_shape_manual(values = c(3, 16)) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.cca.rna.chic.annot, aes(x = umap1, y = umap2, color = experi, shape = experi)) + 
  geom_point(alpha = 0.5) + 
  scale_shape_manual(values = c(3, 16)) + 
  theme_bw() + 
  facet_wrap(~cluster) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(sort(unique(dat.umap.cca.rna.chic.annot$cluster)))

# Show where are discrepancies --------------------------------------------

rename.lst <- list("B cell" = "Bcells",
                   "Eo/Baso prog." = "Baso/Eo", 
                   "Dendritic cells" = "DCs", 
                   "Granulocytes" = "Neutrophils", 
                   "NK cells" = "NKs", 
                   "Erythroblasts" = "Eryths",
                   "T cells" = "Tcells",
                   "Basophils" = "Baso/Eo")

rename.hash <- hash::hash(rename.lst)

dat.umap.cca.rna.chic.annot <- dat.umap.cca.rna.chic.annot %>%
  rowwise() %>%
  mutate(cluster.rename = AssignHash(x = cluster, jhash = rename.hash, null.fill = cluster))


jshuffle <- dat.umap.cca.rna.chic.annot[sample(nrow(dat.umap.cca.rna.chic.annot)), ]

jshuffle <- jshuffle %>%
  rowwise() %>%
  mutate(experi = ifelse(experi == "ChIC", "sortChIC", "scRNAseq"))

ggplot(jshuffle, aes(x = umap1, y = umap2, color = experi, shape = experi)) + 
  geom_point(alpha = 0.5) + 
  scale_shape_manual(values = c(3, 16)) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jshuffle, aes(x = umap1, y = umap2, color = experi, shape = experi)) + 
  geom_point(alpha = 0.5) + 
  scale_shape_manual(values = c(3, 16)) + 
  theme_bw() + 
  facet_wrap(~cluster.rename) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


# Color code --------------------------------------------------------------

inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/metadata/metadata_batch_corrected.arranged_by_lineage.H3K4me1.2021-01-02.txt"
dat.meta <- fread(inf.meta)

clsts.old <- unique(dat.meta$cluster)
colors.old <- unique(dat.meta$clustercol)


clst2color.old <- hash::hash(clsts.old, colors.old)
clst2color.old["Neutrophils"] <- "#D55E00"  # same as Granus
clst2color.old["Baso/Eo"] <- "#696969"
clst2color.old["Gran/Mono prog."] <- "#d77120"
clst2color.old["Ery prog."] <- "#1b7eb5"
clst2color.old["Ery/Mk prog."] <- "#3688b7"
clst2color.old["large pre-B."] <- "#85c3e7"
clst2color.old["pro-B"] <- "#b8d6e7"
clst2color.old["small pre-B."] <- "#6490aa"
clst2color.old["LMPPs"] <- "#ca4890"
clst2color.old["Mono prog."] <- "#bd8151"
clst2color.old["Monocytes"] <- "#ba631e"
clst2color.old["Neutro prog."] <- "#d59563"
clst2color.old["Tcells"] <- "#019f01"

jshuffle.color <- jshuffle %>%
  rowwise() %>%
  mutate(clustercol = AssignHash(cluster.rename, clst2color.old, null.fill = "#7CA4A5"))


ggplot(jshuffle.color, aes(x = umap1, y = umap2, color = clustercol, shape = experi)) + 
  geom_point() + 
  scale_shape_manual(values = c(3, 16)) + 
  scale_color_identity() + 
  facet_wrap(~cluster.rename) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jbreaks <- hash::values(clst2color.old)
jlabels <- base::names(clst2color.old)

ggplot(jshuffle.color, aes(x = umap1, y = umap2, color = clustercol, shape = experi)) + 
  geom_point() + 
  theme_bw() + 
  xlab("CCA-umap1") + 
  ylab("CCA-umap2") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_shape_manual(values = c(3, 16)) + 
  scale_color_identity(name = "Celltype",
                       breaks = jbreaks,
                       labels = jlabels,
                       guide = "legend")

ggplot(jshuffle.color, aes(x = umap1, y = umap2, color = clustercol, shape = experi)) + 
  geom_point() + 
  theme_bw() + 
  xlab("CCA-umap1") + 
  ylab("CCA-umap2") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_shape_manual(values = c(3, 16)) + 
  scale_color_identity(name = "Celltype",
                       breaks = jbreaks,
                       labels = jlabels,
                       guide = "legend")


ggplot(jshuffle.color, aes(x = umap1, y = umap2, color = clustercol, shape = experi)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~experi) + 
  xlab("CCA-umap1") + 
  ylab("CCA-umap2") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_shape_manual(values = c(3, 16)) + 
  scale_color_identity(name = "Celltype",
                       breaks = jbreaks,
                       labels = jlabels,
                       guide = "legend")

ggplot(jshuffle.color, aes(x = umap1, y = umap2, color = clustercol, shape = experi)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~experi) + 
  xlab("CCA-umap1") + 
  ylab("CCA-umap2") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_shape_manual(values = c(3, 16)) + 
  scale_color_identity(name = "Celltype",
                       breaks = jbreaks,
                       labels = jlabels,
                       guide = "legend")

ggplot(jshuffle.color, aes(x = umap1, y = umap2, color = experi)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~cluster.rename) + 
  xlab("CCA-umap1") + 
  ylab("CCA-umap2") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Show each celltype ------------------------------------------------------

jctypes <- sort(unique(jshuffle.color$cluster.rename))

for (jctype in jctypes){
  
  jtmp <- jshuffle.color %>%
    rowwise() %>%
    mutate(is.ctype = cluster.rename == jctype, 
           clustercol.ctype = ifelse(cluster.rename == jctype, clustercol, "grey80"))  %>%
    arrange(is.ctype) 
  m <- ggplot(jtmp, aes(x = umap1, y = umap2, color = clustercol.ctype, shape = experi)) + 
    geom_point(alpha = 1) + 
    ggtitle(paste(jctype)) + 
    theme_bw() + 
    xlab("CCA-umap1") + 
    ylab("CCA-umap2") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_shape_manual(values = c(3, 16)) + 
    scale_color_identity(name = "Celltype",
                         breaks = jbreaks,
                         labels = jlabels,
                         guide = "legend")
  print(m)
  print(m + facet_wrap(~experi))
  
}



dev.off()

# Save outputs ------------------------------------------------------------
 
fwrite(jshuffle.color, file = outmeta, sep = "\t")



