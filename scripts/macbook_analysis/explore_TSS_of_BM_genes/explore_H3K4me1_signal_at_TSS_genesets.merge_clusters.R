# Jake Yeung
# Date of Creation: 2020-03-12
# File: ~/projects/scchic/scripts/macbook_analysis/explore_TSS_of_BM_genes/explore_H3K4me3_signal_at_TSS.R
# 

rm(list=ls())

# load libs ---------------------------------------------------------------


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(reticulate)

library(DescTools)

library(purrr)
library(ggrastr)
library(pracma)

library(hash)
library(JFuncs)

library(scchicFuncs)

library(forcats)
library(ggrepel)

library(RColorBrewer)

library(viridis)

markers.lst <- MarkerToCelltype()

# markers.lst <- markers.lst

reticulate::source_python("/Users/yeung/projects/scchic/scripts/python_functions/parse_dictionary_text.py")


# Preamble ----------------------------------------------------------------


jclusts <- c("Ltf", "Fcrla", "Ccl5", "Prss34", "Cd74", "Siglech", "Car1", "core")
names(jclusts) <- jclusts
cnames.rearranged <- c("Neutrophils", "Bcells", "InnateLymph", "Basophils", "Dendritic", "pDendritic", "Eryth.Sox6", "Eryth.Cdk6", "Eryth.Gfi1b","HSCs.Lrp5", "HSCs.Hlf", "HSCs.Msi2")
cnames.rearranged.full <- c("Neutrophils", "Bcells", "InnateLymph", "Basophils", "Dendritic", "pDendritic", "Eryth", "HSCs.Hlf", "HSCs.other")
cbPalette.all <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
cbPalette <- cbPalette.all[1:length(jclusts)]
names(cbPalette) <- jclusts


pseudos.merge.lst <- list("InnateLymph" = c("H3K4me1-BM_AllMerged.ILC-PrkchPlus_topic18.sorted.100", "H3K4me1-BM_AllMerged.ILC-RoraPlus_topic11.sorted.100"),
                          "Bcells" = c("H3K4me1-BM_AllMerged.Bcells-Cd47_topic29.sorted.100", "H3K4me1-BM_AllMerged.Bcells-Cd83_topic10.sorted.100"),
                          "HSCs.other"  = c("H3K4me1-BM_AllMerged.HSCs-Anxa2_topic22.sorted.100", "H3K4me1-BM_AllMerged.HSCs-Ephb2_topic5.sorted.100"))

# Load annots -------------------------------------------------------------

inf.annot <- "/Users/yeung/data/scchic/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat.celltypes_filt.rds"

dat.annots <- readRDS(inf.annot)
jmark <- "H3K4me1"
jdir <- "greaterthan"
if (jdir == "greaterthan"){
  logfcmin <- 0.5
} else if (jdir == "lessthan"){
  logfcmin <- -0.5
} else {
  warning("jdir must be greaterthan or lessthan:", jdir)
}
pvalmax <- 0.01

# Set output --------------------------------------------------------------

outpdf <- paste0("/Users/yeung/data/scchic/pdfs/marker_genes_Giladi_TSS_signal/heatmap_and_umap_", jmark, "_GiladiMarkerGenes.", jdir, ".pvallmax_", pvalmax, ".logfcmin_", logfcmin, ".mergeclusts.pdf")

# LLoad UMAPs -------------------------------------------------------------

inf.annot <- paste0("/Users/yeung/data/scchic/from_rstudio/pdfs_all/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
load(inf.annot, v=T)

# Load bed ----------------------------------------------------------------

inf.bed <- "/Users/yeung/data/databases/gene_tss/gene_tss/giladi_filtered.seurat/gene_tss_winsize.2.diff_exprs_Giladi_seurat.celltypes_filt.bed"
dat.bed <- fread(inf.bed)

coord.hash <- hash::hash(sapply(dat.bed$V4, function(x) strsplit(x, ";")[[1]][[1]]), sapply(dat.bed$V4, function(x) strsplit(x, ";")[[1]][[2]]))

# Load data  --------------------------------------------------------------

inf <- paste0("/Users/yeung/data/scchic/from_cluster/bigwig_outputs/merged_bams.deeptools_outputs.tss.MAPQ_40.dist_10000.allctypes_from_seurat.", jmark, ".bsize_100/computeMatrix.MAPQ_40.", jmark, ".gene_tss_winsize.2.diff_exprs_Giladi_seurat.celltypes_filt.tab.gz")
assertthat::assert_that(file.exists(inf))

meta <- inf2dic(inf)
mat.all <- fread(inf, header = FALSE, skip = 1, sep = "\t", quote = "")

# name columns 
nbins <- (unique(meta$downstream) + unique(meta$upstream)) / unique(meta[["bin size"]])
nsamps <- length(meta$sample_labels)

# Stack matrix vertically  ------------------------------------------------

mat <- mat.all[, -seq(6)]
coords <- mat.all[, c(seq(6))]
colnames(coords) <- c("chromo", "start", "end", "coord", "meta1", "meta2")
assertthat::assert_that(ncol(mat) == nsamps * nbins)

# split into a list of matrices, then rbind
sampids <- ceiling(seq(ncol(mat)) / nbins)
sampids.uniq <- as.character(unique(sort(sampids)))
names(sampids.uniq) <- meta$sample_labels

colids <- ( (seq(ncol(mat)) - 1) %% nbins ) + 1  # - 1 and +1 so first element is 1, last element is 10

colnames(mat) <- as.character(sampids)

mats.lst <- lapply(sampids.uniq, function(sampid){
  cols.i <- which(colnames(mat) == sampid)
  mat.sub <- mat[, ..cols.i]
  colnames(mat.sub) <- paste("bin", seq(nbins), sep = "")
  return(mat.sub) 
})

# make long dataframe of bed locations
coords.lst <- lapply(sampids.uniq, function(sampid){
  jtmp <- coords[, c(1,2,3,4)]
  jtmp$sampid <- sampid
  return(jtmp)
})  

mats.long <- do.call(rbind, mats.lst)
coords.long <- do.call(rbind, coords.lst)

mat.long.merge <- cbind(mats.long, coords.long)


# Make gene exprs ---------------------------------------------------------

mats.lst.clean <- lapply(mats.lst, function(jmat){
  jmat <- as.matrix(jmat)
  # jmat.win <- Winsorize(jmat, minval = 0, maxval = quantile(jmat, 0.98))
  jmat.win <- jmat
  return(jmat.win)
}) 

samp.remove <- names(mats.lst.clean)[grepl("BM_AllMerged..sorted.100", names(mats.lst.clean))]

for (jsamp in samp.remove){
  mats.lst.clean[[jsamp]] <- NULL
}



for (pseudos.merge.newname in names(pseudos.merge.lst)){
  # merge Eryths
  # https://stackoverflow.com/questions/26018216/calculating-mean-of-multiple-matrices-in-r
  pseudos.merge <- pseudos.merge.lst[[pseudos.merge.newname]]
  mats.lst.clean[[pseudos.merge.newname]] <- Reduce(mats.lst.clean[pseudos.merge], f = '+') / length(pseudos.merge)
  # remove old ones
  for (jname in pseudos.merge){
    mats.lst.clean[[jname]] <- NULL
  }
}

# pseudos.merge <- c("H3K4me3-BM_AllMerged.Eryth-Gfi1b_topic7.sorted.100", "H3K4me3-BM_AllMerged.Eryth-Cdk6_topic9.sorted.100", "H3K4me3-BM_AllMerged.Eryth-Sox6_topic16.sorted.100")
# pseudos.merge.newname <- "Eryth"
# # mats.lst.clean[[pseudos.merge.newname]] <- purrr::reduce(mats.lst.clean[pseudos.merge], .f = "+") / length(pseudos.merge)
# mats.lst.clean[[pseudos.merge.newname]] <- Reduce(mats.lst.clean[pseudos.merge], f = '+') / length(pseudos.merge)
# # remove old ones
# for (jname in pseudos.merge){
#   mats.lst.clean[[jname]] <- NULL
# }
# 
# # merge HSs
# pseudos.merge <- c("H3K4me3-BM_AllMerged.HSCs-Hlf_topic26.sorted.100", "H3K4me3-BM_AllMerged.HSCs-Lrp5_topic14.sorted.100")
# pseudos.merge.newname <- "HSCs"
# # mats.lst.clean[[pseudos.merge.newname]] <- purrr::reduce(mats.lst.clean[pseudos.merge], .f = "+") / length(pseudos.merge)
# mats.lst.clean[[pseudos.merge.newname]] <- Reduce(mats.lst.clean[pseudos.merge], f = '+') / length(pseudos.merge)
# # remove old ones
# for (jname in pseudos.merge){
#   mats.lst.clean[[jname]] <- NULL
# }

gene.exprs <- lapply(mats.lst.clean, function(jmat){
  rowMeans(jmat)
})

cpm.mat <- do.call(cbind, gene.exprs)

# shorten colnames
colnames(cpm.mat) <- gsub(".sorted.100$", "", colnames(cpm.mat))
colnames(cpm.mat) <- gsub(".BM_AllMerged", "", colnames(cpm.mat))
colnames(cpm.mat) <- gsub(paste0("^", jmark, "."), "", colnames(cpm.mat))
colnames(cpm.mat) <- make.names(sapply(colnames(cpm.mat), function(x) strsplit(x, "_")[[1]][[1]]))





rownames(cpm.mat) <- coords$coord

cpm.mat.long <- data.frame(coord = rownames(cpm.mat), gene = sapply(rownames(cpm.mat), function(x) AssignHash(x, coord.hash, NA)), cpm.mat, stringsAsFactors = FALSE)
cpm.mat.long <- tidyr::gather(cpm.mat.long, key = "pseudobulk", value = "cpm", -c(coord, gene)) %>%
  group_by(coord) %>%
  # mutate(cpm = log2(cpm * 10^6 + 1)) %>%
  mutate(zscore = scale(cpm, center = TRUE, scale = TRUE))

# get diff exprs genes ----------------------------------------------------

# jclusts <- as.character(unique(dat.annots$cluster))
# names(jclusts) <- jclusts

ctype.genes <- lapply(jclusts, function(jclust){
  if (jdir == "greaterthan"){
    jsub <- subset(dat.annots, cluster == jclust & p_val_adj <= pvalmax & avg_logFC >= logfcmin) %>%
      arrange(desc(abs(avg_logFC)))
  }  else if (jdir == "lessthan"){
    jsub <- subset(dat.annots, cluster == jclust & p_val_adj <= pvalmax & avg_logFC <= -logfcmin) %>%
      arrange(desc(abs(avg_logFC)))
  }
  jgenes <- jsub$gene
  # filter out ones that are not in rowwnames of matrix
  jgenes <- jgenes[which(jgenes %in% unique(cpm.mat.long$gene))]
})


pdf(outpdf, useDingbats = FALSE)
for (jclust in jclusts){
  jgenes <- ctype.genes[[jclust]]
  
  # select gene with most variability across clusters
  cpm.varfilt.long <- cpm.mat.long %>%
    group_by(coord, gene) %>% 
    summarise(jsd = sd(cpm)) %>%
    group_by(gene) %>%
    filter(jsd == max(jsd))
  coords.keep <- cpm.varfilt.long$coord
  
  cpm.long.filt <- subset(cpm.mat.long, gene %in% jgenes & coord %in% coords.keep) %>%
    ungroup()
  cpm.mat.filt <- as.data.frame(tidyr::spread(cpm.long.filt %>% dplyr::select(-coord, -cpm), key = pseudobulk, value = zscore))
  rownames(cpm.mat.filt) <- cpm.mat.filt$gene
  cpm.mat.filt$gene <- NULL
  
  cpm.mat.filt.ordered <- cpm.mat.filt[jgenes, ]
  heatmap3::heatmap3(cpm.mat.filt.ordered, Rowv = NA, Colv = NA, 
                    main = paste("Clust;", jclust, "(", markers.lst[[jclust]], "). Ngenes:", length(jgenes)), revC = TRUE, labRow = FALSE, 
                    col = colorRampPalette(viridis(12))(1024))
}

# do all genes together?
jgenes <- unlist(ctype.genes)
cpm.long.filt <- subset(cpm.mat.long, gene %in% jgenes & coord %in% coords.keep) %>%
  ungroup()
cpm.mat.filt <- as.data.frame(tidyr::spread(cpm.long.filt %>% dplyr::select(-coord, -cpm), key = pseudobulk, value = zscore))
rownames(cpm.mat.filt) <- cpm.mat.filt$gene
cpm.mat.filt$gene <- NULL
cpm.mat.filt.ordered <- cpm.mat.filt[jgenes, ]


# order based on entropy? 

# manually arrange celltypes


gene.colors <- lapply(jclusts, function(jclust){
  return(rep(cbPalette[[jclust]], length(ctype.genes[[jclust]])))
}) %>%
  unlist()

print(sort(colnames(cpm.mat.filt.ordered)))
print(sort(cnames.rearranged.full))
heatmap3::heatmap3(cpm.mat.filt.ordered[, cnames.rearranged.full], 
                   Rowv = NA, Colv = NA, margins = c(5, 1),
                  main = paste("\n\n\n", jmark, 'signal at TSS of\ncelltype specific BM genes'), 
                  revC = TRUE, labRow = FALSE, 
                  col = colorRampPalette(viridis(12))(1024), 
                  RowSideColors = gene.colors, RowSideLabs = "")


# Plot the clusters in the UMAP again -------------------------------------

# plot UMAP with same colors as heatmap
dat.umap.glm.fillNAs.mod <- dat.umap.glm.fillNAs %>%
  rowwise() %>%
  mutate(cluster = make.names(strsplit(cluster, split = "_")[[1]][[1]])) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = cnames.rearranged))

m.umap <- ggplot(dat.umap.glm.fillNAs.mod, aes(x = umap1, y = umap2, color = cluster)) + geom_point(size = 3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette.all, na.value = "grey95") + ggtitle(jmark)
print(m.umap)

dev.off()



