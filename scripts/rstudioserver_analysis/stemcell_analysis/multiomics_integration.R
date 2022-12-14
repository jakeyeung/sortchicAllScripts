# Jake Yeung
# Date of Creation: 2020-03-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/multiomics_integration.R
# Multi-omics integration of celltypes 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

# get raw counts ----------------------------------------------------------


indir.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2"

count.mat.lst <- lapply(jmarks, function(jmark){
  inf.lda <- file.path(indir.lda, paste0("lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj"))
  load(inf.lda, v=T)
  return(count.mat)
})

# check rownames
rnames.all <- lapply(count.mat.lst, function(x) rownames(x))

lapply(rnames.all, length)

rnames.common <- Reduce(intersect, rnames.all)

count.mat.lst.filt <- lapply(count.mat.lst, function(x){
  x[rnames.common, ]
})
  
lapply(count.mat.lst.filt, dim)


# Load annots  ------------------------------------------------------------

indir.annots <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping")

dat.annots.all <- lapply(jmarks, function(jmark){
  inf.annots <- file.path(indir.annots, paste0("GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData"))
  load(inf.annots, v=T)
  return(dat.umap.glm.fillNAs)
})

# Make pseudobulks  -------------------------------------------------------

ctypes.keep <- c("HSC", "Neutrophil", "Bcell", "Innate", "Eryth")

cnames.keep.lst.all <- lapply(jmarks, function(jmark){
  jsplit <- split(dat.annots.all[[jmark]], dat.annots.all[[jmark]]$cluster)
  cnames.keep <- lapply(jsplit, function(x) x$cell)
})

count.pseudos <- lapply(jmarks, function(jmark){
  exprs.lst <- SumAcrossClusters(count.mat.lst.filt[[jmark]], cnames.keep.lst.all[[jmark]])
  exprs.mat <- as.data.frame(do.call(cbind, exprs.lst))
})

lapply(count.pseudos, function(x) colnames(x))

# filter out relevant celltypes and create "3D" matrix?
h3k4me1.cnames <- list("Neutrophils_topic23" = "Neutro", 
                       "HSCs-Hlf_topic7" = "HSCs",
                       "Eryth_topic27" = "Eryth",
                       "ILC-RoraPlus_topic11" = "NKcells",
                       "Bcells-Cd47_topic29" = "Bcells")
                       # "Bcells-Cd83_topic10" = "Bcells")

h3k4me3.cnames <- list("Bcells_topic13" = "Bcells", 
                       "Eryth-Sox6_topic16" = "Eryth",
                       "HSCs-Hlf_topic26" = "HSCs",
                       "Neutrophils_topic2" = "Neutro",
                       "InnateLymph_topic27" = "NKcells")

h3k27me3.cnames <- list("Bcells_topic16" = "Bcells", 
                       "Eryth-Sox6-_topic6" = "Eryth",
                       "HSCs-Tead1-_topic9" = "HSCs",
                       "Neutrophils_topic22" = "Neutro",
                       "InnateLymph_topic27" = "NKcells")

print(jmarks)
jmarks.cnames <- list(H3K4me1 = h3k4me1.cnames,
                      H3K4me3 = h3k4me3.cnames,
                      H3K27me3 = h3k27me3.cnames)

# Rename column names, collapse names with same names ---------------------

count.mat.pbulk.all <- lapply(jmarks, function(jmark){
  cnames <- jmarks.cnames[[jmark]]
  count.mark <- lapply(names(cnames), function(cname){
    x <- count.pseudos[[jmark]][[cname]]
    assertthat::assert_that(!is.null(x))
    names(x) <- rnames.common
    return(x)
  })
  names(count.mark) <- names(cnames)
  # rename
  names(count.mark) <- sapply(names(cnames), function(x) cnames[[x]])
  count.mat <- do.call(cbind, count.mark)
})


# What to do now?  --------------------------------------------------------

# normalize 
library(DESeq2)

coldat.all <- lapply(count.mat.pbulk.all, function(jmat){
  cdat <- data.frame(pseudobulk = colnames(jmat), stringsAsFactors = FALSE)
  rownames(cdat) <- cdat$pseudobulk
  return(cdat)
})

count.mat.pbulk.all.norm <- lapply(jmarks, function(jmark){
  ds <- DESeqDataSetFromMatrix(count.mat.pbulk.all[[jmark]], colData = coldat.all[[jmark]], design = ~1)
  ds.vst <- assay(vst(ds))
  return(ds.vst)
})

plot(density(unlist(count.mat.pbulk.all.norm$H3K4me1)))
plot(density(unlist(count.mat.pbulk.all.norm$H3K4me3)))
plot(density(unlist(count.mat.pbulk.all.norm$H3K27me3)))

# just do standard PCA 
count.mat.pbulk.all.norm.wide <- lapply(jmarks, function(jmark){
  colnames(count.mat.pbulk.all.norm[[jmark]]) <- paste(jmark, colnames(count.mat.pbulk.all.norm[[jmark]]), sep = "_")
  return(count.mat.pbulk.all.norm[[jmark]])
}) 
count.mat.pbulk.all.norm.wide <- do.call(cbind, count.mat.pbulk.all.norm.wide)

pca.out <- prcomp(t(count.mat.pbulk.all.norm.wide), center = TRUE, scale. = TRUE)

dat.pca <- data.frame(sample = rownames(pca.out$x), pca.out$x, stringsAsFactors = FALSE)

ggplot(dat.pca, aes(x = PC1, y = PC2, label = sample)) + geom_point() + geom_text_repel() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pca, aes(x = PC2, y = PC3, label = sample)) + geom_point() + geom_text_repel() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# plot loadings
dat.loadings <- data.frame(bin = rownames(pca.out$rotation), pca.out$rotation, stringsAsFactors = FALSE) %>%
  arrange(desc(PC3))
  # arrange(pc2)

print(head(dat.loadings))

# show top pc1

# plot a bin across conditions 
jbin <- dat.loadings$bin[[1]]

qplot(x = names(count.mat.pbulk.all.norm.wide[jbin, ]), y = count.mat.pbulk.all.norm.wide[jbin, ], 
      label = names(count.mat.pbulk.all.norm.wide[jbin, ]), geom = "point") + geom_text_repel() + 
  theme_bw() + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

boxplot(count.mat.pbulk.all.norm.wide)
  
# Plot genome-wide summary of neutophils  ---------------------------------

count.sub <- count.mat.pbulk.all.norm.wide[, grepl("Eryth", colnames(count.mat.pbulk.all.norm.wide))]
plot(count.sub[, 1], count.sub[, 2], pch = 20)
plot(count.sub[, 1], count.sub[, 3], pch = 20)
plot(count.sub[, 2], count.sub[, 3], pch = 20)

# filter out TSS? 


# Label each gene by TSS or not -------------------------------------------


jinf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/first_transcript_tss/gene_tss_winsize.10000.first_transcript.bed"
# annotate coord to gene 
jchromos.keep <- c(paste("chr", seq(19), sep = ""), "chrX", "chrY")
bins.annot <- AnnotateCoordsFromList(coords.vec = rownames(count.mat.pbulk.all.norm.wide), inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos.keep)


# Get neutro specific genes -----------------------------------------------


pval.min <- 0.01
fc.min <- 0
inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat.celltypes_filt.rds"
dat.de <- readRDS(inf.de)
neutro.genes <- subset(dat.de, p_val_adj < pval.min & cluster == "Ltf" & avg_logFC > fc.min)$gene

# label genes as cluster specific

jclsts <- as.character(unique(dat.de$cluster))
names(jclsts) <- jclsts

de.genes.lst <- lapply(jclsts, function(jclst){
  subset(dat.de, p_val_adj < pval.min & cluster == jclst & avg_logFC)$gene
})

# Plot H3K4me1 vs H3K4me3, plot distance to nearest TSS  ------------------

jctype <- "HSCs"
jctype <- "Eryth"
jctype <- "Bcells"
jctype <- "H3K4me3"

count.sub <- data.frame(region_coord = rownames(count.mat.pbulk.all.norm.wide), 
                        count.mat.pbulk.all.norm.wide[, grepl(jctype, colnames(count.mat.pbulk.all.norm.wide))], stringsAsFactors = FALSE) %>%
  left_join(., bins.annot$out2.df.closest) %>%
  ungroup() %>%
  mutate(abs.dist.to.tss = abs(dist.to.tss)) %>%
  mutate() %>%
  # mutate(DistCategory = ifelse(abs.dist.to.tss < median(abs.dist.to.tss, na.rm = TRUE), "LessThanMedian", "GreaterThanMedian"))
  mutate(DistCategory = ifelse(abs.dist.to.tss < 1000, "LessThanMedian", "GreaterThanMedian"))

count.sub <- count.sub %>%
  rowwise() %>%
  mutate(is.neutro.specific = gene %in% neutro.genes)

jmed <- median(count.sub$abs.dist.to.tss, na.rm = TRUE)


ggplot(count.sub %>% filter(!is.na(dist.to.tss)) %>% arrange(is.neutro.specific), aes_string(y = paste0(jctype, "_Neutro"), x = paste0(jctype, "_HSCs"), color = "is.neutro.specific")) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  # scale_color_viridis_c(direction = -1) + 
  ggtitle(paste0("", jctype, ", annotate 50 kb bins to dist to TSS"))

ggplot(count.sub %>% filter(!is.na(dist.to.tss)) %>% arrange(is.neutro.specific), aes_string(y = paste0(jctype, "_Neutro"), 
                                                                                             x = paste0(jctype, "_Bcells"), color = "is.neutro.specific")) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  # scale_color_viridis_c(direction = -1) + 
  ggtitle(paste0("", jctype, ", annotate 50 kb bins to dist to TSS"))


ggplot(count.sub %>% filter(!is.na(dist.to.tss)), aes_string(y = paste0("H3K4me1_", jctype), x = paste0("H3K4me3_", jctype), color = "abs.dist.to.tss")) + 
  geom_point(alpha = 0.25) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + ggtitle(paste0("H3K4me3 vs H3K4me1 ", jctype, ", annotate 50 kb bins to dist to TSS"))

ggplot(count.sub %>% filter(!is.na(dist.to.tss)), aes_string(y = paste0("H3K4me1_", jctype), x = paste0("H3K4me3_", jctype), color = "DistCategory")) + 
  geom_point(alpha = 0.25) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  # scale_color_viridis_c(direction = -1) + 
  ggtitle(paste0("H3K4me3 vs H3K4me1 ", jctype, ", annotate 50 kb bins to dist to TSS"), paste0("Median Distance=", jmed)) + facet_wrap(~DistCategory)

ggplot(count.sub %>% filter(!is.na(dist.to.tss)), aes_string(y = paste0("H3K4me1_", jctype), x = paste0("H3K4me3_", jctype), color = "DistCategory")) + 
  geom_point(alpha = 0.25) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste0("H3K4me3 vs H3K4me1 ", jctype, ", annotate 50 kb bins to dist to TSS"), paste0("Median Distance=", jmed))



