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

# indir.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2"
# count.mat.lst <- lapply(jmarks, function(jmark){
#   inf.lda <- file.path(indir.lda, paste0("lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj"))
#   load(inf.lda, v=T)
#   return(count.mat)
# })

indir.mat <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS.AnnotatedGeneRegionsWithPromsEnhs"
count.mat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.mat <- file.path(indir.mat, paste0(jmark, ".countTableTSS.mapq_40.TSS_10000.blfiltered.csv"))
  ReadMatTSSFormat(inf.mat, as.sparse = TRUE, add.coord = TRUE)
})

# check rownames
rnames.all <- lapply(count.mat.lst, function(x) rownames(x))
lapply(rnames.all, length)

rnames.common <- Reduce(intersect, rnames.all)

print(length(rnames.common))

count.mat.lst.filt <- lapply(count.mat.lst, function(x){
  x[rnames.common, ]
})
  
lapply(count.mat.lst.filt, dim)


# Celltype sepcific genes -------------------------------------------------

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


# Define each dot as either TSS or enhancer -------------------------------

coords.vec <- sapply(rnames.common, function(x) strsplit(x, ";")[[1]][[1]])

# # check dupes
# jdupes <- which(duplicated(coords.vec))
# coords.vec[which(coords.vec %in% jdupes)]

# what to do with TSS of genes with exact same TSS? remove it? yeah... they should get same label
# e.g. chrX:135732733-135752733;Armcx5+ chrX:135732733-135752733;Gprasp1
coords.vec <- coords.vec[!duplicated(coords.vec)]

jinf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/first_transcript_tss/gene_tss_winsize.10000.first_transcript.bed"
# annotate coord to gene
jchromos.keep <- c(paste("chr", seq(19), sep = ""), "chrX", "chrY")
bins.annot <- AnnotateCoordsFromList(coords.vec = coords.vec, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos.keep)


jctype <- "H3K27me3"
count.sub <- data.frame(region_coord_full = rownames(count.mat.pbulk.all.norm.wide), 
                        region_coord = sapply(rownames(count.mat.pbulk.all.norm.wide), function(x) strsplit(x, split = ";")[[1]][[1]]),
                        count.mat.pbulk.all.norm.wide[, grepl(jctype, colnames(count.mat.pbulk.all.norm.wide))], stringsAsFactors = FALSE) %>%
  left_join(., bins.annot$out2.df.closest) %>%
  ungroup() %>%
  mutate(abs.dist.to.tss = abs(dist.to.tss)) %>%
  mutate(is.enhancer = grepl("enhancer", region_coord_full))

# plot mark 
counts.sub.long <- data.table::melt(data.frame(region_coord_full = count.sub$region_coord_full, count.sub[, grepl(jctype, colnames(count.sub))]), 
                                    id.vars = "region_coord_full", variable.name = "pseudobulk", value.name = "exprs") %>%
  mutate(is.enhancer = grepl("enhancer", region_coord_full))

ggplot(counts.sub.long, aes(x = exprs, fill = is.enhancer)) + geom_density(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~pseudobulk, ncol = 1) 
  # facet_grid(is.enhancer ~ pseudobulk)


ggplot(count.sub, aes_string(x = paste0("H3K4me3_", jctype), y = paste0("H3K4me1_", jctype), color = "is.enhancer")) + geom_point(alpha = 0.25) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~is.enhancer) + ggtitle(paste0("H3K4me1 vs H3K4me3 on ENCODE annotations. Enhancer (TRUE) vs promoter (FALSE)"))


ggplot(count.sub, aes_string(x = paste0("H3K4me3_", jctype), y = paste0("H3K4me1_", jctype), color = "is.enhancer")) + geom_point(alpha = 0.25) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~is.enhancer)

ggplot(count.sub, aes_string(x = paste0("H3K4me3_", jctype), y = paste0("H3K27me3_", jctype), color = "is.enhancer")) + geom_point(alpha = 0.25) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~is.enhancer)

ggplot(count.sub %>% filter(!is.na(abs.dist.to.tss)) %>% arrange(desc(abs.dist.to.tss)), 
       aes_string(x = paste0("H3K4me3_", jctype), y = paste0("H3K4me1_", jctype), color = "abs.dist.to.tss")) + 
  geom_point(alpha = 0.1) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) 


  
ggplot(count.sub, aes_string(x = "H3K4me3_", jctype)) + geom_density() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(count.sub, aes_string(x = H3K27me3_HSCs)) + geom_density()

# Label annotations  ------------------------------------------------------



