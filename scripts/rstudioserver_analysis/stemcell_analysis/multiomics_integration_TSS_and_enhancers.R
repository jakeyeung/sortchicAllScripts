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

library(DESeq2)

remove.eryths <- TRUE
jdist <- 1000L
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

m2c <- MarkerToCelltype()

# get raw counts ----------------------------------------------------------

# indir.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2"
# count.mat.lst <- lapply(jmarks, function(jmark){
#   inf.lda <- file.path(indir.lda, paste0("lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj"))
#   load(inf.lda, v=T)
#   return(count.mat)
# })

indir.mat <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS.AnnotatedGeneRegionsWithPromsEnhs")
count.mat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.mat <- file.path(indir.mat, paste0(jmark, ".countTableTSS.mapq_40.TSS_", jdist, ".blfiltered.csv"))
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
  subset(dat.de, p_val_adj < pval.min & cluster == jclst & avg_logFC > fc.min)$gene
})

subset(dat.de, gene == "S100a8")

# Load annots  ------------------------------------------------------------

indir.annots <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping")
dat.annots.all <- lapply(jmarks, function(jmark){
  inf.annots <- file.path(indir.annots, paste0("GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData"))
  load(inf.annots, v=T)
  return(dat.umap.glm.fillNAs)
})

indir.cellsizes <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.good_runs")
dat.cellsizes <- lapply(jmarks, function(jmark){
  inf.glmpca <- file.path(indir.cellsizes, paste0("PZ_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_1000.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.winsorize_TRUE.2020-02-11.RData"))
  load(inf.glmpca, v=T)
  return(glm.inits$size.factor)
})

# Make pseudobulks  -------------------------------------------------------

if (!remove.eryths){
  ctypes.keep <- c("HSC", "Neutrophil", "Bcell", "Innate", "Eryth")
} else {
  ctypes.keep <- c("HSC", "Neutrophil", "Bcell", "Innate")
}

cnames.keep.lst.all <- lapply(jmarks, function(jmark){
  jsplit <- split(dat.annots.all[[jmark]], dat.annots.all[[jmark]]$cluster)
  cnames.keep <- lapply(jsplit, function(x) x$cell)
})

count.pseudos <- lapply(jmarks, function(jmark){
  exprs.lst <- SumAcrossClusters(count.mat.lst.filt[[jmark]], cnames.keep.lst.all[[jmark]])
  exprs.mat <- as.data.frame(do.call(cbind, exprs.lst))
})

lapply(count.pseudos, function(x) colnames(x))

if (!remove.eryths){
  # filter out relevant celltypes and create "2D" matrix?
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
} else {
  # filter out relevant celltypes and create "2D" matrix?
  h3k4me1.cnames <- list("Neutrophils_topic23" = "Neutro", 
                         "HSCs-Hlf_topic7" = "HSCs",
                         "ILC-RoraPlus_topic11" = "NKcells",
                         "Bcells-Cd47_topic29" = "Bcells")
  # "Bcells-Cd83_topic10" = "Bcells")
  
  h3k4me3.cnames <- list("Bcells_topic13" = "Bcells", 
                         "HSCs-Hlf_topic26" = "HSCs",
                         "Neutrophils_topic2" = "Neutro",
                         "InnateLymph_topic27" = "NKcells")
  
  h3k27me3.cnames <- list("Bcells_topic16" = "Bcells", 
                          "HSCs-Tead1-_topic9" = "HSCs",
                          "Neutrophils_topic22" = "Neutro",
                          "InnateLymph_topic27" = "NKcells")
}

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

jinf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/first_transcript_tss/gene_tss_winsize.50000.first_transcript.bed"
# annotate coord to gene
jchromos.keep <- c(paste("chr", seq(19), sep = ""), "chrX", "chrY")
bins.annot <- AnnotateCoordsFromList(coords.vec = coords.vec, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos.keep)


jctype <- "HSCs"
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


# ggplot(count.sub, aes_string(x = paste0("H3K4me3_", jctype), y = paste0("H3K4me1_", jctype), color = "is.enhancer")) + geom_point(alpha = 0.25) +
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~is.enhancer) + ggtitle(paste0("H3K4me1 vs H3K4me3 on ENCODE annotations. Enhancer (TRUE) vs promoter (FALSE)"))
# ggplot(count.sub, aes_string(x = paste0("H3K4me3_", jctype), y = paste0("H3K4me1_", jctype), color = "is.enhancer")) + geom_point(alpha = 0.25) +
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~is.enhancer)
# ggplot(count.sub, aes_string(x = paste0("H3K4me3_", jctype), y = paste0("H3K27me3_", jctype), color = "is.enhancer")) + geom_point(alpha = 0.25) +
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~is.enhancer)
# ggplot(count.sub %>% filter(!is.na(abs.dist.to.tss)) %>% arrange(desc(abs.dist.to.tss)), 
#        aes_string(x = paste0("H3K4me3_", jctype), y = paste0("H3K4me1_", jctype), color = "abs.dist.to.tss")) + 
#   geom_point(alpha = 0.1) +
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_viridis_c(direction = -1) 
# ggplot(count.sub, aes_string(x = "H3K4me3_", jctype)) + geom_density() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# # ggplot(count.sub, aes_string(x = H3K27me3_HSCs)) + geom_density()

# Label annotations  ------------------------------------------------------

# check if enhancers are labeled with genes: yes if within a certain distance
plot(density((subset(count.sub, grepl("enhancer", region_coord_full) & !is.na(dist.to.tss))$dist.to.tss)))


# take only proms/enhs labeled with a gene, so we can ask if it is neutro, bcell, or HSC specific (housekeeping?)
count.sub.sub <- subset(count.sub, !is.na(gene))

# label genes as neutro, bcell, hsc specific, or others
head(count.sub.sub)

for (jclst.tmp in jclsts){
  jclst.tmp.name <- paste0("DE_", m2c[[jclst.tmp]])
  assertthat::assert_that(!is.null(jclst.tmp.name))
  count.sub.sub[[jclst.tmp.name]] <- sapply(count.sub.sub$gene, function(x) x %in% de.genes.lst[[jclst.tmp]])
}

subset(count.sub.sub, gene == "S100a8")

# show the expressed genes onto the HSCs graph
ggplot(count.sub.sub, aes(x = H3K4me3_HSCs, y = H3K27me3_HSCs, color = DE_HSCs)) + 
  geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(is.enhancer~DE_HSCs)



# Do everything in long form because it's probably better -----------------

counts.pbulk.long <- data.table::melt(data.frame(region_coord_full = rownames(count.mat.pbulk.all.norm.wide), 
                                               count.mat.pbulk.all.norm.wide, stringsAsFactors = FALSE), 
                                    id.vars = "region_coord_full", variable.name = "pseudobulk", value.name = "exprs") %>%
  ungroup() %>% 
  mutate(biotype = ifelse(grepl("enhancer", region_coord_full), "enhancer", "promoter"),
         region_coord = sapply(region_coord_full, function(x) strsplit(x, split = ";")[[1]][[1]]),
         mark = sapply(as.character(pseudobulk), function(x) strsplit(x, "_")[[1]][[1]]),
         ctype = sapply(as.character(pseudobulk), function(x) strsplit(x, "_")[[1]][[2]])) %>%
  left_join(., bins.annot$out2.df.closest) %>%
  mutate(abs.dist.to.tss = abs(dist.to.tss)) %>%
  filter(!is.na(gene))

# assign gene to celltype specificity

# needs to be a concatenated vector for now? 

gene2ctype <- hash::hash()
for (jclst in jclsts){
  for (jgene in de.genes.lst[[jclst]]){
    gene2ctype[[jgene]] <- c(gene2ctype[[jgene]], jclst)
  }
}
counts.pbulk.long$de.ctype <- sapply(counts.pbulk.long$gene, AssignHash, gene2ctype, null.fill = NA)

counts.pbulk.long <- counts.pbulk.long %>%
  group_by(region_coord_full) %>%
  mutate(de.ctype.choose = sample(de.ctype[[1]], size = 1))

counts.pbulk.long$de.ctype.choose <- sapply(counts.pbulk.long$de.ctype, function(x) sample(x = x, size = 1))

subset(counts.pbulk.long, !is.na(ctype))
dim(subset(counts.pbulk.long, !is.na(ctype)))
subset(counts.pbulk.long, gene == "S100a8")

# eg chr1:34449762-34469762;Ptpn18
subset(counts.pbulk.long, region_coord_full == "chr1:34449762-34469762;Ptpn18+")$de.ctype.choose


# Check celltype specificity  ---------------------------------------------

# on H3K4me3
# jmark <- "H3K27me3"
jmark <- "H3K4me3"
jmark <- "H3K4me1"
jsub <- subset(counts.pbulk.long, mark == jmark)
cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(jsub, aes(x = exprs, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
  facet_grid(ctype ~ de.ctype.choose) + 
  # facet_wrap(~biotype, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist))


# can I plot promoter and enhancer signal into an X-Y plot?
head(jsub)

jsub.sub <- subset(counts.pbulk.long, !is.na(gene), select = c(region_coord_full, exprs, biotype, ctype, gene, de.ctype.choose, mark)) %>%
  group_by(gene, biotype, ctype, de.ctype.choose, mark) %>%
  summarise(exprs = mean(exprs))  %>%
  reshape2::dcast(gene + mark + ctype + de.ctype.choose ~ biotype, value.var = "exprs")

ggplot(jsub.sub %>% filter(mark == "H3K4me3") %>% mutate(de.ctype.choose = ifelse(is.na(de.ctype.choose), "zNA", de.ctype.choose)) %>% arrange(desc(de.ctype.choose)), 
       aes(x = promoter, y = enhancer, color = de.ctype.choose)) + geom_point() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  facet_wrap(mark~ctype, scales = "free") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark)

ggplot(jsub, aes(x = exprs, group = biotype, fill = biotype)) + geom_density(alpha = 0.3) + 
  facet_grid(ctype~de.ctype.choose) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = cbPalette) + ggtitle(paste("Count distribution for", jmark, "genes split by DE category"))

ggplot(jsub, aes(x = exprs, group = biotype, fill = biotype)) + geom_density(alpha = 0.3) + 
  facet_grid(de.ctype.choose ~ ctype) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = cbPalette) + ggtitle(paste("Count distribution for", jmark, "genes split by DE category"))


# Plot H3K4me1 vs H3K27me3  -----------------------------------------------

# the random de.ctype.choose is problematic here!
jctype <- "Neutro"
jsub <- subset(counts.pbulk.long, ctype == jctype)
jsub.wide <- reshape2::dcast(subset(jsub, select = c(region_coord_full, exprs, mark, de.ctype.choose, biotype, gene, de.ctype.choose)), 
                             formula = "region_coord_full + gene + biotype + de.ctype.choose ~ mark", value.var = "exprs")
ggplot(jsub.wide %>% arrange(de.ctype.choose) %>% filter(!is.na(de.ctype.choose)), 
       aes(x = H3K4me3, y = H3K27me3, color = de.ctype.choose)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + ggtitle(jctype) + 
  facet_wrap(~de.ctype.choose)

# Compare one ctype with another ------------------------------------------

jmark <- "H3K4me3"
jmark <- "H3K4me1"


jmark <- "H3K27me3"
jmark <- "H3K4me3"

jmark <- "H3K4me1"

jsub <- subset(counts.pbulk.long, mark == jmark) %>%
  # filter(!ctype %in% c("HSCs", "NKcells")) %>%
  group_by(region_coord_full) %>%
  # mutate(exprs = scale(exprs, center = TRUE, scale = TRUE))
  mutate(exprs = scale(exprs, center = FALSE, scale = FALSE))

jsub.wide <- reshape2::dcast(subset(jsub, select = c(region_coord_full, de.ctype.choose, biotype, gene, de.ctype.choose, ctype, exprs)), 
                             formula = "region_coord_full  + biotype + gene + de.ctype.choose ~ ctype", value.var = "exprs")

# ggplot(jsub.wide %>% arrange(de.ctype.choose) %>% filter(!is.na(de.ctype.choose)), 
#        aes(x = HSCs, y = Bcells, color = biotype)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_manual(values = cbPalette, na.value = "grey85") + ggtitle(jmark) + 
#   facet_wrap(~de.ctype.choose)
# 
# ggplot(jsub.wide %>% arrange(de.ctype.choose), 
#        aes(x = Neutro, y = Bcells, color = biotype)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_manual(values = cbPalette, na.value = "grey85") + ggtitle(jmark) + 
#   facet_wrap(~de.ctype.choose)

# are HSCs in general higher enhancers? 
ggplot(jsub, aes(x = ctype, y = exprs, fill = biotype)) + geom_boxplot() + 
  facet_wrap(~de.ctype.choose) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark)


# Compare promoter and enhancer signal in K4me1 vs K4me3  -----------------

jbtype <- "promoter"
jsub <- counts.pbulk.long %>%
  group_by(region_coord_full, mark) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE),
         logFC = scale(exprs, center = TRUE, scale = FALSE))

# show H3k4me1, H3K4me3, H3K27me3 levels for each gene? 
ggplot(jsub %>% filter(mark != "H3K27me3"), aes(x = ctype, y = exprs, fill = mark)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(biotype ~ de.ctype.choose) + ggtitle(paste0("Bin size:", jdist)) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(jsub, aes(x = ctype, y = zscore, fill = mark)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(biotype ~ de.ctype.choose) + ggtitle(paste0("Bin size:", jdist)) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(jsub %>% filter(!is.na(de.ctype.choose) & mark != "H3K27me3"), aes(x = ctype, y = logFC, fill = mark)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(biotype ~ de.ctype.choose) + ggtitle(paste0("Bin size:", jdist)) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
# ggplot(jsub, aes(x = exprs, group = ctype, fill = ctype)) + geom_density() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_grid(biotype ~ de.ctype.choose) + 
#   ggtitle(paste0("Bin size:", jdist)) + 
#   xlab("Zscore") + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
#   



# Exprs suggests more reads are in H3K4me3 than H3K4me1?? -----------------

# check H3K4me1, H3K4me3, H3K27me3 total sizes

dat.cellsizes.gw <- lapply(jmarks, function(jmark){
  csize <- dat.cellsizes[[jmark]]
  jtmp <- data.frame(cell = names(csize), total.cuts = csize, stringsAsFactors = FALSE) %>%
    left_join(., dat.annots.all[[jmark]])
  jtmp$mark <- jmark
  return(jtmp)
})

ggplot(dat.cellsizes.gw %>% bind_rows(), aes(x = total.cuts, group = interaction(cond, mark), fill = cond)) + geom_density() + 
  facet_grid(mark ~ cond) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + geom_vline(aes(xintercept = median(total.cuts)))

# Plot promoter and enhancer signal in single cells  ----------------------

prom.counts <- lapply(jmarks, function(jmark){
  rows.keep <- rnames.common[!grepl("enhancer", rnames.common)]
  print(head(rows.keep))
  jmat <- count.mat.lst.filt[[jmark]][rows.keep, ]
  dat.tmp <- data.frame(cell = colnames(jmat), prom.cuts = colSums(jmat), stringsAsFactors = FALSE) %>%
    mutate(mark = jmark)
  return(dat.tmp)
}) 

enh.counts <- lapply(jmarks, function(jmark){
  rows.keep <- rnames.common[grepl("enhancer", rnames.common)]
  print(head(rows.keep))
  jmat <- count.mat.lst.filt[[jmark]][rows.keep, ]
  dat.tmp <- data.frame(cell = colnames(jmat), enh.cuts = colSums(jmat), stringsAsFactors = FALSE) %>%
    mutate(mark = jmark)
  return(dat.tmp)
})

dat.cellsizes.gw <- lapply(jmarks, function(jmark){
  x <- dat.cellsizes.gw[[jmark]]
  x <- left_join(x, prom.counts[[jmark]])
  x <- left_join(x, enh.counts[[jmark]])
  return(x)
})

for (jmark in jmarks){
  m <- ggplot(dat.cellsizes.gw[[jmark]], aes(x = umap1, y = umap2, color = prom.cuts / total.cuts)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = 1) + 
    facet_wrap(~mark)
  print(m)
}

for (jmark in jmarks){
  m <- ggplot(dat.cellsizes.gw[[jmark]], aes(x = umap1, y = umap2, color = enh.cuts / total.cuts)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = 1) + 
    facet_wrap(~mark)
  print(m)
}

for (jmark in jmarks){
  m <- ggplot(dat.cellsizes.gw[[jmark]], aes(x = umap1, y = umap2, color = (prom.cuts + enh.cuts) / total.cuts)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = 1) + 
    facet_wrap(~mark)
  print(m)
}

for (jmark in jmarks){
  m <- ggplot(dat.cellsizes.gw[[jmark]], aes(x = umap1, y = umap2, color = (1 - (prom.cuts + enh.cuts) / total.cuts))) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = 1) + 
    facet_wrap(~mark)
  print(m)
}

for (jmark in jmarks){
  m <- ggplot(dat.cellsizes.gw[[jmark]], aes(x = umap1, y = umap2, color = log10(total.cuts))) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = 1) + 
    facet_wrap(~mark)
  print(m)
}

# check genomewide
ggplot(dat.cellsizes.gw %>% bind_rows(), aes(x = prom.cuts, fill = cond)) + 
  geom_density() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(cond ~ mark) + scale_x_log10()

ggplot(dat.cellsizes.gw %>% bind_rows(), aes(x = enh.cuts, fill = cond)) + 
  geom_density() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(cond ~ mark) + scale_x_log10()


# Save objects for downstream exploration ---------------------------------

outrdata <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rsessions/multiomics_integration_proms_enhs_few_celltypes.RData"
save.image(file = outrdata)


  