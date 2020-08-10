# Jake Yeung
# Date of Creation: 2020-06-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/find_TSS_topics_winsizes_for_heatmap_MouseBM.R
# description


rm(list=ls())

library(hash)
library(ggrastr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(preprocessCore)

library(mixtools)

library(scchicFuncs)
library(JFuncs)

library(topicmodels)


library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)
jorg <- "org.Mm.eg.db"
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")




jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

jwinsize <- "10000"
ref.mark <- "H3K4me1"
topnbins <- 2000

# Load DE genes -----------------------------------------------------------

# load this first because it loads a lot of objects, might disuprt things

inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData"
load(inf.de, v=T)


# Load LDA r GLMPCA ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden/jyeung/data"

# load GLMPCA from bins 
# jmark <- "H3K4me1"

jexperi <- "AllMerged"
mergesize <- "1000"
nbins <- "1000"
jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
jpenalty <- 1
ntopics <- 30


out.objs <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.glmpca <- file.path(hubprefix, paste0("scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.good_runs/PZ_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_1000.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.winsorize_TRUE.2020-02-11.RData"))
  inf.lda <- file.path(hubprefix, paste0("scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData"))
  inf.lda.bins <- file.path(hubprefix, paste0("scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj"))
  load(inf.glmpca, v=T)
  load(inf.lda, v=T)
  load(inf.lda.bins, v=T)
  
  out <- list(dat.umap.glm.fillNAs = dat.umap.glm.fillNAs, dat.umap.lda = dat.umap.lda, glm.out = glm.out, out.lda = out.lda)
  return(out)
})

jbins <- out.objs$H3K4me1$out.lda@terms

# get imputed mats

dat.imputes.lst <- lapply(out.objs, function(x){
  tm.result <- topicmodels::posterior(x$out.lda)
  dat.impute <- log2(t(tm.result$topics %*% tm.result$terms) * 10^6)
  return(dat.impute)
})


# Read TSS Signal to figure out which transcript to keep  -----------------


indir.tss <- file.path(hubprefix, "scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS")
assertthat::assert_that(dir.exists(indir.tss))

tss.out <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.tss <- file.path(indir.tss, paste0(jmark, ".countTableTSS.mapq_40.TSS_", jwinsize, ".blfiltered.csv"))
  mat.tss <- ReadMatTSSFormat(inf.tss)
  return(list(mat.tss = mat.tss, tss.exprs = rowSums(mat.tss)))
})

tss.exprs.lst.unfilt <- lapply(tss.out, function(x) x$tss.exprs)
tss.mats.singlecell.unfilt <- lapply(tss.out, function(x) x$mat.tss)


# exprs.vec <- tss.exprs.lst$H3K4me1
lapply(jmarks, function(jmark){
  plot(density(tss.exprs.lst.unfilt[[jmark]]), main = jmark)
})


tss.mat.ref <- CollapseRowsByGene(count.mat = tss.mats.singlecell.unfilt[[ref.mark]], as.long = FALSE, track.kept.gene = TRUE)
tss.keep <- rownames(tss.mat.ref)

tss.exprs.lst <- lapply(tss.exprs.lst.unfilt, function(exprs.vec){
  jkeep <- names(exprs.vec) %in% tss.keep
  return(exprs.vec[jkeep])
})

print("Dimensions of TSS raw keeping all TSS")
lapply(tss.mats.singlecell.unfilt, dim)
tss.mats.singlecell <- lapply(tss.mats.singlecell.unfilt, function(tss.mat){
  jkeep <- rownames(tss.mat) %in% tss.keep
  return(tss.mat[jkeep, ])
})

print("Dimensions of TSS after keeping one TSS for each gene, defined by highest expression in H3K4me3")
lapply(tss.mats.singlecell, dim)

# Get common rows ---------------------------------------------------------

lapply(tss.exprs.lst.unfilt, length)

tss.all <- lapply(tss.exprs.lst, function(exprs.lst){
  names(exprs.lst)
}) %>%
  unlist() %>%
  unique()

tss.common <- lapply(tss.exprs.lst, function(exprs.lst){
  names(exprs.lst) 
}) %>%
  Reduce(f = intersect, .)

# get ensembl names ? 
genes.common <- sapply(tss.common, function(x) strsplit(x, ";")[[1]][[2]])
ens.common <- Gene2Ensembl.ZF(genes.common, return.original = TRUE, species = "mmusculus")

g2e.hash2 <- hash(genes.common, ens.common)

# create tss, genes, ens dat
genes.annot <- data.frame(bin = tss.common, gene = genes.common, ens = ens.common, stringsAsFactors = FALSE)




# Annotate bins to gene  --------------------------------------------------

# use same winsize (10kb as the TSS analysis)
# take any mark
# inf.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/zebrafish/gene_tss.CodingOnly.winsize_10000.species_drerio.bed"
inf.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.", jwinsize, ".bed")
assertthat::assert_that(file.exists(inf.annot))
annot.out <- AnnotateCoordsFromList.GeneWise(coords.vec = jbins, inf.tss = inf.annot, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = jorg, chromos.keep = jchromos)

annot.regions <- annot.out$out2.df
annot.regions <- subset(annot.regions, select = c(dist.to.tss, region_coord, gene, tssname))


# Filter bins for only TSS's that are good  -------------------------------

annot.regions.filt <- subset(annot.regions, tssname %in% tss.common)
annot.regions.filt$ens <- sapply(annot.regions.filt$gene, function(g) AssignHash(g, jhash = g2e.hash2, null.fill = g))

print(head(annot.regions.filt))

# g2e.annot <- hash(annot.out$regions.annotated$SYMBOL, annot.out$regions.annotated$ENSEMBL)
g2e.annot <- hash(annot.regions.filt$gene, annot.regions.filt$ens)
r2g.annot <- hash(annot.regions.filt$region_coord, annot.regions.filt$gene)
g2tss.annot <- hash(genes.annot$gene, genes.annot$bin)

plot(density(annot.regions.filt$dist.to.tss))


# Find celltype-specific topics  --------------------------------------------


out.lda <- out.objs[[ref.mark]]$out.lda


cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(out.objs$H3K4me1$dat.umap.glm.fillNAs, aes(x = umap1, y = umap2, color = cluster)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)


# browse /hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/BM_LDA_downstream_topics_celltypes_Giladi.UnenrichedAllMerged.KeepBestPlates2

# ertryth, bcell, granu, hsc for H3K4me3
ctypes.vec <- c("Eryth", "Bcell", "Granu", "HSPCs")
topics.vec <- c("topic13", "topic10", "topic23", "topic7")
names(topics.vec) <- ctypes.vec

# # check topics are correct
# jctype.check <- "Eryth"
# jctype.check <- "Bcell"
# jctype.check <- "Granu"
# jctype.check <- "HSPCs"
# jtop.check <- topics.vec[[jctype.check]]
# dat.umap.glm.check <- out.objs$H3K27me3$dat.umap.glm.fillNAs
# ctypes.check <- names(sort(tm.result$topics[, jtop.check], decreasing = TRUE)[1:100])
# dat.umap.glm.check$cellcheck <- sapply(dat.umap.glm.check$cell, function(x) x %in% ctypes.check)
# ggplot(dat.umap.glm.check, aes(x = umap1, y = umap2, color = cellcheck)) + 
#   geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   ggtitle(paste(jctype.check, jtop.check))


tm.result <- posterior(out.lda)
tm.result <- AddTopicToTmResult(tm.result, jsep = "")

# Convert topic regions to genes  -----------------------------------------

topbins.lst <- lapply(topics.vec, function(jtop){
  jvec <- sort(tm.result$terms[jtop, ], decreasing = TRUE)
  return(names(jvec)[1:topnbins])
})

topgenes.lst <- lapply(topbins.lst, function(jbins){
  jvec <- sapply(jbins, AssignHash, r2g.annot)
  jvec <- gsub(pattern = "Hoxa11", replacement = "Hoxa9", jvec)
  return(jvec)
})

toptss.lst <- lapply(topgenes.lst, function(jgenes){
  tss <- sapply(jgenes, AssignHash, g2tss.annot)
  # remove NA
  tss <- tss[which(!is.na(tss))]
})

lapply(toptss.lst, length)

# Make TSS into 2bp bins  -------------------------------------------------

print(head(toptss.lst$Eryth))

tss.bed.lst <- lapply(toptss.lst, function(tss.vec){
  jcoords <- sapply(tss.vec, function(x) strsplit(x, ";")[[1]][[1]])
  jtx <- sapply(tss.vec, function(x) strsplit(x, ";")[[1]][[2]])
  bed.tmp <- data.frame(chromo = sapply(jcoords, GetChromo), 
                        Start = as.numeric(sapply(jcoords, GetStart)),
                        End = as.numeric(sapply(jcoords, GetEnd)), 
                        tx = jtx, 
                        stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(midpt = (Start + End) / 2,
           Start2 = midpt - 1,
           End2 = midpt + 1)
  return(bed.tmp)
})



# Write to output ---------------------------------------------------------

outdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations/MouseBMFromTopics.refmark_", ref.mark, ".", topnbins)
dir.create(outdir)
assertthat::assert_that(dir.exists(outdir))


for (ctype in ctypes.vec){
  print(ctype)
  fname <- paste0("MouseBM_TSS_FromTopics.refmark_", ref.mark, ".", ctype, ".bsize_2.bed")
  outf <- file.path(outdir, fname)
  outdat <- tss.bed.lst[[ctype]] %>%
    dplyr::select(chromo, Start2, End2, tx)
  print(head(outdat))
  print(outf)
  fwrite(outdat, file = outf, sep = "\t", col.names = FALSE, na = "NA", quote = FALSE)
}

