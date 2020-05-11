# Jake Yeung
# Date of Creation: 2020-05-08
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/redo_multiomics_plot_genesets_on_LDA_or_GLMPA_umap.R
# 


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


library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)
jorg <- "org.Mm.eg.db"
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")




jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks


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

jwinsize <- "10000"

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


# go with the K4me3 definition...

ref.mark <- "H3K4me3"
jthres <- 275  # maybe not exactly at hump? what about tissuespecific stuff? rare celltypes? complicated from the bulk 
plot(density(tss.exprs.lst.unfilt[[ref.mark]]))
abline(v = jthres)

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

g2e.annot <- hash(annot.out$regions.annotated$SYMBOL, annot.out$regions.annotated$ENSEMBL)

plot(density(annot.regions.filt$dist.to.tss))

# Plot expression of different gene sets onto UMAP  -----------------------

# try for K4me3 first?

# get a gene set



jmark <- "H3K4me3"
jctype <- "Erythroblast"

jctypes <- names(de.ens.sorted.stringent)

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/stemcell_analysis/plot_exprs_genesets_UMAP_and_averages/genesets_umaps_and_avgs.", Sys.Date(), ".pdf")
pdf(file = outpdf, useDingbats = FALSE)
for (jctype in jctypes){
  for (jmark in jmarks){
    # jtitle <- c(paste(jmark, jctype), paste("Ngenes=", length(jgset.bins)))
    jtitle <- paste(jmark, jctype, "Ngenes=", length(jgset.bins))
    
    print(jtitle)
    
    jgset <- as.character(de.ens.sorted.stringent[[jctype]])
    jgset.bins <- subset(annot.regions.filt, ens %in% jgset)$region_coord
    # get average imputed exprs
    exprs.vec <- colMeans(dat.imputes.lst[[jmark]][jgset.bins, ])
    exprs.dat <- data.frame(cell = names(exprs.vec), exprs = exprs.vec, stringsAsFactors = FALSE)
    print(head(exprs.dat))
    
    jmerge <- left_join(out.objs[[jmark]]$dat.umap.glm.fillNAs, exprs.dat)
    jmerge.summary <- jmerge %>%
      group_by(cluster) %>%
      summarise(ncells = length(cell),
                exprs.mean = mean(exprs),
                exprs.sd = sd(exprs))
    
    m.umap <- ggplot(jmerge, aes(x = umap1, y = umap2, color = exprs)) + geom_point() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      scale_color_viridis_c() + 
      ggtitle(jtitle)
    m.umap.rev <- ggplot(jmerge, aes(x = umap1, y = umap2, color = exprs)) + geom_point() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      scale_color_viridis_c(direction = -1) + 
      ggtitle(jtitle)
    m.avg <- ggplot(jmerge.summary %>% 
                      filter(!is.na(cluster)), 
                    aes(x = forcats::fct_reorder(.f = cluster, .x = exprs.mean, .fun = median, .desc = TRUE), y = exprs.mean, ymin = exprs.mean - exprs.sd, ymax = exprs.mean + exprs.sd)) + 
      geom_point() + geom_errorbar() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
      xlab("") + ylab("Imputed scChIC signal") + ggtitle(jtitle)
    print(m.umap)
    print(m.umap.rev)
    print(m.avg)
  }
}
dev.off()



