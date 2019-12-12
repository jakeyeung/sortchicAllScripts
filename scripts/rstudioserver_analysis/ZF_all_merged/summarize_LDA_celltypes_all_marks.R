# Jake Yeung
# Date of Creation: 2019-12-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/summarize_LDA_celltypes_all_marks.R
# All marks celltypes 



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(scchicFuncs)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(ggrepel)

ProcessLDA <- function(out.lda, jmark){
  print(jmark)
  tm.result <- posterior(out.lda)
  topics.mat <- tm.result$topics
  colnames(topics.mat) <- paste("topic", colnames(topics.mat), sep = "_")
  
  # Plot UMAP  --------------------------------------------------------------
  
  umap.out <- umap(topics.mat, config = jsettings)
  
  dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(experi = ClipLast(cell))
  
  m.umap.split <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~experi) + ggtitle(jmark)
  
  dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
  
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  
  dat.merge <- left_join(dat.umap.long, dat.var)
  
  # Get entropy for each topics ---------------------------------------------
  
  # load bins
  terms.mat.tmp <- tm.result$terms
  # rownames(terms.mat.tmp) <- rownames(terms.mat.tmp)
  annot.out <- AnnotateBins2(terms.mat = tm.result$terms, top.thres = 0.995, inf.tss = inf.annot, txdb = TxDb.Drerio.UCSC.danRer11.refGene, annodb = "org.Dr.eg.db", chromos.keep = jchromos)
  annot.out$terms.annot <- annot.out$terms.annot %>%
    rowwise() %>%
    mutate(termgene = ifelse(is.na(termgene), "", termgene))
  annot.out$terms.annot$topic <- sapply(annot.out$terms.annot$topic, function(x) paste("topic_", x, sep = ""))
  annot.out$terms.annot$gene <- sapply(annot.out$terms.annot$termgene, function(x){
    jsplit <- strsplit(x, ";")[[1]]
    if (length(jsplit) > 1){
      return(jsplit[[2]])
    } else {
      return("")
    }
  })
  
  topics.sum <- OrderTopicsByEntropy(tm.result, jquantile = 0.99)
  topics.sum$topic <- gsub("^X", "topic_", topics.sum$topic)
  outlst <- list()
  return(list(tm.result = tm.result, annot.pout = annot.out, topics.sum = topics.sum, dat.merge = dat.merge))
}


# Load LDA ----------------------------------------------------------------


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks




# Constants ---------------------------------------------------------------



jbin <- "FALSE"
jsuff <- "AllMerged"
# outdir <- "/home/jyeung/data/from_rstudioserver/scchic/pdfs/"
outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/LDA_downstream_ZF"

jwin <- 100000L

zscore.cutoff <- "top_500"

# load Chloe data
jchromos <- paste("chr", seq(25), sep = "")
inf.annot <- paste0("/home/jyeung/hpc/databases/gene_tss/gene_tss.winsize_", jwin, ".species_drerio.bed")
assertthat::assert_that(file.exists(inf.annot))


normtype <- "_cpmnorm"
normtype2 <- "_removethrombo"
jdate <- "2019-12-10"

inf.WKM <- paste0("/home/jyeung/hpc/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/from_macbook/Baron_et_al_pseudobulk_Zebrafish_WKM", normtype, ".", jdate, ".rds")
assertthat::assert_that(file.exists(inf.WKM))

inf.WKMnothrombo <- paste0("/home/jyeung/hpc/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/from_macbook/Baron_et_al_pseudobulk_Zebrafish_WKM", normtype, ".", jdate, ".rds")
assertthat::assert_that(file.exists(inf.WKMnothrombo))

# UMAP settings
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# load tx data ------------------------------------------------------------

# from make_tx_dataset_zebrafish_WKM.R
dat.bulk <- readRDS(inf.WKM)
dat.bulk$celltype <- as.factor(dat.bulk$celltype)

dat.bulk.nothrombo <- readRDS(inf.WKMnothrombo)
dat.bulk.nothrombo$celltype <- as.factor(dat.bulk.nothrombo$celltype)


# Process -----------------------------------------------------------------

jbin <- TRUE
jmarks.act <- jmarks[1:2]
infs.act <- lapply(jmarks.act, function(jmark){
    inf <- paste0("/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_ZFWKM_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.ZF_", 
                  jsuff, "_", jmark, 
                  ".TAcutoff_0.5.countscutoff_500.binfilt_cellfilt.2019-12-05.K-30.binarize.", jbin, 
                  "/ldaOut.ZF_", jsuff, "_", jmark, ".TAcutoff_0.5.countscutoff_500.binfilt_cellfilt.2019-12-05.K-30.Robj")
    assertthat::assert_that(file.exists(inf))
    return(inf)
})

jbin <- FALSE
jmarks.rep <- jmarks[3:4]
infs.rep <- lapply(jmarks.rep, function(jmark){
  inf <- paste0("/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_ZFWKM_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.ZF_", 
                jsuff, "_", jmark, 
                ".TAcutoff_0.5.countscutoff_500.binfilt_cellfilt.2019-12-05.K-30.binarize.", jbin, 
                "/ldaOut.ZF_", jsuff, "_", jmark, ".TAcutoff_0.5.countscutoff_500.binfilt_cellfilt.2019-12-05.K-30.Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

infs <- c(infs.act, infs.rep)

out.ldas <- lapply(infs, function(inf){
  load(inf, v=T)
  return(out.lda)
})
 
system.time(
  dat.outputs <- lapply(jmarks, function(jmark) ProcessLDA(out.ldas[[jmark]], jmark))
)

outdir <- "~/data/from_rstudioserver/scchic/ZF_objs"
# dir.create(outdir)
save(dat.outputs, file = file.path(outdir, "dat_across_marks_outputs.binactive_nobinrepress.RData"))

# plot dat.outputs

# plot stem cell topic?
# sc.top.k4me1 <- dat.outputs$H3K4me1$topics.sum$topic




