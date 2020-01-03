# Jake Yeung
# Date of Creation: 2019-12-29
# File: ~/projects/scchic/scripts/rstudioserver_analysis/intestinal_all_merged/1-analyze_LDA_downstream.R
# Analyze LDA downstream

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(hash)
library(igraph)
library(umap)

library(scchicFuncs)

# Load LDA ----------------------------------------------------------------

prefix <- "Unenriched"
jmark <- "k4me1"
inf <- paste0("/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestines.2019-12-22/lda_outputs.mat.Scraped.", prefix, ".", jmark, ".TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.2019-12-22.K-30.binarize.FALSE/ldaOut.mat.Scraped.", prefix, ".", jmark, ".TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.2019-12-22.K-30.Robj")
assertthat::assert_that(file.exists(inf))
load(inf, v=T)
tm.result <- posterior(out.lda)

# jmark <- "k27me3"
jmark <- "k4me1"
jprefix <- "Unenriched"
inf.topics <- paste0("/home/jyeung/hpc/intestinal_scchic/LDA_outputs/gensim/intestines_all/intestines.Scraped.", prefix, ".", jmark, ".eta_auto.alpha_auto.ntopics_30.mapq_40.npasses_100.lda_model.doc_to_topics.csv")
inf.terms <- paste0("/home/jyeung/hpc/intestinal_scchic/LDA_outputs/gensim/intestines_all/intestines.Scraped.", prefix, ".", jmark, ".eta_auto.alpha_auto.ntopics_30.mapq_40.npasses_100.lda_model.topic_to_terms.csv")
inf.cnames <- paste0("/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/quality_controls_intestines.OK/mat.Scraped.", prefix, ".", jmark, ".TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.2019-12-22.colnames")
inf.rnames <- paste0("/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/quality_controls_intestines.OK/mat.Scraped.", prefix, ".", jmark, ".TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.2019-12-22.rownames")
tm.result <- GetTmResultFromGensim(inf.topics = inf.topics, inf.terms = inf.terms, inf.cellnames = inf.cnames, inf.binnames = inf.rnames)

# Run UMAP  ---------------------------------------------------------------


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(tm.result$topics, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(is.stem = grepl(pattern = "stemcells", x = cell),
         is.lgr5 = grepl(pattern = "Lgr5", x = cell),
         is.stem = is.stem | is.lgr5,
         experi = ClipLast(cell, jsep = "_"))
unique(dat.umap.long$experi)

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  # facet_wrap(~experi)




# Find celltypes  ---------------------------------------------------------





# Do varaince  ------------------------------------------------------------

dat.impute.log <- t(log2(tm.result$topics %*% tm.result$terms))
jchromos <- c(paste("chr", seq(19), sep = ""), "chrX", "chrY")

dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.merge <- left_join(dat.umap.long, dat.var)

ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~is.stem)

ggplot(dat.merge, aes(x = cell.var.within.sum.norm, fill = is.stem)) + geom_density(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

