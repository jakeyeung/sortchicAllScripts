# Jake Yeung
# Date of Creation: 2020-04-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/4-GLMPCA_check_promoters_enhancers.H3K4me1.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)

library(hash)

library(hash)
library(igraph)
library(umap)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)
library(GenomicRanges)


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"


jmark <- "H3K4me1"
jwin <- 50000L
jchromos <- paste("chr", seq(25), sep = "")
keepn <- 150
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/celltyping/from_louvain"
outf <- file.path(outdir, paste0("cell_to_cluster_table.", jmark, ".", Sys.Date(), ".txt"))
outpdf <- file.path(outdir, paste0("cell_to_cluster_table.", jmark, ".", Sys.Date(), ".pdf"))

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")



# Load pseudobulk data  ---------------------------------------------------

normtype <- "_cpmnorm"
jdate <- "2019-12-10"
inf.WKM <- file.path(hubprefix, paste0("jyeung/data/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/from_macbook/Baron_et_al_pseudobulk_Zebrafish_WKM", normtype, ".", jdate, ".rds"))
assertthat::assert_that(file.exists(inf.WKM))

# from make_tx_dataset_zebrafish_WKM.R
dat.bulk <- readRDS(inf.WKM)
dat.bulk$celltype <- as.factor(dat.bulk$celltype)


# Load GLMPCA make final table for clustering -----------------------------


# Load GLMPCA output ------------------------------------------------------

inmain <- paste0(hubprefix, "/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZF_AllMerged.winsize_50000")
assertthat::assert_that(dir.exists(inmain))

inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/GLMPCA_outputs.ZF_AllMerged/ZF_", jmark, ".AllMerged.ZF_AllMerged.GLMPCA_var_correction.mergebinsize_1000.binskeep_500.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.5.winsorize_TRUE.2020-04-13.RData")
assertthat::assert_that(file.exists(inf))
load(inf, v=T)


tm.result.glm <- list(topics = glm.out$factors, terms = t(glm.out$loadings))
dat.impute.glm <- t(as.matrix(tm.result.glm$topics) %*% as.matrix(tm.result.glm$terms))

# Load LDA ----------------------------------------------------------------

infname <- paste0("lda_outputs.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.Robj")
inf <- file.path(inmain, infname)
assertthat::assert_that(file.exists(inf))

load(inf, v=T)

tm.result.lda <- posterior(out.lda)

tm.result.lda <- AddTopicToTmResult(tm.result.lda, jsep = "")

dat.umap.before <- DoUmapAndLouvain(tm.result.lda$topics, jsettings) %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_")) %>%
  dplyr::rename(louvain.before = louvain)

dat.umap <- DoUmapAndLouvain(glm.out$factors, jsettings) %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_")) %>%
  left_join(., subset(dat.umap.before, select = c(cell, louvain.before)))

m.before <- ggplot(dat.umap.before, aes(x = umap1, y = umap2, color = louvain.before)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = cbPalette) 

m.after.colorbybefore <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain.before)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = cbPalette) 

m.after <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = cbPalette) 

m.after.plates <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)  + facet_wrap(~plate)


pdf(outpdf, width = 1020/72, height = 815/72, useDingbats = FALSE)

JFuncs::multiplot(m.before, m.after.colorbybefore, m.after, cols = 3)
print(m.after.plates)
# Annotate bins -----------------------------------------------------------

inf.annot <- file.path(hubprefix, paste0("jyeung/data/databases/gene_tss/gene_tss.winsize_", jwin, ".species_drerio.bed"))
assertthat::assert_that(file.exists(inf.annot))


annot.out <- AnnotateBins2(terms.mat = tm.result.lda$terms, top.thres = 0.995, inf.tss = inf.annot, txdb = TxDb.Drerio.UCSC.danRer11.refGene, annodb = "org.Dr.eg.db", chromos.keep = jchromos, skip.split = TRUE)
annot.out$terms.annot <- annot.out$terms.annot %>%
  rowwise() %>%
  mutate(termgene = ifelse(is.na(termgene), "", termgene))
# annot.out$terms.annot$topic <- sapply(annot.out$terms.annot$topic, function(x) paste("topic_", x, sep = ""))
annot.out$terms.annot$gene <- sapply(annot.out$terms.annot$termgene, function(x){
  jsplit <- strsplit(x, ";")[[1]]
  if (length(jsplit) > 1){
    return(jsplit[[2]])
  } else {
    return("")
  }
})


annot.out$terms.annot <- annot.out$terms.annot %>% 
  group_by(topic) %>%
  mutate(rnk.rev = rank(weight))

# maybe 0 distance is too stringent?? 

print(head(annot.out$terms.annot))

jcoord <- "chr5:70750000-70800000"
subset(annot.out$regions.annotated, region_coord == jcoord)


# Check celltypes ---------------------------------------------------------

jtopics.lst <- list("monocytes" = c("topic16"),
                    "eryth" = c("topic14"),
                    "lymph" = c("topic11", "topic29"),
                    "HSCs" = c("topic24"),
                    "thrombocytes" = c("topic8"))

# name louvains to a celltype
louv2ctype <- hash("louvain1" = "lymph1",
                   "louvain2" = "lymph2",
                   "louvain3" = "monocyte",
                   "louvain4" = "Unknown",
                   "louvain5" = "eryth",
                   "louvain6" = "thrombo",
                   "louvain7" = "HSC")

merged.coords.exprs.lst <- lapply(names(jtopics.lst), function(jname){
  print(jname)
  jtopic <- jtopics.lst[[jname]]
  
  # show genes in monocytes
  jsub.terms <- subset(annot.out$terms.annot, topic %in% jtopic) %>%
    dplyr::filter(rnk <= keepn)
  
  jgenes <- jsub.terms$gene
  jcoords <- jsub.terms$term
  
  merged.coords.exprs <- data.frame(cell = colnames(dat.impute.glm), exprs = colSums(dat.impute.glm[jcoords, ]), stringsAsFactors = FALSE)
  colnames(merged.coords.exprs)[2] <- jname
  
  # plot on UMAP 
  dat.merge <- left_join(dat.umap, merged.coords.exprs)
  
  m <- ggplot(dat.merge, aes_string(x = "umap1", y = "umap2", color = jname)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c() + ggtitle(paste(jname, jtopic))
  print(m)
  return(merged.coords.exprs)
})

# Write output assigning cell to cluster ----------------------------------

print(m.after.plates)


dat.umap.pretty <- dat.umap %>%
  mutate(louvain = paste("louvain", louvain, sep = ""))


dat.umap.pretty$cluster <- sapply(dat.umap.pretty$louvain, function(l) louv2ctype[[l]])

# write to output
dat.umap.pretty <- dat.umap.pretty %>%
  dplyr::select(cell, cluster, umap1, umap2, louvain, plate)

# plot pretty
m.final <- ggplot(dat.umap.pretty, aes(x = umap1, y = umap2, color = cluster)) + geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + 
  facet_wrap(~plate)
print(m.final)
dev.off()


fwrite(dat.umap.pretty, file = outf, sep = "\t", row.names = FALSE, col.names = TRUE)



