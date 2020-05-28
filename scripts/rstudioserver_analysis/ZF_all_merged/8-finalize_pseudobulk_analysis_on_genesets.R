# Jake Yeung
# Date of Creation: 2020-05-01
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/8-finalize_pseudobulk_analysis_on_genesets.R
# Gene sets 

rm(list=ls())

jstart <- Sys.time()

set.seed(0)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(JFuncs)
library(scchicFuncs)

library(topicmodels)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)
library(GenomicRanges)


# write to output




# Constants ---------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
jvarcutoffs <- c(0.75, 2, 1, 0.5)
names(jvarcutoffs) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/pseudobulk_analysis.manual_DE_genes.fix_bin_to_gene_assignment"
outf <- file.path(outdir, paste0("finalize_geneset_analysis_pseudobulk_and_singlecells.", Sys.Date(), ".pdf"))

pdf(outf, useDingbats = FALSE)

# Load GLMPCA -------------------------------------------------------------

glm.imputes.lst <- lapply(jmarks, function(jmark){
  inf.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/GLMPCA_outputs.ZF_AllMerged.imputevarfilt.lessstringent/ZF_", jmark, ".AllMerged.ZF_AllMerged.imputevarfilt.lessstringent.GLMPCA_var_correction.mergebinsize_1000.binskeep_500.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.5.winsorize_TRUE.2020-04-29.RData")
  load(inf.glmpca, v=T)
  # do it on the GLMPCA? 
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- 0.1
  jsettings$random_state <- 123
  dat.umap.glm <- DoUmapAndLouvain(glm.out$factors, jsettings = jsettings)
  glm.impute <- t(as.matrix(glm.out$factors) %*% as.matrix(t(glm.out$loadings)))
  return(glm.impute)
})

# Load LDA  ---------------------------------------------------------------

inmain <- paste0(hubprefix, "/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZF_AllMerged.winsize_50000.imputevarfilt.lessstringent")
assertthat::assert_that(dir.exists(inmain))

tm.result.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jvarcutoff <- jvarcutoffs[[jmark]]
  infname <- paste0("lda_outputs.counts_table_var_filt.", jmark, ".imputevar_", jvarcutoff, ".K-30.binarize.FALSE/ldaOut.counts_table_var_filt.", jmark, ".imputevar_", jvarcutoff, ".K-30.Robj")
  inf <- file.path(inmain, infname)
  assertthat::assert_that(file.exists(inf))
  
  load(inf, v=T)
  count.mat <- as.matrix(count.mat)
  tm.result <- posterior(out.lda)
  colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = "_")
  rownames(tm.result$terms) <- paste("topic", rownames(tm.result$terms), sep = "_")
  return(tm.result)
})

# get imputes

dat.imputes.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  tm.result <- tm.result.lst[[jmark]]
  pmat <- t(tm.result$topics %*% tm.result$terms)
  dat.impute.log <- log2(pmat * 10^6)
  # pmat.odds <- pmat / (1 - pmat)
  # dat.impute.log <- log2(pmat.odds)
  # dat.impute.log2 <- log2(pmat)
  return(dat.impute.log)
})

# load annots 

annot.glmpca.filt.lst <- lapply(jmarks, function(jmark){
  
  # filter by previously defined celltypes? 
  inf.annot.louv <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/LDA_downstream/LDA_downstream_ZF.2020-04-23.imputevarfilt.lessstringent/ZF_LDA_output.", jmark, ".keepn_150.final.ClusterTables.txt")
  assertthat::assert_that(file.exists(inf.annot.louv))
  
  
  inf.annot.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt")
  assertthat::assert_that(file.exists(inf.annot.glmpca))
  
  annot.louv <- fread(inf.annot.louv)
  annot.louv$clusterplate <- paste(annot.louv$cluster, annot.louv$plate, "_")
  
  annot.glmpca <- fread(inf.annot.glmpca)
  annot.glmpca.filt <- subset(annot.glmpca, cell %in% annot.louv$cell) %>%
    rowwise() %>%
    mutate(clusterplate = paste(cluster, plate, sep = "_")) %>%
    mutate(cluster = ifelse(cluster %in% c("lymph1", "lymph2"), "lymph", cluster)) %>%   # rename lymph1 and lymph2 into lymph
    ungroup() %>%
    filter(cluster != "Unknown")  # remove the small cluster Unknown
  annot.glmpca.filt <- left_join(annot.glmpca.filt, subset(annot.louv, select = c(cell, var.imputed)))
  return(annot.glmpca.filt)
  
})


# Annotate bins -----------------------------------------------------------


# take any mark
inf.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/zebrafish/gene_tss.CodingOnly.winsize_10000.species_drerio.bed"
assertthat::assert_that(file.exists(inf.annot))
jchromos <- paste("chr", seq(25), sep = "")
annot.out <- AnnotateCoordsFromList.GeneWise(rownames(dat.imputes.lst[[1]]), inf.tss = inf.annot, txdb = TxDb.Drerio.UCSC.danRer11.refGene, annodb = "org.Dr.eg.db", chromos.keep = jchromos)
# annot.out <- AnnotateCoordsFromList.GeneWise2(rownames(dat.imputes.lst[[1]]), inf.tss = inf.annot, txdb = TxDb.Drerio.UCSC.danRer11.refGene, annodb = "org.Dr.eg.db", chromos.keep = jchromos)

# where is chr11:25408153-25418153;gata1a;1

jsub <- subset(annot.out$out2.df, seqnames == "chr11" & start > 25000000 & end < 25500000)


annot.regions <- annot.out$out2.df

annot.regions <- subset(annot.regions, select = c(dist.to.tss, region_coord, gene, tssname))

# add ensembl
g2e.dat <- data.frame(gene = as.character(annot.out$regions.annotated$SYMBOL), ens = as.character(annot.out$regions.annotated$ENSEMBL), stringsAsFactors = FALSE)
annot.regions <- left_join(annot.regions, g2e.dat)

# annot.regions.sub <- subset(annot.regions, distanceToTSS < 10000)
# why NAs?

genes.na <- unique(subset(annot.regions, is.na(ens))$gene)
ens.na <- Gene2Ensembl.ZF(genes.na)
jhash.na <- hash::hash(genes.na, ens.na)

annot.regions <- annot.regions %>%
  rowwise() %>%
  mutate(ens = ifelse(is.na(ens), AssignHash(gene, jhash.na), ens))

# annot.regions$ens <- sapply(annot.r)

# fill the NAs 


# Load pseudobulks ChIC  -------------------------------------------------------


pbulk.long.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.pbulk.chic <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/rdata_robjs/WKM_scChICseq_pseudobulk_downsampled.ens.", jmark, ".RData")
  load(inf.pbulk.chic, v=T)
  return(pbulk.long)
})

pbulk.chic.mat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.pbulk.chic <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/rdata_robjs/WKM_scChICseq_pseudobulk_downsampled.ens.", jmark, ".RData")
  load(inf.pbulk.chic, v=T)
  return(pbulk.filt.ds)
})


# Create one dataframe ----------------------------------------------------

ctypes.keep <- c("eryth", "HSC", "lymph", "monocyte")

# ctypes: eryth, HSC, lymph, monocytes
pbulk.k4me1 <- subset(pbulk.long.lst$H3K4me1) %>%
  ungroup() %>%
  mutate(mark = "H3K4me1") %>%
  filter(pbulk %in% ctypes.keep)

pbulk.k4me3 <- subset(pbulk.long.lst$H3K4me3) %>%
  ungroup() %>%
  mutate(mark = "H3K4me3") %>%
  filter(pbulk %in% ctypes.keep)

pbulk.k27me3 <- subset(pbulk.long.lst$H3K27me3) %>%
  ungroup() %>%
  mutate(mark = "H3K27me3",
         pbulk = gsub("eryth2", "eryth", pbulk),
         pbulk = gsub("HSC2", "HSC", pbulk)) %>%
  filter(pbulk %in% ctypes.keep)

# combine it all?
pbulk.merge <- bind_rows(pbulk.k4me1, pbulk.k4me3, pbulk.k27me3) %>%
  mutate(pbulk = gsub("monocyte", "granu", pbulk))

# recalculate log2FC and zscores...
pbulk.merge <- pbulk.merge %>%
  group_by(gene, ens, mark) %>%
  mutate(log2FC = log2cuts - mean(log2cuts),
         log2zscore = log2FC / sd(log2cuts))

pbulk.merge$mark <- factor(pbulk.merge$mark, levels = c("H3K4me1", "H3K4me3", "H3K27me3"))

ggplot(subset(pbulk.merge), aes(x = pbulk, y = log2cuts, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(pbulk.merge), aes(x = pbulk, y = log2FC, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(pbulk.merge), aes(x = pbulk, y = log2FC, fill = mark)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(pbulk.merge), aes(x = pbulk, y = log2zscore, fill = mark)) + 
  geom_boxplot() + geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Map gene name to ensembl ------------------------------------------------

g2e <- hash(pbulk.merge$gene, pbulk.merge$ens)  # lots of duplicates but should be OK?


# Define sets of genes from literature ------------------------------------

genesets <- list()

inf.kob.erythmyeloid <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Kobayashi_2019_zf_genelists/erythmyeloid_genes_list.txt"
kobayashi.genes.erythmyeloid <- fread(inf.kob.erythmyeloid, header = FALSE)$V1
jens.erythmyeloid <- kobayashi.genes.erythmyeloid

inf.kob <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Kobayashi_2019_zf_genelists/hspc_genes_list.txt"
kobayashi.genes <- fread(inf.kob, header = FALSE)$V1
jens.hspc <- kobayashi.genes
# HSPCs
jgenes.choose <- c("meis1b", "pmp22b", "ahnak", "krt8", "anxa5b", "mrc1b", "pa2g4a", "npm1a", "ahcy", "adh5", "fabp3", "myb")
# add kobayashi?
jgenes.choose <- c(jgenes.choose, kobayashi.genes)
jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)

jens.choose <- jens.choose[which(!is.na(jens.choose))]

genesets[["HSC"]] <- jens.choose

# lymphs
jgenes.choose <- c("pax5", "cd79a", "bhlhe40", "cd83", "cxcr4a", "cd74b", "cd74a", "cd37", "zfp36l1a")  # lymphs
jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)

genesets[["lymph"]] <- jens.choose

jens.choose <- jens.choose[which(!is.na(jens.choose))]

# monocytes
jgenes.choose1 <- c("adam8a", "odc1", "lta4h", "thy1", "scpp8", "illr4", "timp2b", "mmp9", "mmp13a", "scinlb")  # monocytes
jgenes.choose2 <- c("cpa5", "lyz", "lect2l", "npsn", "sms", "abcb9", "ch25hl2", "papss2b", "hsd3b7", "cfd")  # neutros
jgenes.choose <- c(jgenes.choose1, jgenes.choose2)
jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)

genesets[["granu"]] <- jens.choose

jens.choose <- jens.choose[which(!is.na(jens.choose))]

# eryths
jgenes.choose <- c("rhag", "prdx2", "epor", "gata1a", "tspo", "slc4a1a", "sptb", "cahz", "hbba1", "hbba2", "alas2", "epb41b", "nt5c2l1")
jens.choose <- sapply(jgenes.choose, AssignHash, g2e, null.fill = NA)

# ba1 is hbb-like, maybe another name? chr3:55,098,010-55,099,718 use hbba1 and hbba2 instead
# jsubset(annot.out$out2.df, seqnames == "chr3" & start > 55000000 & end < 55500000)

genesets[["eryth"]] <- jens.choose

jens.choose <- jens.choose[which(!is.na(jens.choose))]

# add a random set?

genesets[["random"]] <- sample(unique(pbulk.merge$ens), size = 500)

# Show genesets -----------------------------------------------------------

lapply(genesets, length)


# Show boxplots for gene sets ---------------------------------------------

jctypes <- names(genesets)

for (jctype in jctypes){
  jens.choose <- genesets[[jctype]]
  
  jtitle <- paste("Geneset:", jctype, "Ngenes:", length(jens.choose))
  
  m.exprs <- ggplot(subset(pbulk.merge, ens %in% jens.choose), aes(x = pbulk, y = log2cuts, fill = mark)) + 
    geom_boxplot() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  
  m.log2fc <- ggplot(subset(pbulk.merge, ens %in% jens.choose), aes(x = pbulk, y = log2FC, fill = mark)) + 
    geom_boxplot() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  
  m.zscore <- ggplot(subset(pbulk.merge, ens %in% jens.choose), aes(x = pbulk, y = log2zscore, fill = mark)) + 
    geom_boxplot() + geom_hline(yintercept = 0, linetype = "dotted") + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle) + ylab("Zscore")
  
  print(m.exprs)
  print(m.log2fc)
  print(m.zscore)
}


# Show DE genes in GLMPCA or LDA  -----------------------------------------

# check single cells

jmark.tmp <- "H3K4me3"
jclst <- "HSC"
jclst <- "lymph"
jclst <- "granu"
jclst <- "eryth"

jens.choose <- jens.erythmyeloid
print(names(genesets))

jclsts <- names(genesets)

rnames.glm <- rownames(glm.imputes.lst[[jmark.tmp]])

for (jclst in jclsts){
  for (jmark.tmp in jmarks){
    (jens.choose <- genesets[[jclst]])
    
    jsub <- subset(annot.regions, ens %in% jens.choose)
    
    regions_coord.matched <- unique(jsub$region_coord)
    ngenes <- length(regions_coord.matched)
    regions_coord.matched.glm <- which(rnames.glm %in% unique(jsub$region_coord))
    ngenes.glm <- length(regions_coord.matched.glm)
    
    exprs.dat.lda <- data.frame(cell = colnames(dat.imputes.lst[[jmark.tmp]]), exprs = colMeans(dat.imputes.lst[[jmark.tmp]][regions_coord.matched, ]))
    exprs.dat.glm <- data.frame(cell = colnames(glm.imputes.lst[[jmark.tmp]]), exprs = colMeans(glm.imputes.lst[[jmark.tmp]][regions_coord.matched.glm, ]))
    
    # plojt on UMAP 
    jmerge.lda <- left_join(annot.glmpca.filt.lst[[jmark.tmp]], exprs.dat.lda)
    jmerge.glm <- left_join(annot.glmpca.filt.lst[[jmark.tmp]], exprs.dat.glm)
    
    # plot the sumaries
    jmerge.lda.sum <- jmerge.lda %>%
      group_by(cluster) %>%
      summarise(ncells = length(cell),
                exprs.mean = mean(exprs),
                exprs.sd = sd(exprs))
    
    jmerge.glm.sum <- jmerge.glm %>%
      group_by(cluster) %>%
      summarise(ncells = length(cell),
                exprs.mean = mean(exprs),
                exprs.sd = sd(exprs))
    
    jtitle1 <- paste(jmark.tmp, "LDA:", jclst, "ngenes:", ngenes)
    jtitle2 <- paste(jmark.tmp, "GLMPCA:", jclst, "ngenes:", ngenes)
    
    m.lda <- ggplot(jmerge.lda, aes(x = umap1, y = umap2, color = exprs)) + 
      geom_point() + theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
      scale_color_viridis_c() + ggtitle(jtitle1)
    print(m.lda)
    m.glmpca <- ggplot(jmerge.glm, aes(x = umap1, y = umap2, color = exprs)) + geom_point() + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
      scale_color_viridis_c() + ggtitle(jtitle2)
    print(m.glmpca)
    
    m.lda.rev <- ggplot(jmerge.lda, aes(x = umap1, y = umap2, color = exprs)) + 
      geom_point() + theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
      scale_color_viridis_c(direction = -1) + ggtitle(jtitle1)
    print(m.lda.rev)
    m.glmpca.rev <- ggplot(jmerge.glm, aes(x = umap1, y = umap2, color = exprs)) + geom_point() + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
      scale_color_viridis_c(direction = -1) + ggtitle(jtitle2)
    print(m.glmpca.rev)
    
    m.avg1 <- ggplot(jmerge.lda.sum %>% 
                      filter(!is.na(cluster)), 
                    aes(x = forcats::fct_reorder(.f = cluster, .x = exprs.mean, .fun = median, .desc = TRUE), y = exprs.mean, ymin = exprs.mean - exprs.sd, ymax = exprs.mean + exprs.sd)) + 
      geom_point() + geom_errorbar() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
      xlab("") + ylab("Imputed scChIC signal") + ggtitle(jtitle1)
    
    m.avg2 <- ggplot(jmerge.glm.sum %>% 
                      filter(!is.na(cluster)), 
                    aes(x = forcats::fct_reorder(.f = cluster, .x = exprs.mean, .fun = median, .desc = TRUE), y = exprs.mean, ymin = exprs.mean - exprs.sd, ymax = exprs.mean + exprs.sd)) + 
      geom_point() + geom_errorbar() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
      xlab("") + ylab("Imputed scChIC signal") + ggtitle(jtitle2)
    print(m.avg1)
    print(m.avg2)
  }
}
dev.off()

# m <- ggplot(jmerge.lda %>% filter(cluster!="eryth"), aes(x = umap1, y = umap2, color = exprs)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_viridis_c() + ggtitle(jmark.tmp, paste(jclst, "ngenes:", ngenes))
# print(m)



print(Sys.time() - jstart)

