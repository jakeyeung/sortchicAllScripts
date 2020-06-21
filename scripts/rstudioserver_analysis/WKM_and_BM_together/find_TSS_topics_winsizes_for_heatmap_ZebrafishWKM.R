# Jake Yeung
# Date of Creation: 2020-06-18
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/find_TSS_topics_winsizes_for_heatmap_ZebrafishWKM.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(scchicFuncs)
library(JFuncs)

# Contants ----------------------------------------------------------------

# ref.mark <- "H3K4me3"
ref.mark <- "H3K4me1"
ctypes <- c("HSPCs", "Bcells", "Granus", "Eryths")
names(ctypes) <- ctypes

# # H3K4me1
# jtops.lst <- list(c("topic26"), c("topic13", "topic22"), c("topic17"), c("topic25"))

# H3K4me3

if (ref.mark == "H3K4me3"){
  (jtops.lst <- list(c("topic1"), c("topic24", "topic19"), c("topic16"), c("topic18")))
  # (jtops.lst <- list(c("topic1"), c("topic24", "topic19"), "topic16", "topic18"))
} else if (ref.mark == "H3K4me1"){
  jtops.lst <- list(c("topic26"), c("topic13", "topic22"), c("topic17"), c("topic25"))
} else {
  stop(ref.mark, "not yet coded")
}

names(jtops.lst) <- ctypes

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jvarcutoffs <- c(0.75, 2, 1, 0.5)
names(jvarcutoffs) <- jmarks
jmark <- "H3K4me1"

# topnbins <- 500
topnbins.vec <- c(500, 1000, 2000, 5000)

# if (ref.mark == "H3K4me3"){
#   outmain <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations"
# } else {
#   outmain <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations/ZFrefmark_", ref.mark)
#   dir.create(outmain)
# }
outmain <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations"


hubprefix <- "/home/jyeung/hub_oudenaarden"
for (topnbins in topnbins.vec){
  
  print(topnbins)
  
  if (ref.mark == "H3K4me3"){
    outdir <- file.path(outmain, paste0("ZebrafishWKMFromTopics.", topnbins))
  } else {
    outdir <- file.path(outmain, paste0("WKMrefmark_", ref.mark, ".FromTopics.", topnbins))
  }
  dir.create(outdir)
  
  
  out.objs <- lapply(jmarks, function(jmark){
    
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
    
    # plot uMAP 
    inmain <- paste0(hubprefix, "/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZF_AllMerged.winsize_50000.imputevarfilt.lessstringent")
    assertthat::assert_that(dir.exists(inmain))
    
    jvarcutoff <- jvarcutoffs[[jmark]]
    infname <- paste0("lda_outputs.counts_table_var_filt.", jmark, ".imputevar_", jvarcutoff, ".K-30.binarize.FALSE/ldaOut.counts_table_var_filt.", jmark, ".imputevar_", jvarcutoff, ".K-30.Robj")
    # infname <- paste0("lda_outputs.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.Robj")
    inf <- file.path(inmain, infname)
    assertthat::assert_that(file.exists(inf))
    
    annot.louv.annot <- left_join(annot.louv, subset(annot.glmpca.filt, select = c(cell, cluster)), by = "cell")
    
    load(inf, v=T)
    out <- list(annot.glmpca.filt = annot.glmpca.filt, annot.louv = annot.louv.annot, out.lda = out.lda)
  })
  
  
  
  # Annotate bins -----------------------------------------------------------
  
  
  
  # annotate bins
  
  # inf.annot1 <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/zebrafish/gene_tss.CodingOnly.winsize_20000.species_drerio.bed"
  inf.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/zebrafish.SelectedTSSfromWKM/gene_tss.SelectedFromWKM.winsize_50000.species_drerio.bed"
  assertthat::assert_that(file.exists(inf.annot))
  
  jchromos <- paste("chr", seq(25), sep = "")
  
  rnames.common <- Reduce(intersect, lapply(out.objs, function(x) x$out.lda@terms))
  # rnames.all <- lapply(out.objs, function(x) x$out.lda@terms)
  
  tss.dat <- fread(inf.annot, col.names = c("seqnames", "start", "end", "tssname"))
  tss.dat$gene <- sapply(tss.dat$tssname, function(x) strsplit(x, ";")[[1]][[2]])
  
  # tss.dat1 <- fread(inf.annot1, col.names = c("seqnames", "start", "end", "tssname"))
  # tss.dat1$gene <- sapply(tss.dat1$tssname, function(x) strsplit(x, ";")[[1]][[2]])
  # annots.tss.gr <- makeGRangesFromDataFrame(tss.dat, keep.extra.columns = TRUE)
  # annots.tss.gr <- makeGRangesFromDataFrame(tss.dat1, keep.extra.columns = TRUE)
  
  
  annot.out <- AnnotateCoordsFromList(rnames.common, inf.tss = inf.annot, txdb = TxDb.Drerio.UCSC.danRer11.refGene, annodb = "org.Dr.eg.db", chromos.keep = jchromos)
  
  # assign regions to either promoter or not 
  
  # creaete hash or coord to distance??
  
  # dhash <- hash(annot.out$regions.annotated$region_coord, annot.out$regions.annotated$distanceToTSS)
  
  r2g <- hash(annot.out$out2.df.closest$region_coord, annot.out$out2.df.closest$gene)
  g2tss <- hash(annot.out$out2.df.closest$gene, annot.out$out2.df.closest$tssname)
  
  
  # Load zebrafish TSS  -----------------------------------------------------
  
  
  jdate <- "2020-06-09"
  indir.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/setup_matrix_for_poisson_regression.likeBM.redo_count_tables"
  assertthat::assert_that(dir.exists(indir.tss))
  
  jprefix <- file.path(indir.tss, paste0("integrated_analysis.", jdate, ".UseTSSfromH3K4me3.likeBM"))
  
  infrdata <- paste0(jprefix, ".RData")
  assertthat::assert_that(file.exists(infrdata))
  
  load(infrdata, v=T)
  
  
  
  # Plot umap  --------------------------------------------------------------
  
  head(out.objs$H3K4me1$annot.glmpca.filt)
  head(out.objs$H3K4me1$annot.louv)
  
  ggplot(out.objs$H3K4me1$annot.louv, aes(x = umap1, y = umap2, color = cluster.y)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point() 
  
  # Get top bins by celltype specific topics -----------------------------------------------------------
  
  
  
  # Get top bins for each celltype ------------------------------------------
  
  tm.result.lst <- lapply(out.objs, function(jout){
    tm.result <- topicmodels::posterior(jout$out.lda)
    tm.result <- AddTopicToTmResult(tm.result, jsep = "")
  })
  
  bins.lst <- lapply(ctypes, function(ctype){
    jtops <- jtops.lst[[ctype]]
    jbins.weights <- tm.result.lst[[ref.mark]]$terms[jtops, ]
    
    if (length(jtops) > 1){
      jbins.weights <- apply(jbins.weights, MARGIN = 2, FUN = max)
    }
    jbins.weights <- sort(jbins.weights, decreasing = TRUE)
    jbins <- names(jbins.weights)
    jbins.filt <- jbins[1:topnbins]
    return(jbins.filt)
  })
  
  print("Bin length")
  print(lapply(bins.lst, length))
  
  # bins 2 genes
  genes.lst <- lapply(bins.lst, function(jbins){
    jgenes <- sapply(jbins, AssignHash, r2g, null.fill = NA)
    # remove NAs
    jgenes <- jgenes[!is.na(jgenes)]
    # jgenes <- gsub("clul1", "yes1", jgenes)
  })
  
  print("Gene length, no NAs")
  print(lapply(genes.lst, length))
  
  # genes 2 TSS 
  
  tss.lst <- lapply(genes.lst, function(jgenes){
    tss.vec <- sapply(jgenes, AssignHash, g2tss, null.fill = NA)
    return(tss.vec)
  })
  
  jcoords.lst <- lapply(tss.lst, function(tss.vec){
    sapply(tss.vec, function(x) strsplit(x, split = ";")[[1]][[1]])
  })
  
  tx.lst <- lapply(tss.lst, function(tss.vec){
    sapply(tss.vec, function(x) strsplit(x, split = ";")[[1]][[2]])
  })
  
  # Create bed file  --------------------------------------------------------
  
  beds.lst <- lapply(ctypes, function(ctype){
    bed.tmp <- data.frame(chromo = sapply(jcoords.lst[[ctype]], JFuncs::GetChromo), 
                          Start = sapply(jcoords.lst[[ctype]], JFuncs::GetStart, returnAsInt = TRUE),
                          End = sapply(jcoords.lst[[ctype]], JFuncs::GetEnd, returnAsInt = TRUE),
                          gene = tx.lst[[ctype]],
                          stringsAsFactors = FALSE) %>%
      rowwise() %>%
      mutate(Midpt = as.integer((Start + End) / 2),
             Start2 = Midpt - 1,
             End2 = Midpt + 1)
    bed.tmp <- subset(bed.tmp, select = c(chromo, Start2, End2, gene)) %>%
      arrange(Start2, End2)
    bed.tmp <- bed.tmp[gtools::mixedorder(bed.tmp$chromo), ]
    return(bed.tmp)
  })
  
  for (ctype in ctypes){
    outf <- file.path(outdir, paste0("ZebrafishWKM_TSS_FromTopics.", ctype, ".bsize_2.bed"))
    bed.out <- beds.lst[[ctype]]
    fwrite(bed.out, file = outf, quote = FALSE, sep = "\t", col.names = FALSE, na = "NA")
  }
}




