# Jake Yeung
# Date of Creation: 2019-01-09
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/summarize_histone_marks_compare_metacell.R
# Summarize histone marks, compare with metacell
# downloading /hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisHiddenDomains_1000/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.TRUE to tmp dir

rm(list=ls())

library(topicmodels)
library(dplyr)
library(ggplot2)
library(umap)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rGREAT)
library(hash)

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")

GetBins <- function(jchromo, midpt, jname, Q, winsize = 100000L){
  jleft <- midpt - winsize / 2
  jright <- midpt + winsize / 2
  Qsub <- subset(Q, chromo == jchromo & coord > jleft & coord < jright)
  Qsub$name <- jname
  return(as.data.frame(Qsub))
}

Vectorize(AssignHash <- function(x, jhash, null.fill = NA){
  # assign hash key to hash value, handle NULLs
  # null.fill = "original", returns original value x into jhash
  x.mapped <- jhash[[as.character(x)]]
  if (is.null(x.mapped)){
    if (is.na(null.fill)){
      x.mapped <- null.fill
    } else if (as.character(null.fill) == "original"){
      x.mapped <- x
    } else {
      x.mapped <- null.fill
    }
  }
  return(x.mapped)
}, vectorize.args = "x")

GetGOData <- function(out.tb.lst, ontology, jdecreasing=TRUE, order.by="Hyper_Fold_Enrichment"){
  topics <- which(lapply(1:length(out.tb.lst), function(i) !is.null(out.tb.lst[[i]][[ontology]])) == TRUE)
  GOdata <- lapply(topics, function(i) out.tb.lst[[i]][[ontology]])
  GOdata <- lapply(1:length(GOdata), function(i) GOdata[[i]][order(GOdata[[i]][order.by], 
                                                                   decreasing = jdecreasing),])
  GOdata <- lapply(1:length(GOdata), function(i) GOdata[[i]][1:top,])
  
  # annotate each GOdata by topic number
  for (i in topics){
    GOdata[[i]]$topic <- i
  }
  GOdata <- dplyr::bind_rows(GOdata) %>%
    mutate(topic = factor(topic, levels = sort(unique(topic))),
           onto = ontology)
  return(GOdata)
}



# Parameters to get files -------------------------------------------------




outdir <- "~/Dropbox/scCHiC_figs/FIG4_BM/LDA_outputs/test"
dir.create(outdir)


jchip <- "H3K9me3"
jchip <- "H3K4me3"
jchip <- "H3K27me3"
jchip <- "H3K4me1"

# do all
jchips <- c("H3K27me3", "H3K4me1", "H3K9me3", "H3K4me3")
# jchips <- c("H3K4me3", "H3K4me1")
  for (jchip in jchips){
  
  jdist <- 1000L
  jmean <- 1
  jmin <- 100L
  jmax <- 500000L
  # binarize <- "TRUE";  jtops <- "5_7_10_12_15_20_25_30"
  binarize <- "FALSE"; jtops <- "15_20_25_30_35"
  
  jdir <- paste0('/tmp/ldaAnalysisHiddenDomains_', jdist, '/lda_outputs.meanfilt_', jmean, '.cellmin_', jmin, '.cellmax_', jmax, '.binarize.', binarize)
  inf <- file.path(jdir, paste0('lda_out_meanfilt.PZ-BM-', jchip, '.CountThres0.K-', jtops, '.Robj'))
  infbase <- basename(inf)
  infbase <- strsplit(infbase, ".Robj")[[1]][[1]]
  
  # 0.98 threshold 
  # inf.GREAT <- file.path(jdir, "downstream", paste0(infbase, ".GREAT.Robj"))
  # 0.96 threshold 
  inf.GREAT <- file.path(jdir, "downstream", paste0(infbase, ".GREAT.0.96.Robj"))
  
  inf.mc <- file.path(paste0("~/Dropbox/scCHiC_figs/FIG4_BM/", jchip, ".datadir_mc_f.Rda"))
  
  assertthat::assert_that(file.exists(inf))
  assertthat::assert_that(file.exists(inf.GREAT))
  assertthat::assert_that(file.exists(inf.mc))
  
  load(inf, v=T)
  load(inf.GREAT, v=T)
  
  lapply(topic.regions, length)
  
  
  # what's the sparsity?
  print(1 - Matrix::nnzero(count.mat) / length(count.mat))
  
  regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                  TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                  annoDb='org.Mm.eg.db'))
  regions.annotated$region_coord <- rownames(regions.range)
  
  # if already loaded the GREAT object, then don't need to choose best K
  
  if (length(out.lda) > 1){
    out.lda.lst <- out.lda
    # we did multicore, so which one to choose?
    Kvec <- sapply(out.lda.lst, function(x) x@k)
    best.K <- Kvec[which.max(sapply(out.lda.lst, function(x) x@loglikelihood))]
    # plot loglikelihood
    par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    plot(Kvec, sapply(out.lda.lst, function(x) x@loglikelihood), 'o')
    kchoose <- best.K
    out.lda <- out.lda.lst[[which(Kvec == kchoose)]]
    print(paste("Likelihood: ", out.lda@loglikelihood))
  } else {
    kchoose <- out.lda@k
  }
  
  tm.result <- posterior(out.lda)
  
  # print(head(tm.result$topics))
  # cluster the terms, but first get interesting terms
  top.regions <- unique(unlist(topic.regions))  # topic.regions defined by threshold, 98 or 96th percentile of top weights in each column of the betas matrix
  terms.mat <- t(tm.result$terms)[top.regions, ]
  
  # Compare with metacell ---------------------------------------------------
  
  mc.out <- GetMCObjects(inf = inf.mc)
  mc_index <- mc.out$mc_index; mc_colors <- mc.out$mc_colors
  barcodes <- read.table("data/barcode_summaries/barcodes/maya_384NLA.bc", col.names = "Barcodes", stringsAsFactors = FALSE)
  cellhash <- hash(rownames(barcodes), unlist(barcodes))
  cellhash.bc <- hash(unlist(barcodes), paste("cell", rownames(barcodes), sep = ""))
  
  cellnames <- unname(out.lda@documents)
  experinames <- unique(sapply(cellnames, function(x) paste(strsplit(x, "-")[[1]][-6], collapse = "-"), USE.NAMES = FALSE))
  # rename technical replicates 
  replhash <- hash()
  mice <- c("m1", "m2")
  for (m in mice){
    replhash <- GetReplHash(experinames, m, replhash)
  }
  # add cell barcode
  cellnames.new <- sapply(cellnames, MakeNewCellName, replhash, cellhash.bc, add.m=FALSE, USE.NAMES = FALSE)
  cellnames.hash <- hash(cellnames.new, cellnames)
  
  # create new mc_index
  mc_index.newnames <- mc_index
  names(mc_index.newnames) <- sapply(names(mc_index), function(x) cellnames.hash[[x]])
  mc_index.newnames.filt <- mc_index.newnames[which(names(mc_index.newnames) != "NULL")]
  
  cells.clstr.hash <- hash(names(mc_index.newnames.filt), mc_index.newnames.filt)
  cells.color.hash <- hash(as.character(seq(mc_colors)), mc_colors)
  
  # filter for common cells between the two analyses
  topics.mat <- tm.result$topics
  cells.lda <- unname(rownames(topics.mat))
  cells.metacell <- names(mc_index.newnames.filt)
  
  cells.common <- intersect(cells.lda, cells.metacell)  # 997 cells left
  
  topics.mat <- topics.mat[cells.common, ]
  cells.indx <- sapply(rownames(topics.mat), function(x) cells.clstr.hash[[x]])
  cells.rgb <- sapply(cells.indx, function(x) cells.color.hash[[as.character(x)]])
  
  
  # Check Hox cluster genes -------------------------------------------------
  
  # use same region as AvO, but can also filter with regions.annotated
  
  #HoxA chr6 midpoint 52.2 Mb
  Q <- mc.out$Q
  chromos <- c("chr6", "chr11", "chr15", "chr2")
  midpoints <- c(52200000L, 96300000L, 103000000L, 74700000L)
  jnames <- c("HoxA", "HoxB", "HoxC", "HoxD")
  
  # chromos <- c("chr7")
  # midpoints <- c(103830000L)
  # jnames <- c("Hbb")
  
  i <- 1
  Qsub <- subset(Q, chromo == chromos[[i]] & coord > midpoints[[i]]-50000 & coord < midpoints[[i]]+50000)
  Qchr.lst <- GetBins(chromos[[i]], midpoints[[i]], jnames[[i]], Q, 100000L)
  Qchr.lst <- bind_rows(mapply(GetBins, chromos, midpoints, jnames, MoreArgs = list(Q = Q, winsize = 100000), SIMPLIFY = FALSE))
  
  # create GR
  Qgr <- GRanges(seqnames = Qchr.lst$chromo, 
                 ranges = IRanges(start = Qchr.lst$left, end = Qchr.lst$right), 
                 name = Qchr.lst$name)
  genes.keep <- c("HoxA", "HoxB", "HoxC", "HoxD")
  genes.keep.str <- paste(jnames, collapse = ",")
  Qgr <- subset(Qgr, name %in% genes.keep)
  
  # Assign to peaks in the LDA ----------------------------------------------
  regions.intersect <- subsetByOverlaps(regions.range, Qgr)
  peaks.keep <- names(sort(terms.mat[, 5], decreasing = TRUE)[1:50])
  row.i <- which(rownames(count.mat) %in% peaks.keep)
  jcounts.total <- Matrix::colSums(count.mat[, cells.common])
  jcounts <- Matrix::colSums(count.mat[row.i, cells.common])
  jcounts.norm <- jcounts / jcounts.total
  # jcounts.norm <- jcounts.total
  jcol.counts <- ColorsByCounts(jcounts.norm, nbreaks = 100, colvec = c("lightblue", "blue", "darkblue"))
  hoxcounts.hash <- hash(jcounts.norm, names(jcounts))
  
  
  # Process GREAT output ----------------------------------------------------
  
  (jontos <- unique(unlist(lapply(out.tb.lst, function(x) names(x)))))
  
  # no need cutoff, GREAT output already filtered
  # cutoff <- 0.05
  # foldchange <- 1.5
  
  # ontology <- jontos[[8]]
  # ontology <- jontos[[14]]
  # ontology <- jontos[[8]]  # do them all
  order.by <- "Hyper_Fold_Enrichment"
  jdecreasing <- TRUE
  top <- 3
  # 
  
  GOdata.all <- bind_rows(lapply(jontos, function(onto) GetGOData(out.tb.lst, onto, jdecreasing, order.by)))
  
  
  
  # Plotting everything now ----------------------------------------------------------------
  

  # # UMAP settings 
  # nn <- 5
  # jmetric <- 'euclidean'
  # # jmetric <- 'cosine'
  # jmindist <- 0.1
  # custom.settings <- umap.defaults
  # custom.settings$n_neighbors <- nn
  # custom.settings$metric <- jmetric
  # custom.settings$min_dist <- jmindist

  nn=5
  nnterms=15
  jmetric='euclidean' 
  jmindist=0.1
  custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist)
  custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)
  # # use different settings for peak weights for H3K4me1 and H3K4me3 for some reason
  # if (jchip %in% c("H3K4me1", "H3K4me3")){
  #   custom.settings.terms <- GetUmapSettings(nn=15, jmetric='euclidean', jmindist=0.01)
  # } else {
  #   custom.settings.terms <- custom.settings
  # }
  

  dat.umap <- umap(topics.mat, config = custom.settings)
  rownames(dat.umap$layout) <- rownames(topics.mat)
  jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)
  
  # # run umap again on all cells or just common between MetaCell and LDA?
  # dat.umap.all <- umap(tm.result$topics, config = custom.settings)
  # rownames(dat.umap.all$layout) <- rownames(tm.result$topics)
  pdf(file.path(outdir, paste0("lda_output", jchip, ".K.", kchoose, ".pdf")), useDingbats = FALSE, paper = "A4r")
  
  # compare metacell and lda
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = "Color by MetaCell", col = cells.rgb, asp = 1, xlab = "UMAP Dim 1", ylab = "UMAP Dim 2")
  
  # plot topics soft clustering weights
  jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "blue", "darkblue"))
  nb.col <- 5
  nb.row <- ceiling(kchoose / nb.col)
  par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
  mapply(function(jcol.rgb, jtopic){
    plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
  }, jcol.rgbs, seq(kchoose))
  
  # plot topics, plot by Hox counts (size normalized)
  # Color by normalized UMI counts
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, 
       main = paste("Normalized", genes.keep.str, "counts"), col = jcol.counts, asp = 0.75)
  
  # count median counts across metacells
  out.dat <- data.frame(X1 = dat.umap$layout[, 1], 
                        X2 = dat.umap$layout[, 2], 
                        cellname = rownames(dat.umap$layout), stringsAsFactors = FALSE)
  out.dat$mc <- sapply(out.dat$cellname, AssignHash, cells.clstr.hash)
  out.dat$hoxcounts <- jcounts.norm
  # Total sum needs to be normalized? What does metacell say here?
  out.dat %>%
    group_by(mc) %>%
    summarize(hoxcounts.median = median(hoxcounts)) %>%
    ggplot() + geom_bar(aes(x = as.character(mc), y = hoxcounts.median), stat = "identity") + 
    ggtitle(paste("Median counts in MetaCell for ", genes.keep.str))
  
  # plot peaks soft clustering weights. Should be about K clusters visually. Overlaps are probably peaks that are important in multiple clusters

  dat.umap.terms <- umap(terms.mat, config = custom.settings.terms)
  # downsample rows for plotting purposes
  downsamp.i <- sample(seq(nrow(dat.umap.terms$layout)), size = round(0.1 * nrow(dat.umap.terms$layout)), replace = FALSE)
  jcol.rgb.terms <- lapply(seq(kchoose), ColorsByGamma, terms.mat[downsamp.i, ], c("pink", "red", "darkred"))
  
  par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
  mapply(function(jcol.rgb.term, jtopic){
    plot(dat.umap.terms$layout[downsamp.i, 1], dat.umap.terms$layout[downsamp.i, 2], 
         col = jcol.rgb.term, pch = 20, asp = 0.75,
         main = paste("Peak Weights, T", jtopic))
  }, jcol.rgb.terms, seq(kchoose))
  
  # Show GREAT gene enrichment
  
  ontos <- c("MSigDB Immunologic Signatures", "MSigDB Predicted Promoter Motifs", "GO Biological Process")
  # pdf("/tmp/test.pdf", paper = "A4r")
  for (jonto in ontos){
    m <- ggplot(subset(GOdata.all, onto == jonto), 
                aes(x = topic, y = stringr::str_wrap(name, 100), size = Hyper_Fold_Enrichment, color = -log10(Hyper_Adjp_BH))) +
      geom_point() + theme_bw() +
      theme(aspect.ratio=1, panel.grid.minor = element_blank(), 
            axis.text.y =element_text(size=5, angle=0), 
            axis.text.x = element_text(size=7, angle=0),
            legend.position = "bottom",
            legend.text = element_text(size=7),
            legend.title = element_text(size=7)) +
      xlab("") + ylab("")
    print(m)
  }
  dev.off()
  
# Print data tables for downstream analysis
  
# tm.result is the output of matrix decomposition by LDA. 
# Each cell is a distribution over topics (tm.result$topics), each topic is a distribution over peaks (tm.result$terms). 
# topic.regions is top 96%ile regions for each topic. Basically cutoff to assign peak to topic
# regions.annotated: annotations for each region 
# GOdata.all: Output of GREAT analysis for each topic
save(tm.result, topic.regions, regions.annotated, GOdata.all, file = file.path(outdir, paste0("output_data.", jchip, ".Rda")))
}



# Cool, can you see this signature in the betas? --------------------------

# (x <- terms(out.lda, 1000))
# 
# # can we get the Hox regions in one of the topic regions?
# hox.filt <- lapply(topic.regions, function(x) x[which(x %in% peaks.keep)])
# 
# peaks.label <- sapply(rownames(terms.mat), function(x) ifelse(x %in% peaks.keep, x, NA))
# # where do hox genes weigh in topic 3 and 5?
# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# plot(log10(terms.mat[, 5]), log10(terms.mat[, 3]), pch = 20)
# text(log10(terms.mat[, 5]), log10(terms.mat[, 3]), labels = peaks.label)
# 
# mat.sub <- terms.mat[, c(3, 5)]
# mat.sub <- mat.sub[order(rowSums(mat.sub^2), decreasing = TRUE), ]

# GREAT -------------------------------------------------------------------

# load(inf.GREAT, v=T)





