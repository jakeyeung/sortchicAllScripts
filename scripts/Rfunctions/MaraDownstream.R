# Jake Yeung
# Date of Creation: 2019-02-10
# File: ~/projects/scchic/scripts/Rfunctions/MaraDownstream.R
# MARA downstream functions

LoadMARA <- function(mdir, bc.inf = "data/barcode_summaries/barcodes/maya_384NLA.bc"){
  act.mat <- fread(file.path(mdir, "Activities"), header = FALSE)
  se.mat <- fread(file.path(mdir, "StandardError"), header = FALSE)
  cnames <- unlist(fread(file.path(mdir, "Colnames"), header = FALSE), use.names = FALSE)
  zscores <- fread(file.path(mdir, "Zscores"), header = FALSE)
  colnames(zscores) <- c("motif", "zscore")
  zscores <- zscores %>% arrange(desc(zscore))
  
  barcodes <- read.table(bc.inf, col.names = "Barcodes", stringsAsFactors = FALSE)
  cellhash <- hash(rownames(barcodes), unlist(barcodes))
  cellhash.bc <- hash(unlist(barcodes), paste("cell", rownames(barcodes), sep = ""))
  
  cnames.new <- sapply(cnames, function(x){
    jmark <- strsplit(x, "\\.")[[1]][[1]]
    jtiss <- strsplit(x, "\\.")[[1]][[2]]
    jmouse <- strsplit(x, "\\.")[[1]][[3]]
    jrep <- paste0("rep", strsplit(x, "\\.")[[1]][[4]])
    jcell <- cellhash.bc[[strsplit(x, "\\.")[[1]][[6]]]]
    xnew <- paste(jtiss, jmark, jmouse, jrep, jcell, sep = "_")
  })
  
  colnames(act.mat) <- c("motif", cnames.new)
  colnames(se.mat) <- c("motif", cnames.new)
  
  act.long <- tidyr::gather(act.mat, -motif, key = "cell", value = "activity")
  
  return(list(act.mat = act.mat, se.mat = se.mat, zscores = zscores, act.long = act.long))
}


LoadLDA <- function(inf){
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  out.lda <- ChooseBestLDA(out.lda)
  (kchoose <- out.lda@k)
  return(out.lda)
}

DoUmapFromLDA <- function(out.lda, custom.settings){
  tm.result <- topicmodels::posterior(out.lda)
  topics.mat <- tm.result$topics
  # terms.mat <- tm.result$terms
  
  # custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)
  # custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)
  
  dat.umap <- umap(topics.mat, config = custom.settings)
  rownames(dat.umap$layout) <- rownames(topics.mat)
  dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout))
  return(dat.umap.long)
}

PlotMotifInUmap <- function(jmotif, dat.merged, zscores, jmark, jsize = 1){
  dat.sub <- subset(dat.merged, motif == jmotif)
  jzscore <- signif(subset(zscores, motif == jmotif)$zscore, digits = 2)
  jtitle <- paste(jmotif, "zscore:", jzscore, jmark)
  m <- ggplot(dat.sub, aes(x = umap1, y = umap2, color = activity)) + geom_point(size = jsize) + 
    ggtitle(jtitle) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_gradient2(low = muted("blue"), mid = "grey95", high = muted("red"), midpoint = mean(dat.sub$activity))
  return(m)
}

