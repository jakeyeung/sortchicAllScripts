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

PlotMotifInUmap <- function(jmotif, dat.merged, zscores, jmark, jsize = 1, colvec = c("blue", "gray95", "red")){
  dat.sub <- subset(dat.merged, motif == jmotif)
  jzscore <- signif(subset(zscores, motif == jmotif)$zscore, digits = 2)
  jtitle <- paste(jmotif, "zscore:", jzscore, jmark)
  # add ordering statistic
  # https://stackoverflow.com/questions/11805354/geom-point-put-overlapping-points-with-highest-values-on-top-of-others
  dat.sub <- dat.sub %>%
    group_by(motif) %>%
    mutate(orderrank = rank(activity, ties.method = "first")) %>%
    arrange(orderrank)
  m <- ggplot(dat.sub, aes(x = umap1, y = umap2, color = activity, order = orderrank)) + geom_point(size = jsize) +
  # m <- ggplot(dat.sub, aes(x = umap1, y = umap2, color = activity)) + geom_point(size = jsize) + 
    ggtitle(jtitle) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  if (length(colvec) == 3){
    midpt <- min(dat.sub$activity) + (max(dat.sub$activity) - min(dat.sub$activity)) / 2
    m <- m + 
      scale_color_gradient2(low = colvec[[1]], mid = colvec[[2]], high = muted(colvec[[3]]), midpoint = midpt)
  } else if (length(colvec) == 2){
    m <- m + 
      scale_color_gradient(low = colvec[[1]], high = muted(colvec[[2]]), space = "Lab")
  } else {
    warning(paste("Colvec must be length 2 or 3"))
    print(colvec)
  }
  return(m)
}

