PlotGTrack <- function(x.long, jstart, jend, mart.obj, gen = "mm10", chr = "chr7", jheight = "auto"){
  datheight <- 10
  trackheight <- 2
  itrack <- IdeogramTrack(genome=gen, chromosome=chr)
  gtrack <- BiomartGeneRegionTrack(genome=gen, chromosome=chr, start = jstart, end = jend, collapseTranscript="meta",
                                   # symbol = "Sox6", 
                                   stacking = "squish",
                                   biomart = mart.obj,
                                   name="Genes",showId=T, max.height = 3)
  
  x.sum <- x.long %>%
    group_by(louvain, coord) %>%
    # summarise(exprs = log2(mean(exprs * 1000000 + 1))) %>%
    summarise(exprs = mean(exprs)) %>%
    rowwise() %>%
    mutate(start = as.numeric(GetStart(coord)),
           end = as.numeric(GetEnd(coord)),
           seqnames = GetChromo(coord))
  
  dat <- makeGRangesFromDataFrame(x.sum, keep.extra.columns = TRUE)
  
  # add one track first, then another
  # louvains <- sort(unique(hash::values(clstr)))
  louvains <- sort(unique(x.sum$louvain))
  
  ymin <- min(x.sum$exprs)
  ymax <- max(x.sum$exprs)
  
  dlst <- list()
  dlst[["itrack"]] <- itrack
  dlst[["gtrack"]] <- gtrack
  for (l in louvains){
    print(paste("Iterating cluster", l))
    dlst[[as.character(l)]] <- DataTrack(range=subset(dat, louvain == l, select = c(exprs)), ylim = c(ymin,ymax),
                                         start = jstart, 
                                         end = jend, 
                                         chromo = chr,
                                         genome = gen, name = paste("Clstr", l),
                                         type="mountain")
  }
  if (jheight != "auto"){
    jsizes <- rep(jheight, length(dlst))
    jsizes[1:2] <- c(0.25, 1)
  } else {
    jsizes <- NULL
  }
  plotTracks(dlst, from = jstart, to = jend, collapseTranscripts="meta", sizes = jsizes)
  # plotTracks(dlst, collapseTranscripts="meta", from = jstart, to = jend)
}

PlotImputedPeaks <- function(tm.result, peaks.keep, jchip, show.plot=TRUE, return.plot.only=FALSE, use.count.mat=NULL, usettings=NULL, gname = ""){
  # regions.sub <- subset(regions.annotated, abs(distanceToTSS) < jdist & grepl(jgene, SYMBOL)) %>%
  #   mutate(coord = paste(paste(seqnames, start, sep = ":"), end, sep = "-"),
  #          coord.kb = paste(paste(seqnames, start / 1000, sep = ":"), end / 1000, sep = "-"))
  # if (nrow(regions.sub) == 0){
  #   warning(paste("empty regions.sub for jgene", jgene, "jchip", jchip))
  #   m <- ggplot()
  #   return(m)
  # }
  # peaks.keep <- regions.sub$coord
  # peaks.keep.kb <- regions.sub$coord.kb
  peaks.keep.kb <- peaks.keep
  
  if (is.null(use.count.mat)){
    mat.norm <- t(tm.result$topics %*% tm.result$terms)
  } else {
    mat.norm <- use.count.mat
  }
  
  jlab <- peaks.keep
  
  # run on imputed matrix
  row.i <- which(rownames(mat.norm) %in% peaks.keep)
  
  if (length(row.i) > 1){
    # print(paste("Merging", length(row.i)), "rows")
    jcounts.norm <- Matrix::colSums(mat.norm[row.i, ])
  } else {
    jcounts.norm <- mat.norm[row.i, ]
  }
  # print(range(jcounts.norm))
  jcol.counts <- ColorsByCounts(jcounts.norm, nbreaks = 100, colvec = c("lightblue", "darkblue"))
  # plot topics soft clustering weights
  topics.mat <- tm.result$topics
  
  if (is.null(usettings)){
    nn=5
    jmetric='euclidean' 
    jmindist=0.1
    custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist)
  } else {
    custom.settings <- usettings
  }

  dat.umap <- umap(topics.mat, config = custom.settings)
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  jpeaks.str <- paste(peaks.keep.kb, collapse = ",")
  jmain <- paste0(jchip, " ", jlab, ' npeaks ', length(row.i), "\n", gname)
  # prepare plot object
  dat <- data.frame(umap1 = dat.umap$layout[, 1], 
                    umap2 = dat.umap$layout[, 2], 
                    jcol = jcol.counts,
                    counts.norm = jcounts.norm)
  m <- ggplot(dat, aes(x = umap1, y = umap2, col = jcol.counts)) + geom_point() + scale_color_identity() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=7), legend.position = "bottom") + 
    ggtitle(jmain)
  
  if (show.plot){
    print(m)
  }
  
  # plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, 
  #      main = jmain, col = jcol.counts, asp = 0.75)
  if (!return.plot.only){
    return(list('m' = m))
  } else {
    return(m)
  }
}
