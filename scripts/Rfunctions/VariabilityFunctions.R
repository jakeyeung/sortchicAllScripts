
GetDiffRelToCell <- function(imputed.dat, gstr, trajs, trajname, dat.umap.long){
  hsc.cell <- (trajs[[jmark]][[trajname]] %>% arrange(lambda) %>% dplyr::top_n(-1))$cell[[1]]
  jsub.all <- MatToLong(imputed.dat, gstr = gstr, cells.vec=NULL)
  jsub.hsc <- jsub.all %>% filter(cell == hsc.cell)
  # plot by reference to stem cell 
  jsub.hsc.ref <- jsub.hsc %>% rename(exprs.ref = exprs) %>% select(-cell, -start, -end, -pos, -chromo)
  jsub.ref <- left_join(jsub.all, jsub.hsc.ref)
  # do the difference over pseudotime?? 
  jsub.ref$exprs.diff <- log2(jsub.ref$exprs) - log2(jsub.ref$exprs.ref)
  jsub.ref.sum <- jsub.ref %>%
    group_by(cell) %>%
    summarise(exprs.diff.med = median(abs(exprs.diff)))
  # join to UMAP 
  jsub.ref.merge <- left_join(jsub.ref.sum %>% dplyr::select(cell, exprs.diff.med), dat.umap.long) %>%
    mutate(label = gstr)
  return(jsub.ref.merge)
}

GetCellSd <- function(dat.mat, grep.str, log2.scale = TRUE, fn = sd){
  # calculate standard deviation from matrix
  row.filt.indx <- grepl(grep.str, rownames(dat.mat))
  if (log2.scale){
    cell.sd.df <- data.frame(cell = colnames(dat.mat[row.filt.indx, ]),
                             cell.sd = apply(log2(dat.mat[row.filt.indx, ]), 2, fn))
  } else {
    cell.sd.df <- data.frame(cell = colnames(dat.mat[row.filt.indx, ]),
                             cell.sd = apply(dat.mat[row.filt.indx, ], 2, fn))
  }
  cell.sd.df$label <- grep.str
  return(cell.sd.df)
}

MatToLong <- function(imputed.dat, gstr, cells.vec = NULL){
  if (!is.null(cells.vec)){
    jsub <- as.data.frame(imputed.dat[grepl(gstr, rownames(imputed.dat)), cells.vec])
    colnames(jsub) <- cells.vec
  } else {
    jsub <- as.data.frame(imputed.dat[grepl(gstr, rownames(imputed.dat)), ])
  }
  if (nrow(jsub) == 0){
    print(paste("Warning: grepstr", gstr, "found no matches, returning empty dataframe"))
    return(data.frame(NULL))
  }
  jsub$coord <- rownames(jsub)
  jsub$start <- as.numeric(sapply(jsub$coord, GetStart))
  jsub$end <- as.numeric(sapply(jsub$coord, GetStart))
  jsub$pos <- jsub$start + (jsub$end - jsub$start) / 2
  jsub$chromo <- sapply(jsub$coord, GetChromo)
  jsub <- gather(jsub, key = "cell", value = "exprs", c(-coord, -start, -end, -pos, -chromo))
  jsub <- jsub %>% arrange(desc(pos))
  return(jsub)
}

MergeSdWithPseudotime <- function(dat.umap.long.trajs, tm.result.lst, trajs, jmark, jtraj, grep.strs, jscale=TRUE, jfn = mad){
  imputed.dat <- t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics)
  dat.umap.long <- dat.umap.long.trajs[[jmark]]
  cell.sd.df.long <- lapply(grep.strs, function(grep.str){
    return(GetCellSd(imputed.dat, grep.str, log2.scale = jscale, fn = jfn))
  }) %>%
    bind_rows()
  dat.umap.filt <- left_join(dat.umap.long, cell.sd.df.long)
  # add a trajectory
  dat.umap.filt <- left_join(trajs[[jmark]][[jtraj]] %>% dplyr::select(cell, lambda), dat.umap.filt, by = "cell")
  return(dat.umap.filt)
}


CalculateACF <- function(jsub.hsc, jstep = 20000, jtype = "correlation", jmain = "Title", show.plot = TRUE, maxlag = "full"){
  # impute missing positions with minimum value
  # impute missing bins with minimum value
  jcells <- unique(jsub.hsc$cell)
  pos.start <- min(jsub.hsc$pos)
  pos.end <- max(jsub.hsc$pos)
  # print(head(jsub.hsc))
  # print(tail(jsub.hsc))
  # print(paste(pos.start, pos.end))
  # jstep <- 20000
  pos.vec <- seq(pos.start, pos.end, jstep)
  jsub.impute.vec <- data.frame(pos = rep(pos.vec, length(jcells)), cell = rep(jcells, each = length(pos.vec)))
  jsub.impute.vec <- left_join(jsub.impute.vec, jsub.hsc %>% dplyr::select(c(chromo, pos, cell, exprs)))
  # jsub.impute.vec$exprs[which(is.na(jsub.impute.vec$exprs))] <- min(jsub.hsc$exprs)
  
  if (maxlag == "full"){
    maxlag <- nrow(jsub.impute.vec)
  }
  
  acf.out <- acf(log2(jsub.impute.vec$exprs), type = jtype, lag.max = nrow(jsub.impute.vec), na.action = na.pass, main = jmain, plot = show.plot)
  acf.out$lag.stepadj <- acf.out$lag * jstep
  return(acf.out)
}


SummarizeACF <- function(acf.out.hsc.lst){
  acf.out.hsc.long <- lapply(names(acf.out.hsc.lst), function(g){
    lst <- acf.out.hsc.lst[[g]]
    dat.tmp <- data.frame(acfval = lst$acf, 
                          dx = lst$lag.stepadj, 
                          label = g)
    return(dat.tmp)
  }) %>% bind_rows() 
  
  acf.out.hsc.sum <- acf.out.hsc.long %>%
    filter(dx < pos.max) %>%
    group_by(dx) %>%
    summarise(acfval = mean(acfval))
  
  # get error bar for each chromose
  acf.hsc.ci <- acf.out.hsc.long %>%
    group_by(label) %>%
    summarise(ci = ci.interval / sqrt(length(acfval))) %>%
    ungroup() %>%
    summarise(ci = mean(ci))
  acf.hsc.ci <- acf.hsc.ci$ci[[1]]
  return(list(acf.out.sum = acf.out.hsc.sum, acf.ci = acf.hsc.ci))
}

