FitGetPval <- function(jsub, jform){
  jfit <- lm(formula = jform, data = jsub)
  pval <- summary(jfit)$coefficients[2, "Pr(>|t|)"]
  slope.val <- summary(jfit)$coefficients[2, "Estimate"]
  slope.se <- summary(jfit)$coefficients[2, "Std. Error"]
  return(data.frame(pval = pval, slope.val = slope.val, slope.se = slope.se))
}

SumSqrDev <- function(x){
  return( sum( (x - mean(x)) ^ 2 ))
}

BinTrajectory <- function(trajs.spring.lst, jtraj, nearest = 0.1){
  round.int <- 1 / nearest
  trajs.sum <- lapply(trajs.spring, function(x) x[[jtraj]] %>% mutate(traj = jtraj)) %>%
    bind_rows() %>%
    rowwise() %>%
    mutate(mark = strsplit(cell, "_")[[1]][[2]]) %>%
    left_join(., mat.sub.merge) %>%
    rowwise() %>%
    mutate(lambda.bin = floor(lambda * round.int) / round.int) %>%
    group_by(traj, lambda.bin, mark, coord, pos) %>%
    summarise(exprs = mean(exprs)) %>%
    return(trajs.sum)
}


GetTrajSum <- function(tm.result.lst, trajs.mixed, jmarks, jstr, jpseudo, jfac, jtraj){
  mat.sub.merge <- lapply(jmarks, function(jmark) GetMatSub(tm.result.lst, jmark, jstr, jpseudo, jfac) %>% mutate(mark = jmark)) %>% 
    bind_rows()
  trajs.long <- lapply(trajs.mixed, function(x) x[[jtraj]]) %>%
    bind_rows() %>%
    rowwise() %>%
    mutate(mark = strsplit(cell, "_")[[1]][[2]]) %>%
    left_join(., mat.sub.merge) %>%
    rowwise() %>%
    mutate(lambda.bin = floor(lambda * 10) / 10)
  trajs.sum <- trajs.long %>%
    group_by(lambda.bin, mark, coord, pos) %>%
    summarise(exprs = mean(exprs)) %>%
    mutate(chromo = jstr)
  return(trajs.sum)
}

FitSlope <- function(dat.sub){
  # fit exprs to trajectory to find slope
  fit <- lm(formula = exprs ~ lambda.bin, data = dat.sub)
  slope <- fit$coefficients[['lambda.bin']]
  int <- fit$coefficients[['(Intercept)']]
  pval <- summary(fit)$coefficients["lambda.bin", "Pr(>|t|)"]
  return(data.frame(slope = slope, int = int, pval = pval))
}

HexPlot <- function(data, mapping, ...){
  p <- ggplot(data, mapping) + geom_hex() + geom_density2d() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(p)
}

MakeLowerPlot <- function(data, mapping, ..., jbins = 150){
  m <- ggplot(data = data, mapping = mapping) + 
    geom_point_rast(data = data %>% sample_frac(size = 0.1), alpha = 0.3) + 
    geom_density_2d() 
  # geom_point(data = data %>% sample_frac(size = 0.1), alpha = 0.1) + 
  return(m)
}




GetCellsAlongTraj <- function(trajs.spring, jmark, ctype, n.sections = 3, thres = 0.05, jfrom = 0, jto = 1, show.plots = TRUE){
  lambda.cutoffs <- seq(from = jfrom, to = jto, length.out = n.sections)
  cells.vec <- lapply(lambda.cutoffs, function(jlambda) (trajs.spring[[jmark]][[ctype]] %>% arrange(abs(lambda - jlambda)))$cell[[1]])
  if (show.plots){
    m <- ggplot(dat.trajs.long %>% filter(mark == jmark) %>% mutate(highlighted = cell %in% c(cell0, cell1, cell2)), aes(x = X1, y = X2, color = highlighted)) +
      geom_point() +
      theme_bw() +
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
    # plot chromosome 15 over time
    m <- ggplot(jmerge %>% filter(cell %in% c(cell0, cell1, cell2)), aes(x = pos, y = exprs, group = lambda, color = lambda)) + 
      geom_line() + theme_bw() + facet_wrap(~lambda, ncol = 1) + 
      theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
      ggtitle(jmark, ctype) 
    print(m)
  }
  return(unlist(cells.vec))
}


GetMatSub <- function(tm.result.lst, jmark, gstr, jpseudo, jfac, cells.vec = NULL){
  imputed.dat <- log2((t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics) * jfac) + jpseudo)
  mat.sub <- MatToLong(imputed.dat, gstr = gstr, cells.vec = cells.vec) %>% dplyr::select(-start, -end)
  # mat.sub <- MatToLong(imputed.dat, gstr = gstr, cells.vec = NULL)
  return(mat.sub)
}


MakeVariabilityPlots <- function(jmark, trajname, tm.result.lst, dat.umap.long.trajs, 
                                 jcol = c("gray85", "gray50", "darkblue"), grep.strs = paste("chr", c(seq(21)), ":", sep = ""), jalpha = 0.5, pseudo = 0, jscale = 1, 
                                 mdpt.sd = 1, ms.sd = c(0, 3), mdpt.fc = 0.75, lims.fc = c(0, 3), 
                                 jsize.facet = 0.2, gw.size.facet = 2, lagmax = 2500, ci.interval = 1.96,
                                 chromogstr="chr15:",
								 pdfout = FALSE){
  
  # jmark <- "H3K27me3"
  # trajname <- "myeloid"

  if (!is.null(pdfout)){
    pdf(pdfout, useDingbats=FALSE)
  }
  
  names(grep.strs) <- grep.strs  # SummarizeACF() uses a named list that fails if you don't do this 
  
  print(paste("Making plots for", jmark, trajname))
  
  imputed.dat <- t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics)
  dat.umap.long <- dat.umap.long.trajs[[jmark]]
  
  cell.sd.df.long <- lapply(grep.strs, function(grep.str){
    return(GetCellSd(jscale * (imputed.dat + pseudo), grep.str, log2.scale = TRUE, fn = sd))
  }) %>%
    bind_rows()
  
  dat.umap.filt <- left_join(dat.umap.long, cell.sd.df.long)
  
  m.chr <- ggplot(dat.umap.filt, aes(x = umap1, y = umap2, color = cell.sd)) + 
    geom_point(size = jsize.facet)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.sd, limits = lims.sd) + 
    ggtitle(jmark, paste0(deparse(substitute(sd)), " across chromosome"))
  print(m.chr)
  
  # what about genome wide
  cell.sd.genomewide <- GetCellSd(imputed.dat, "", log2.scale=TRUE)
  dat.umap.filt.gw <- left_join(dat.umap.long, cell.sd.genomewide)
  
  m.gw <- ggplot(dat.umap.filt.gw, aes(x = umap1, y = umap2, color = cell.sd)) + 
    geom_point(size = gw.jsize.facet)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.sd, limits = lims.sd) + 
    ggtitle(jmark, paste0(deparse(substitute(sd)), " genome wide"))
  
  print(m.gw)
  
  
  # Highlight differences for two representative cells on a representative chromosome 
  hsc.cell <- (trajs[[jmark]][[trajname]] %>% arrange(lambda) %>% dplyr::top_n(-1))$cell[[1]]
  diff.cell <- (trajs[[jmark]][[trajname]] %>% arrange(lambda) %>% dplyr::top_n(1))$cell[[1]]
  
  gstr <- paste0(chromogstr)
  jsub <- MatToLong(imputed.dat, gstr, cells.vec = c(hsc.cell, diff.cell))
  m.spatial <- ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs))) + 
    geom_line(alpha = jalpha) + 
    facet_wrap(~cell) + 
    ggtitle(paste(jmark, gstr)) + 
    xlab("MB") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  m.spatial.merged <- ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs), group = cell, color = cell)) + 
    geom_line(alpha = jalpha) + 
    ggtitle(paste(jmark, gstr)) + 
    xlab("MB") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  m.spatial.log2fc <- ggplot(jsub %>% group_by(pos) %>% summarise(exprs = diff(log2(exprs))), aes(x = pos / 10^6, y = exprs)) + 
    geom_line(alpha = jalpha) + 
    ggtitle(paste(jmark, gstr)) + 
    xlab("MB") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ylab("log2 Fold Change")
  
  print(m.spatial)
  print(m.spatial.merged)
  print(m.spatial.log2fc)
  
  
  # spatial pattern?
  
  jsub.hsc <- jsub %>% filter(cell == hsc.cell) %>% arrange(pos)
  jsub.myeloid <- jsub %>% filter(cell == diff.cell) %>% arrange(pos)
  
  # jmain <- "hsc"
  acf.out.hsc <- CalculateACF(jsub.hsc, jstep = jstep, jtype = "correlation", jmain = paste(jmark, trajname, "Prog Cell", gstr), show.plot = TRUE)
  acf.out.hsc <- CalculateACF(jsub.hsc, jstep = jstep, jtype = "partial", jmain = paste(jmark, trajname, "Prog Cell", gstr), show.plot = TRUE)
  acf.out.myeloid <- CalculateACF(jsub.myeloid, jstep = jstep, jtype = "correlation", jmain = paste(jmark, trajname, "Diff Cell", gstr), show.plot = TRUE)
  acf.out.myeloid <- CalculateACF(jsub.myeloid, jstep = jstep, jtype = "partial", jmain = paste(jmark, trajname, "Diff Cell", gstr), show.plot = TRUE)
  
  # do it genome wide
  jsub.hsc.lst <- lapply(grep.strs, function(g) MatToLong(imputed.dat, g, cells.vec = c(hsc.cell)))
  jsub.myeloid.lst <- lapply(grep.strs, function(g) MatToLong(imputed.dat, g, cells.vec = c(diff.cell)))
  
  # plot chromosome over space for all chromosomes
  
  m.hsc.chromo.all <- ggplot(jsub.hsc.lst %>% bind_rows(), aes(x = pos / 10^6, y = log2(exprs))) + geom_line() + facet_wrap(~chromo, scales = "free_x", ncol = 7) + 
    theme_bw(8) + theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.x = element_text(size = 4)) + 
    xlab("MB") + ylab("log2(imputed counts)") + ggtitle(paste(jmark, trajname, "Prog Cell"))
  print(m.hsc.chromo.all)
  
  m.myeloid.chromo.all <- ggplot(jsub.myeloid.lst %>% bind_rows(), aes(x = pos / 10^6, y = log2(exprs))) + geom_line() + facet_wrap(~chromo, scales = "free_x", ncol = 7) + 
    theme_bw(8) + theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.x = element_text(size = 4)) + 
    xlab("MB") + ylab("log2(imputed counts)") + ggtitle(paste(jmark, trajname, "Diff Cell"))
  print(m.myeloid.chromo.all)
  

  acf.out.hsc.lst <- lapply(jsub.hsc.lst, function(jsub.hsc) CalculateACF(jsub.hsc, jstep = jstep, jtype = "correlation", jmain = paste(jmark, trajname, "Prog Cell", gstr), show.plot = FALSE))
  acf.out.myeloid.lst <- lapply(jsub.myeloid.lst, function(jsub.hsc) CalculateACF(jsub.hsc, jstep = jstep, jtype = "correlation", jmain = paste(jmark, trajname, "Diff Cell", gstr), show.plot = FALSE))
  
  # average out the plots for different lags 
  # plot all chromosomes
  par(mfrow=c(1, 2), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0, pty = "s")
  for (i in seq(length(acf.out.hsc.lst))){
    plot(acf.out.hsc.lst[[i]]$lag.stepadj / 10^6, acf.out.hsc.lst[[i]]$acf, main = paste(jmark, trajname, "Prog", grep.strs[[i]]), type = "h", xlab = "Step Size (MB)", ylab = "Autocorrelation")
    abline(h = ci.interval / sqrt(length(acf.out.hsc.lst[[i]]$lag.stepadj)), lty = "dotted", col = "blue")
    abline(h = -ci.interval / sqrt(length(acf.out.hsc.lst[[i]]$lag.stepadj)), lty = "dotted", col = "blue")
    abline(h = 0, cex = 2)

    plot(acf.out.myeloid.lst[[i]]$lag.stepadj / 10^6, acf.out.myeloid.lst[[i]]$acf, main = paste(jmark, trajname, "Diff", grep.strs[[i]]), type = "h", xlab = "Step Size (MB)", ylab = "Autocorrelation")
    abline(h = ci.interval / sqrt(length(acf.out.myeloid.lst[[i]]$lag.stepadj)), lty = "dotted", col = "blue")
    abline(h = -ci.interval / sqrt(length(acf.out.myeloid.lst[[i]]$lag.stepadj)), lty = "dotted", col = "blue")
    abline(h = 0, cex = 2)
  }
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  
  # summarize across chromosomes
  acf.hsc.sum.lst <- SummarizeACF(acf.out.hsc.lst)
  m.acf.hsc.gw <- ggplot(acf.hsc.sum.lst$acf.out.sum, aes(x = dx / 10^6, y = acfval)) + geom_area() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Lag (MB)") + ylab("ACF avg over chromosomes") + 
    geom_hline(yintercept = acf.hsc.sum.lst$acf.ci, linetype = "dashed", col = "blue") + 
    geom_hline(yintercept = -acf.hsc.sum.lst$acf.ci, linetype = "dashed", col = "blue") + 
    ggtitle(paste(jmark, trajname, "Prog Genome-wide"))
  print(m.acf.hsc.gw)
  
  acf.myeloid.sum.lst <- SummarizeACF(acf.out.myeloid.lst)
  m.acf.myeloid.gw <- ggplot(acf.myeloid.sum.lst$acf.out.sum, aes(x = dx / 10^6, y = acfval)) + geom_area() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Lag (MB)") + ylab("ACF avg over chromosomes") + 
    geom_hline(yintercept = acf.myeloid.sum.lst$acf.ci, linetype = "dashed", col = "blue") + 
    geom_hline(yintercept = -acf.myeloid.sum.lst$acf.ci, linetype = "dashed", col = "blue") + 
    ggtitle(paste(jmark, trajname, "Diff Genome-wide"))
  print(m.acf.myeloid.gw)

  multiplot(m.acf.hsc.gw, m.acf.myeloid.gw, cols = 1)
  
  
  # Plot the median log2 fold change relative to HSC cell: for one chromo
  
  jsub.ref.merge <- lapply(grep.strs, function(gstr) GetDiffRelToCell(imputed.dat, gstr = gstr, trajs, trajname = trajname, dat.umap.long = dat.umap.long, jmark = jmark)) %>%
    bind_rows() 
  m.mad <- ggplot(jsub.ref.merge, aes(x = umap1, y = umap2, color = exprs.diff.med)) + geom_point(size = jsize.facet) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.fc, limits = lims.fc) + facet_wrap(~label) + ggtitle(jmark, "Median log2 fold change with HSC")
  print(m.mad)
  
  # do genome-wide?
  jsub.ref.merge.gw <- lapply(c(""), function(gstr) GetDiffRelToCell(imputed.dat, gstr = gstr, trajs, trajname = trajname, dat.umap.long = dat.umap.long, jmark = jmark)) %>%
    bind_rows() 
  m.mad.gw <- ggplot(jsub.ref.merge.gw, aes(x = umap1, y = umap2, color = exprs.diff.med)) + geom_point(size = gw.jsize.facet) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.fc, limits = lims.fc) + facet_wrap(~label) + ggtitle(jmark, "Median log2 fold change with HSC")
  print(m.mad.gw)
  
  print(range(jsub.ref.merge.gw$exprs.diff.med))
  # plot along pseudotime? 
  traj.sub <- trajs[[jmark]][[trajname]]
  # add exprs.diff.med
  traj.sub <- left_join(traj.sub, jsub.ref.merge.gw %>% dplyr::select(cell, exprs.diff.med), by = c("cell"))
  
  m.mad.traj <- ggplot(traj.sub, aes(x = lambda, y = exprs.diff.med)) + geom_point(alpha = 0.1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Pseudotime") + ylab("Median Log2 FC from Prog Cell") + 
    ggtitle(jmark, paste(trajname, "Genome-wide"))
  print(m.mad.traj)

  m.mad.traj.fixscale <- ggplot(traj.sub, aes(x = lambda, y = exprs.diff.med)) + geom_point(alpha = 0.1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ylim(lims.fc) + 
    xlab("Pseudotime") + ylab("Median Log2 FC from Prog Cell") + 
    ggtitle(jmark, paste(trajname, "Genome-wide"))
  print(m.mad.traj.fixscale)

  if (!is.null(pdfout)){
  	dev.off()
  }
  
  # return something that can be used later for integrated analysis 
  return(jsub.ref.merge.gw)
}

GetDiffRelToCell2 <- function(jmerge.long, cell.ref){
  # add reference cell expression as exprs.ref column, then calculate median log2 fold change
  lambda.dat <- jmerge.long %>%
    group_by(cell) %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    dplyr::select(cell, lambda)
  
  jsub <- jmerge.long %>% filter(cell == cell.ref) %>%
    dplyr::select(coord, mark, ctype, exprs) %>%
    dplyr::rename(exprs.ref = exprs)
  
  jmerge.med.diff <- left_join(jmerge.long, jsub) %>%
    group_by(mark, ctype, cell) %>%
    mutate(exprs.diff = exprs - exprs.ref) %>% 
    summarise(exprs.diff.med = median(abs(exprs.diff))) %>%
    left_join(., lambda.dat)  # add back lambda info, which is lost
  return(jmerge.med.diff)
}



GetDiffRelToCell <- function(imputed.dat, gstr, trajs, trajname, dat.umap.long, jmark){
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

TotalVar <- function(x){
  # sum of squares
  return((x - mean(x)) ^ 2)
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
  assertthat::assert_that(!is.null(names(acf.out.hsc.lst)))
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

