# Jake Yeung
# Date of Creation: 2019-04-03
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/variance_over_pseudotime_plots.R
# Plot the outputs 

rm(list=ls())

tstart <- Sys.time()

library(dplyr)
library(ggplot2)
library(tidyr)
library(JFuncs)

source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/FourierFunctions.R")
source("scripts/Rfunctions/VariabilityFunctions.R")

MakeVariabilityPlots <- function(jmark, trajname, tm.result.lst, dat.umap.long.trajs, 
                                 jcol = c("gray80", "gray50", "darkblue"), grep.strs = paste("chr", c(seq(21)), ":", sep = ""), jalpha = 0.5, pseudo = 0, jscale = 1, 
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
    theme_bw() + theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    xlab("MB") + ylab("log2(imputed counts)") + ggtitle(paste(jmark, trajname, "Prog Cell"))
  print(m.hsc.chromo.all)
  
  m.myeloid.chromo.all <- ggplot(jsub.myeloid.lst %>% bind_rows(), aes(x = pos / 10^6, y = log2(exprs))) + geom_line() + facet_wrap(~chromo, scales = "free_x", ncol = 7) + 
    theme_bw() + theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    xlab("MB") + ylab("log2(imputed counts)") + ggtitle(paste(jmark, trajname, "Diff Cell"))
  print(m.myeloid.chromo.all)
  

  acf.out.hsc.lst <- lapply(jsub.hsc.lst, function(jsub.hsc) CalculateACF(jsub.hsc, jstep = jstep, jtype = "correlation", jmain = paste(jmark, trajname, "Prog Cell", gstr), show.plot = FALSE))
  acf.out.myeloid.lst <- lapply(jsub.myeloid.lst, function(jsub.hsc) CalculateACF(jsub.hsc, jstep = jstep, jtype = "correlation", jmain = paste(jmark, trajname, "Diff Cell", gstr), show.plot = FALSE))
  
  # average out the plots for different lags 
  # plot all chromosomes
  par(mfrow=c(3,7), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0, pty = "s")
  for (i in seq(length(acf.out.hsc.lst))){
    plot(acf.out.hsc.lst[[i]]$lag.stepadj / 10^6, acf.out.hsc.lst[[i]]$acf, main = paste(jmark, trajname, "Prog", grep.strs[[i]]), type = "h", xlab = "Step Size (MB)", ylab = "Autocorrelation")
    abline(h = ci.interval / sqrt(length(acf.out.hsc.lst[[i]]$lag.stepadj)), lty = "dotted", col = "blue")
    abline(h = -ci.interval / sqrt(length(acf.out.hsc.lst[[i]]$lag.stepadj)), lty = "dotted", col = "blue")
    abline(h = 0, cex = 2)
  }
  for (i in seq(length(acf.out.myeloid.lst))){
    plot(acf.out.myeloid.lst[[i]]$lag.stepadj / 10^6, acf.out.hsc.lst[[i]]$acf, main = paste(jmark, trajname, "Diff", grep.strs[[i]]), type = "h", xlab = "Step Size (MB)", ylab = "Autocorrelation")
    abline(h = ci.interval / sqrt(length(acf.out.hsc.lst[[i]]$lag.stepadj)), lty = "dotted", col = "blue")
    abline(h = -ci.interval / sqrt(length(acf.out.hsc.lst[[i]]$lag.stepadj)), lty = "dotted", col = "blue")
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

  if (!is.null(pdfout)){
  	dev.off()
  }
  
  # return something that can be used later for integrated analysis 
  return(jsub.ref.merge.gw)
}

# Load constants ----------------------------------------------------------


jcol <- c("gray80", "gray50", "darkblue")
grep.strs <- paste("chr", c(seq(21)), ":", sep = "")

jalpha <- 0.5

pseudo <- 0
jscale <- 1

mdpt.sd <- 1
lims.sd <- c(0, 3)
mdpt.fc <- 0.75
lims.fc <- c(0, 3)

jsize.facet <- 0.2
gw.jsize.facet <- 2

do.plots <- TRUE

jstep <- 20000
# jtype <- "correlation"

pos.max <- 50 * 10^6
jstep <- 20000
lagmax <- pos.max / jstep

ci.interval <- 1.96  # corresponds to zscore 1.96 = alpha = 0.05 = 95% confidence interval

# Load data ---------------------------------------------------------------

inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.RData"

load(inf, v=T)

if (do.plots){
  pdf("~/data/scchic/pdfs/variance_over_pseudotime_plots.pdf", useDingbats = FALSE)
}


# Do H3K27me3 -------------------------------------------------------------

# jmark <- "H3K27me3"
# trajname <- "myeloid"

# print(paste(jmark, trajname))

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
trajnames <- c("bcell", "myeloid", "myeloid", "myeloid")
outmain <- "/Users/yeung/data/scchic/pdfs"
pdfsouts <- lapply(jmarks, function(jmark) file.path(outmain, paste0("variance_over_pseudotime_plots.", jmark, ".pdf")))

# test on one 
# out <- MakeVariabilityPlots(jmarks[[1]], trajnames[[1]], tm.result.lst, dat.umap.long.trajs)

out <- parallel::mcmapply(function(jmark, trajname, pdfout) MakeVariabilityPlots(jmark, trajname, tm.result.lst, dat.umap.long.trajs, pdfout = pdfout), jmarks, trajnames, pdfouts, SIMPLIFY = FALSE, mc.cores = 4)
# out <- mapply(function(jmark, trajname) MakeVariabilityPlots(jmark, trajname, tm.result.lst, dat.umap.long.trajs), jmarks, trajnames, SIMPLIFY = FALSE)

print(out)

if (do.plots){
  dev.off()
}



print(Sys.time() - tstart)
