# Jake Yeung
# Date of Creation: 2019-04-11
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/variance_analysis_pretty_figures.R
# Finalize variance analysis with pretty figures

rm(list=ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(JFuncs)
library(tidyr)
library(GGally)
library(ggrastr)
library(gridExtra)

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")



# Load data  --------------------------------------------------------------

# inf.trajs <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-09.RData"
inf.trajs <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-11.RData"
assertthat::assert_that(file.exists(inf.trajs))
load(inf.trajs, v=T)

head(trajs.spring[[4]][[1]])


inf.dat <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
assertthat::assert_that(file.exists(inf.dat))
load(inf.dat, v=T)


# Get constants -----------------------------------------------------------


make.plots <- TRUE
make.plots <- FALSE

colvec <- c("gray85", "gray50", "blue")  
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0")

nsecs <- 5
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

# do.plots <- TRUE

jstep <- 20000
# jtype <- "correlation"

pos.max <- 50 * 10^6
jstep <- 20000
lagmax <- pos.max / jstep

ci.interval <- 1.96  # corresponds to zscore 1.96 = alpha = 0.05 = 95% confidence interval


jsize <- 2


jstr <- "chr15:"

jthres <- 0.05

# Explore the data  -------------------------------------------------------

jmark <- "H3K4me1"
jmark <- "H3K27me3"
jmark <- "H3K9me3"

# Prepare dat -------------------------------------------------------------

jpseudo <- 0
jfac <- 10^7


# imputed.dat <- log2(t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics) + jpseudo) * jfac
# mat.sub <- MatToLong(imputed.dat, gstr = "chr15:", cells.vec = NULL) %>% dplyr::select(-start, -end)

mat.sub.merge <- lapply(jmarks, function(jmark) GetMatSub(tm.result.lst, jmark, jstr, jpseudo, jfac) %>% mutate(mark = jmark)) %>% 
  bind_rows()




# Explore -----------------------------------------------------------------

if (make.plots){
  pdf(paste0("~/data/scchic/pdfs/variance_over_pseudotime_plots_primetime.pseudo_", jpseudo, ".fac_", jfac, ".", Sys.Date(), ".pdf"), useDingbats = FALSE)
}

# out <- MakeVariabilityPlots(jmark, jname, tm.result.lst, dat.umap.long.trajs, pdfout = pdfout)


# get spatial distirbution for chromosome 15 for 3 cells along each trajectory 


# plot UMAP with pseudotime 

# plot one mark across the 3 trajectories

ctypes <- c("granu", "eryth", "lymphoid")
names(ctypes) <- ctypes

for (jmark in jmarks){
  
  m <- ggplot(dat.trajs.long %>% filter(mark == jmark), aes(x = X1, y = X2)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    # scale_color_manual(values = cbPalette) + 
    geom_path(data = trajs.spring[[jmark]]$granu %>% mutate(X1 = umap1, X2 = umap2), inherit.aes = FALSE, aes(x = X1, y = X2, color = lambda), size = jsize)  + 
    geom_path(data = trajs.spring[[jmark]]$eryth %>% mutate(X1 = umap1, X2 = umap2), inherit.aes = FALSE, aes(x = X1, y = X2, color = lambda), size = jsize)  + 
    geom_path(data = trajs.spring[[jmark]]$lymphoid %>% mutate(X1 = umap1, X2 = umap2), inherit.aes = FALSE, aes(x = X1, y = X2, color = lambda), size = jsize)  
  print(m)
  
  cells.lst <- lapply(ctypes, function(ctype) GetCellsAlongTraj(trajs.spring, jmark, ctype, n.sections = 5, thres = 0.2, show.plots = FALSE))
  jmerge.long <- lapply(ctypes, function(ctype) left_join(trajs.spring[[jmark]][[ctype]] %>% dplyr::select(cell, lambda), mat.sub.merge) %>% mutate(ctype = ctype, mark = jmark)) %>%
    bind_rows() 
  
  jsub <- dat.trajs.long %>% filter(mark == jmark)
  # plot location of cells on umap 
  m <- ggplot(mapping = aes(x = X1, y = X2)) + 
    geom_point(data = jsub, size = 1, color = "gray80", alpha = 0.5) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_path(data = trajs.spring[[jmark]]$granu %>% mutate(X1 = umap1, X2 = umap2), inherit.aes = FALSE, aes(x = X1, y = X2), size = jsize)  + 
    geom_path(data = trajs.spring[[jmark]]$eryth %>% mutate(X1 = umap1, X2 = umap2), inherit.aes = FALSE, aes(x = X1, y = X2), size = jsize)  + 
    geom_path(data = trajs.spring[[jmark]]$lymphoid %>% mutate(X1 = umap1, X2 = umap2), inherit.aes = FALSE, aes(x = X1, y = X2), size = jsize)   + 
    geom_point(data = jsub %>% filter(cell %in% unlist(cells.lst)), size = 3, color = "blue", alpha = 1) + 
    ggtitle("Location of cells for all 3 trajectories")
  print(m)
  
  mlst <- lapply(seq(length(cells.lst)), function(i){
    cells <- cells.lst[[i]]
    jname <- names(cells.lst)[[i]]
    m <- ggplot(mapping = aes(x = X1, y = X2)) + 
      geom_point(data = jsub, size = 1, color = "gray80", alpha = 0.5) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      geom_path(data = trajs.spring[[jmark]]$granu %>% mutate(X1 = umap1, X2 = umap2), inherit.aes = FALSE, aes(x = X1, y = X2), size = jsize)  + 
      geom_path(data = trajs.spring[[jmark]]$eryth %>% mutate(X1 = umap1, X2 = umap2), inherit.aes = FALSE, aes(x = X1, y = X2), size = jsize)  + 
      geom_path(data = trajs.spring[[jmark]]$lymphoid %>% mutate(X1 = umap1, X2 = umap2), inherit.aes = FALSE, aes(x = X1, y = X2), size = jsize)   + 
      geom_point(data = jsub %>% filter(cell %in% cells), size = 3, color = "blue", alpha = 1) + 
      ggtitle(paste("Location of cells for:", jname))
  })
  print(mlst)
  
  jmerge.sub <- jmerge.long %>%
    filter(cell %in% unlist(cells.lst)) %>%
    group_by(mark, coord, ctype) %>%
    mutate(lambda.rank = rank(lambda))
  
  # plot chromosome 15 over time
  m <- ggplot(jmerge.sub, aes(x = pos / 10^6, y = exprs, group = lambda.rank, color = lambda)) + 
    geom_line() + theme_bw() + facet_grid(lambda.rank ~ ctype) + 
    theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
    ggtitle(jmark)  +  
    xlab("MB") + 
    ylab("log2(Counts Per Million)") + 
    scale_color_gradient2(low = colvec[[1]], mid = colvec[[2]], high = scales::muted(colvec[[3]]), lim = c(0, 1), midpoint = 0.5)
  print(m)
  
  # calculate variance over time?
  take.nth.cell <- 3
  jmerge.med.diff.lst <- lapply(ctypes, function(jctype){
    cell.ref.lst <- jmerge.long %>%
      filter(ctype == jctype) %>%
      group_by(cell) %>%
      filter(row_number() == 1) %>%
      arrange(lambda) %>%
      dplyr::select(cell)
    cell.ref.lst <- cell.ref.lst$cell[[take.nth.cell]]
    print(cell.ref.lst)
    # cell.ref.lst <- unique(subset(jmerge.long, lambda > 0.25 & lambda < 0.3 & ctype == jctype)$cell)[[1]]
    jmerge.diff.lst <- lapply(cell.ref.lst, function(cell.ref){
      jtmp <- GetDiffRelToCell2(jmerge.long %>% filter(ctype == jctype), cell.ref)
      jtmp$cell.ref <- cell.ref
      return(jtmp)
    }) %>%
      bind_rows()
  }) %>%
    bind_rows()
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0")
  m.exprs.diff.med <- ggplot(jmerge.med.diff.lst, aes(x = lambda, y = exprs.diff.med, color = ctype)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark) + scale_color_manual(values = cbPalette)
  print(m.exprs.diff.med)
  m.exprs.diff.med2 <- ggplot(jmerge.med.diff.lst, aes(x = lambda, y = exprs.diff.med, color = ctype)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark) + scale_color_manual(values = cbPalette) + facet_wrap(~ctype)
  print(m.exprs.diff.med2)
}


  


# plot one trajectoory across the 4 marks

ctype <- "granu"

for (ctype in ctypes){
  cells.lst <- lapply(jmarks, function(jmark) GetCellsAlongTraj(trajs.spring, jmark, ctype, n.sections = nsecs, thres = jthres, show.plots = FALSE))
  jmerge.long <- lapply(jmarks, function(jmark) left_join(trajs.spring[[jmark]][[ctype]] %>% dplyr::select(cell, lambda), mat.sub.merge) %>% mutate(ctype = ctype, mark = jmark)) %>%
    bind_rows()
  
  # check points 
  
  mlst <- lapply(seq(length(cells.lst)), function(i){
    cells <- cells.lst[[i]]
    jmark <- names(cells.lst)[[i]]
    jsub <- dat.trajs.long %>% filter(mark == jmark)
    m <- ggplot(mapping = aes(x = X1, y = X2)) + 
      geom_point(data = jsub, size = 1, color = "gray80", alpha = 0.5) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      geom_path(data = trajs.spring[[jmark]][[ctype]] %>% mutate(X1 = umap1, X2 = umap2), inherit.aes = FALSE, aes(x = X1, y = X2), size = jsize)   + 
      geom_point(data = jsub %>% filter(cell %in% cells), size = 3, color = "blue", alpha = 1) + 
      ggtitle(paste("Location of cells for:", jmark))
  })
  print(mlst)
  
  
  # add rank
  jmerge.filt <- jmerge.long %>%
    filter(cell %in% unlist(cells.lst)) %>%
    group_by(mark, coord, ctype) %>%
    mutate(lambda.rank = rank(lambda)) 
  
  jmerge.filt$mark <- factor(as.character(jmerge.filt$mark), levels = jmarks)
  
  m <- ggplot(jmerge.filt,
              aes(x = pos / 10^6, y = exprs, group = lambda.rank, color = lambda)) + 
    geom_line() + theme_bw() + facet_grid(lambda.rank ~ mark) + 
    theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
    ggtitle(ctype)  +  
    xlab("MB") + 
    ylab("log2(Counts Per Million)") + 
    scale_color_gradient2(low = colvec[[1]], mid = colvec[[2]], high = scales::muted(colvec[[3]]), lim = c(0, 1), midpoint = 0.5)
  print(m)
  
  # do it on select cells, but whole genome
  mat.gw.cellfilt <- mapply(function(jmark, jcells.vec) GetMatSub(tm.result.lst, jmark, "", jpseudo, jfac, cells.vec = jcells.vec) %>% mutate(mark = jmark) %>% bind_rows(), jmarks, cells.lst, SIMPLIFY = FALSE)
  
  mat.gw.cellfilt <- mat.gw.cellfilt %>% bind_rows() 
  
  jmerge.gw.long <- lapply(jmarks, function(jmark) left_join(trajs.spring[[jmark]][[ctype]] %>% filter(cell %in% cells.lst[[jmark]]) %>% dplyr::select(cell, lambda), 
                                                             mat.gw.cellfilt) %>% mutate(ctype = ctype, mark = jmark)) %>%
    bind_rows() %>%
    group_by(mark, coord, ctype, pos, chromo) %>%
    mutate(lambda.rank = rank(lambda))
  
  # correlation across the marks ? 
  jmerge.wide <- tidyr::spread(jmerge.gw.long %>% 
                          ungroup() %>%
                          mutate(markrank = paste(mark, lambda.rank, sep = "_")) %>%
                          dplyr::select(-lambda, -lambda.rank, -cell, -mark),
                        key = markrank, value = exprs)
  print(jmerge.wide)
  lambda.ranks <- seq(nsecs)
  
  jstrs <- lapply(lambda.ranks, function(l) return(paste0("^H3.*_", l)))
  jplots <- lapply(jstrs, function(jstr){
    cols.keep <- grep(jstr, colnames(jmerge.wide))
    m <- ggpairs(jmerge.wide, columns = cols.keep, 
                 # lower = list(continuous = wrap("points", alpha = 0.05))) + 
                 lower = list(continuous = MakeLowerPlot)) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                         axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
      ggtitle(paste(ctype, "genomewide development stage", substr(jstr, nchar(jstr), nchar(jstr))))
    return(m)
  })
  print(jplots)
}

dev.off()
