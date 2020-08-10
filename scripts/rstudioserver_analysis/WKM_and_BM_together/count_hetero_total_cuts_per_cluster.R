# Jake Yeung
# Date of Creation: 2020-06-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/count_cells_per_cluster.R
# Count cells per cluster 


rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)


# Indir -------------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

indir.clst <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/cell_to_cluster_tables_merged"
indir.counts <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/offsets_hetero_and_totalcuts"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/offsets_hetero_and_totalcuts.countsByCluster"
dir.create(outdir)
# outf <- file.path(outdir, paste0("hetero_and_totalcuts_by_cluster.txt"))
outpdf <- file.path(outdir, paste0("hetero_and_totalcuts_by_cluster.pdf"))

pdf(outpdf, useDingbats = FALSE)

cell.clusters <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(indir.clst, paste0("BM_cell_to_clusters.", jmark, ".txt"))
  print(inf)
  clst <- fread(inf)
  return(clst)
})

cell.summary <- lapply(jmarks, function(jmark){
  jdat <- cell.clusters[[jmark]]
  jdat <- jdat %>% 
    group_by(cluster) %>%
    summarise(ncell = length(cell))
  jdat$mark <- jmark
  return(jdat)
})


load(file.path(indir.counts, "MouseBM_HeteroTotalCounts_50kb_bins.2020-06-16.RData"), v=T)


# Add clusters ------------------------------------------------------------

dat.merge.lst <- lapply(jmarks, function(jmark){
  dat.merge <- left_join(dat.ncuts.hetero.total[[jmark]], cell.clusters[[jmark]]) %>%
    left_join(., cell.summary[[jmark]])
})

dat.merge.filt.lst <- lapply(dat.merge.lst, function(jdat){
  subset(jdat, !is.na(cluster))
})

m.lst <- lapply(jmarks, function(jmark){
  ggplot(dat.merge.filt.lst[[jmark]], aes(x = cluster, y = ncuts.total)) + 
    geom_boxplot() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark) + scale_y_log10()
})
print(m.lst)

m.lst <- lapply(jmarks, function(jmark){
  ggplot(dat.merge.filt.lst[[jmark]], aes(x = cluster, y = ncuts.hetero)) + 
    geom_boxplot() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark) + scale_y_log10()
})
print(m.lst)

m.lst <- lapply(jmarks, function(jmark){
  ggplot(dat.merge.filt.lst[[jmark]], aes(x = cluster, y = ncuts.hetero / ncuts.total)) + 
    geom_boxplot() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark) + scale_y_log10()
})
print(m.lst)

# Pseudobulk avgs ---------------------------------------------------------

dat.merge.pbulk.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.merge.filt.lst[[jmark]]
  jdat %>%
    group_by(cluster) %>%
    summarise(ncuts.total = sum(ncuts.total),
              ncuts.hetero = sum(ncuts.hetero)) %>%
    left_join(., cell.summary[[jmark]])
})


m.pbulk.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.merge.pbulk.lst[[jmark]]
  ggplot(jdat, aes(x = cluster, y = ncuts.total)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_col() + 
    ggtitle(jmark)
})
print(m.pbulk.lst)

m.pbulk.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.merge.pbulk.lst[[jmark]]
  ggplot(jdat, aes(x = cluster, y = ncuts.hetero)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_col() + 
    ggtitle(jmark)
})
print(m.pbulk.lst)

m.pbulk.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.merge.pbulk.lst[[jmark]]
  ggplot(jdat, aes(x = cluster, y = ncuts.hetero / ncuts.total)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_col() + 
    ggtitle(jmark)
})
print(m.pbulk.lst)

dev.off()


# Write table  ------------------------------------------------------------

# outf <- file.path(outdir, paste0("hetero_and_totalcuts_by_cluster.txt"))
lapply(jmarks, function(jmark){
  jdat <- dat.merge.pbulk.lst[[jmark]]
  jclsts <- jdat$cluster; names(jclsts) <- jclsts
  lapply(jclsts, function(jclst){
    outf <- file.path(outdir, paste0("hetero_and_totalcuts.", jclst, ".", jmark, ".txt"))
    outf2 <- file.path(outdir, paste0("hetero_and_totalcuts.", jclst, ".", jmark, ".NoCname.txt"))
    fwrite(jdat %>% filter(cluster == jclst), file = outf, sep = "\t")
    fwrite(jdat %>% filter(cluster == jclst), file = outf2, sep = "\t", col.names = FALSE)
  })
})

# 
# N <- seq(1000)
# 
# y1 <- 2 / (2 * N)
# y2 <- 3 / (4 * N - 1)
# 
# plot(N, y1, type = "l")
# lines(N, y2, col = 'blue')
# 
# plot(y2 - y1)
# 

