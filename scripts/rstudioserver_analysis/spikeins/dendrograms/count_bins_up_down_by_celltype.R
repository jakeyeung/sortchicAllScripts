# Jake Yeung
# Date of Creation: 2021-03-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/count_bins_up_down_by_celltype.R
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# Load metas --------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps/count_bins_up_down_cellfate_independent.", Sys.Date(), ".pdf")

make.plots <- FALSE

if (make.plots){
  pdf(outpdf, useDingbats = FALSE)
}

indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage")

dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  fread(file.path(indir.meta, fname)) %>%
    rowwise() 
}) 

jmarktest <- "H3K4me1"
cluster2col <- hash(dat.metas[[jmarktest]]$cluster, dat.metas[[jmarktest]]$clustercol)


# Load dynoamic bins ------------------------------------------------------


infs.de.lst <- lapply(jmarks, function(jmark){
  inf.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.jmidbug_fixed.", jmark, ".2021-02-19.txt")
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

dat.de.bins.lst <- lapply(infs.de.lst, function(jinf){
  fread(jinf)
})


# Get mat -----------------------------------------------------------------

inf.rds <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps/variance_by_clusters_by_HighBins.2021-03-09.1000000.log10filt_-1.methods.check_other_norms.add_DE_bins.rds"
dat.mat.lst <- readRDS(inf.rds)



# Get celltype-specific bins ----------------------------------------------


dat.mat.log2.diff.lst <- lapply(dat.mat.lst, function(dat.mat){
  dat.mat.log2 <- log2(dat.mat)
  dat.mat.log2.diff <- sweep(dat.mat.log2, MARGIN = 1, STATS = dat.mat.log2[, "HSPCs"], FUN = "-")
  cnames.keep <- colnames(dat.mat.log2.diff) != "HSPCs"
  return(dat.mat.log2.diff[, cnames.keep])
})

ctypes.k4 <- colnames(dat.mat.log2.diff.lst$H3K4me1)[colnames(dat.mat.log2.diff.lst$H3K4me1) != "HSPCs"]
names(ctypes.k4) <- ctypes.k4

ctypes.k9 <- colnames(dat.mat.log2.diff.lst$H3K9me3)[colnames(dat.mat.log2.diff.lst$H3K9me3) != "HSPCs"]
names(ctypes.k9) <- ctypes.k9

ctype.spec.up <- lapply(jmarks, function(jmark){
  if (jmark != "H3K9me3"){
    ctypes <- ctypes.k4
  } else {
    ctypes <- ctypes.k9
  }
  print(jmark)
  print(ctypes)
  lapply(ctypes, function(ctype){
    cnames.keep <- colnames(dat.mat.log2.diff.lst[[jmark]]) != ctype
    jmat <- dat.mat.log2.diff.lst[[jmark]][, cnames.keep]
    jref <- dat.mat.log2.diff.lst[[jmark]][, ctype]
    ctype.up.any <- jref > 0
    others.down.all <- apply(jmat, MARGIN = 1, FUN = function(jrow) all(jrow < 0))
    ctype.up.vec <- ctype.up.any & others.down.all
    return(names(which(ctype.up.vec)))
  })
})

ctype.spec.down <- lapply(jmarks, function(jmark){
  if (jmark != "H3K9me3"){
    ctypes <- ctypes.k4
  } else {
    ctypes <- ctypes.k9
  }
  print(jmark)
  lapply(ctypes, function(ctype){
    cnames.keep <- colnames(dat.mat.log2.diff.lst[[jmark]]) != ctype
    jmat <- dat.mat.log2.diff.lst[[jmark]][, cnames.keep]
    jref <- dat.mat.log2.diff.lst[[jmark]][, ctype]
    ctype.down.any <- jref < 0
    others.up.all <- apply(jmat, MARGIN = 1, FUN = function(jrow) all(jrow > 0))
    ctype.down.vec <- ctype.down.any & others.up.all
    return(names(which(ctype.down.vec)))
  })
})

for (jmark in jmarks){
  print(jmark)
  print(lapply(ctype.spec.up[[jmark]], length))
}


for (jmark in jmarks){
  print(jmark)
  print(lapply(ctype.spec.down[[jmark]], length))
}


# Get mats ----------------------------------------------------------------


dat.long.log.lst <- lapply(dat.mat.lst, function(dat.mat){
  dat.mat <- log2(dat.mat)
  hspc.col <- colnames(dat.mat) == "HSPCs"
  dat.mat.hspcs <- data.frame(bin = rownames(dat.mat), reference.exprs = dat.mat[, hspc.col], stringsAsFactors = FALSE)
  dat.mat.otherhspcs <- dat.mat[, !hspc.col]
  jdat <- dat.mat.otherhspcs %>%
    melt()
  colnames(jdat) <- c("bin", "celltype", "log2exprs")
  jdat.merge <- left_join(jdat, dat.mat.hspcs)
  return(jdat.merge)
})  

dat.long.log.filt.lst <- lapply(jmarks, function(jmark){
  subset(dat.long.log.lst[[jmark]], bin %in% dat.de.bins.lst[[jmark]]$CoordOriginal)
})


# Count up or down relative to HSPCs --------------------------------------

dat.long.log.filt.counts.lst <- lapply(dat.long.log.filt.lst, function(jdat){
  jdat %>%
    rowwise() %>%
    mutate(is.up = log2exprs - reference.exprs > 0) %>%
    group_by(bin) %>%
    mutate(all.up = all(is.up),
           all.down = all(!is.up))
})


dat.sum <- lapply(jmarks, function(jmark){
  jdat <- dat.long.log.filt.counts.lst[[jmark]]
  jdat <- jdat %>%
    group_by(celltype, is.up) %>%
    summarise(bincount = length(bin)) %>%
    mutate(mark = jmark) %>%
    group_by(celltype)  %>%
    mutate(binfrac = bincount / sum(bincount)) %>%
    rowwise() %>%
    mutate(bincount = ifelse(!is.up, -bincount, bincount)) 
}) %>%
  bind_rows()


dat.sum.all.up <- lapply(jmarks, function(jmark){
  jdat <- dat.long.log.filt.counts.lst[[jmark]]
  jdat <- jdat %>%
    group_by(celltype, all.up) %>%
    summarise(bincount.all = length(bin)) %>%
    mutate(mark = jmark) 
}) %>%
  bind_rows() %>%
  filter(all.up)


dat.sum.all.down <- lapply(jmarks, function(jmark){
  jdat <- dat.long.log.filt.counts.lst[[jmark]]
  jdat <- jdat %>%
    group_by(celltype, all.down) %>%
    summarise(bincount.all = length(bin)) %>%
    mutate(mark = jmark) 
}) %>%
  bind_rows() %>%
  filter(all.down)

dat.sum.all.merge <- rbind(subset(dat.sum.all.up, all.up, select = -all.up) %>% mutate(status = "AllUp"), 
                           subset(dat.sum.all.down, all.down, select = -all.down) %>% mutate(status = "AllDown"))


ggplot(dat.sum, aes(x = celltype, y = bincount, fill = is.up)) + 
  geom_col() + 
  theme_bw() + 
  facet_wrap(~mark, scales = "free") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.sum.all.merge, aes(x = celltype, y = bincount.all, fill = status)) + 
  geom_col() + 
  theme_bw() + 
  facet_wrap(~mark, scales = "free") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# What fraction of bins are up in all cell types?  ------------------------

jcheck.up <- left_join(subset(dat.sum, is.up), dat.sum.all.up) %>%
  mutate(status = "up") %>%
  dplyr::select(-all.up)
jcheck.down <- left_join(subset(dat.sum, !is.up), dat.sum.all.down) %>%
  mutate(status = "down") %>%
  dplyr::select(-all.down) 

jcheck.down.rev <- left_join(subset(dat.sum, !is.up), dat.sum.all.down) %>%
  mutate(status = "down") %>%
  dplyr::select(-all.down) %>%
  mutate(bincount = -bincount)

jcheck.merge <- rbind(jcheck.up, jcheck.down) %>%
  group_by(celltype, mark) 

jcheck.merge.rev <- rbind(jcheck.up, jcheck.down.rev) %>%
  group_by(celltype, mark)  %>%
  summarise(bincount = sum(bincount),
            bincount.all = sum(bincount.all))
  
jcheck.merge.rev$clustercol <- sapply(jcheck.merge.rev$celltype, function(x) cluster2col[[x]])
jcheck.merge$mark <- factor(jcheck.merge$mark, levels = jmarks)
jcheck.merge.rev$mark <- factor(jcheck.merge.rev$mark, levels = jmarks)
jcheck.merge$clustercol <- sapply(jcheck.merge$celltype, function(x) cluster2col[[x]])


ggplot(jcheck.up, aes(y = bincount.all / bincount, x = celltype, fill = mark)) + 
  facet_wrap(~mark) +
  geom_col() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

ggplot(jcheck.merge, aes(y = bincount.all / bincount, x = celltype, fill = clustercol)) + 
  scale_fill_identity() + 
  facet_wrap(~mark, scales = "free_x", nrow = 1) +
  geom_col() + 
  theme_bw() + 
  xlab("") + 
  geom_hline(yintercept = 0) + 
  ylab("Fraction of bins that change in\ncelltype-independent manner") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")


ggplot(jcheck.merge, aes(y = bincount, x = celltype, fill = clustercol)) + 
  scale_fill_identity() + 
  facet_wrap(~mark, scales = "free", nrow = 1) +
  geom_col() + 
  theme_bw() + 
  xlab("") + 
  geom_hline(yintercept = 0) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")


ggplot(jcheck.merge, aes(y = bincount.all, x = celltype, fill = clustercol)) + 
  scale_fill_identity() + 
  facet_grid(status~mark, scales = "free_x") +
  geom_col() + 
  theme_bw() + 
  xlab("") + 
  geom_hline(yintercept = 0) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")


ggplot(jcheck.merge.rev, aes(y = bincount.all / bincount, x = celltype, fill = clustercol)) + 
  scale_fill_identity() + 
  facet_wrap(~mark, scales = "free_x", nrow = 1) +
  geom_col() + 
  theme_bw() + 
  xlab("") + 
  ylab("Fraction of bins that change in\ncelltype-independent manner") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")



# Annotate each celltypespecific bin  -------------------------------------


dat.sum.bymark <- split(x = dat.sum, f = dat.sum$mark)

dat.long.log.filt.counts.lst.byctype <- lapply(dat.long.log.filt.counts.lst, function(jdat){
  split(x = jdat, f = jdat$celltype)
})

# annotate ctype up 
dat.long.log.filt.counts.lst.byctype.added <- lapply(jmarks, function(jmark){
  if (jmark != "H3K9me3"){
    ctypes <- ctypes.k4
  } else {
    ctypes <- ctypes.k9
  }
  print(jmark)
  print(ctypes)
  lapply(ctypes, function(jctype){
    ctype.cname <- paste0("ctypeup")
    ctype.cname2 <- paste0("ctypedown")
    jdat.annot <- dat.long.log.filt.counts.lst.byctype[[jmark]][[jctype]]  %>%
      rowwise() %>%
      mutate(ctypeup = bin %in% ctype.spec.up[[jmark]][[jctype]],
             ctypedown = bin %in% ctype.spec.down[[jmark]][[jctype]])
    cnames.out <- c("bin", ctype.cname, ctype.cname2)
    jdat.annot.sub <- jdat.annot[, cnames.out]
    jdat.joined <- left_join(dat.long.log.filt.counts.lst.byctype[[jmark]][[jctype]], jdat.annot.sub)
    jdat.joined$mark <- jmark
    return(jdat.joined)
  })
})




# Count for each celltype -------------------------------------------------

count.sum.long.up <- lapply(jmarks, function(jmark){
  if (jmark != "H3K9me3"){
    ctypes <- ctypes.k4
  } else {
    ctypes <- ctypes.k9
  }
  lapply(ctypes, function(jctype){
    counts.sum <- subset(dat.long.log.filt.counts.lst.byctype.added[[jmark]][[jctype]], is.up) %>%
      group_by(celltype, mark, is.up) %>% 
      summarise(count.any = length(which(is.up)),
                count.all = length(which(all.up)),
                count.ctype = length(which(ctypeup)))
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  mutate(status = "up")


count.sum.long.down <- lapply(jmarks, function(jmark){
  if (jmark != "H3K9me3"){
    ctypes <- ctypes.k4
  } else {
    ctypes <- ctypes.k9
  }
  lapply(ctypes, function(jctype){
    counts.sum <- subset(dat.long.log.filt.counts.lst.byctype.added[[jmark]][[jctype]], !is.up) %>%
      group_by(celltype, mark, is.up) %>% 
      summarise(count.any = length(which(!is.up)),
                count.all = length(which(all.down)),
                count.ctype = length(which(ctypedown)))
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  mutate(status = "down") %>%
  mutate(count.any = -count.any)
 

count.sum.long.merge <- rbind(count.sum.long.up, count.sum.long.down)
count.sum.long.merge$clustercol <- sapply(count.sum.long.merge$celltype, function(x) cluster2col[[x]])
count.sum.long.merge$mark <- factor(count.sum.long.merge$mark, levels = jmarks)


count.sum.long.up$clustercol <- sapply(count.sum.long.up$celltype, function(x) cluster2col[[x]])
count.sum.long.down$clustercol <- sapply(count.sum.long.down$celltype, function(x) cluster2col[[x]])

ggplot(count.sum.long.up, aes(x = celltype, y = count.all / count.any, fill = clustercol)) + 
  geom_col() + 
  facet_wrap(~mark, nrow = 1, scales = "free_x") + 
  scale_fill_identity() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")


ggplot(count.sum.long.down, aes(x = celltype, y = count.all / count.any, fill = clustercol)) + 
  geom_col() + 
  facet_wrap(~mark, nrow = 1, scales = "free_x") + 
  scale_fill_identity() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")



ggplot(count.sum.long.up, aes(x = celltype, y = count.ctype / count.any, fill = clustercol)) + 
  geom_col() + 
  facet_wrap(~mark, nrow = 1, scales = "free_x") + 
  scale_fill_identity() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")

ggplot(count.sum.long.down, aes(x = celltype, y = count.ctype / count.any, fill = clustercol)) + 
  geom_col() + 
  facet_wrap(~mark, nrow = 1, scales = "free_x") + 
  scale_fill_identity() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")


ggplot(count.sum.long.merge, aes(x = celltype, y = count.all / count.any, fill = clustercol)) + 
  geom_col() + 
  facet_wrap(~mark, nrow = 1, scales = "free_x") + 
  scale_fill_identity() + 
  geom_hline(yintercept= 0) + 
  ylab("Fraction of bins that are cell fate independent") + 
  ylim(c(-1, 1)) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")


ggplot(count.sum.long.merge, aes(x = celltype, y = count.ctype / count.any, fill = clustercol)) + 
  geom_col() + 
  facet_wrap(~mark, nrow = 1, scales = "free_x") + 
  scale_fill_identity() + 
  geom_hline(yintercept= 0) + 
  ylim(c(-1, 1)) + 
  ylab("Fraction of bins that are celltype specific") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")

if (make.plots){
  dev.off()
}

