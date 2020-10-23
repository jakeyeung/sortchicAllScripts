# Jake Yeung
# Date of Creation: 2020-09-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse/annotate_celltypes_from_glmpca.R
# Quick and dirty


rm(list=ls())

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

# Load glmpca -------------------------------------------------------------


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"
jsuffix <- "TSS.bsize_10000"



# Load data  --------------------------------------------------------------

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_mouse_spikein_VAN5046_varfilt"
assertthat::assert_that(dir.exists(indir))

jvarfilt <- 1
jmaxiter <- 1000
jtol <- as.character("1e-6")


dat.umap.lst <- lapply(jmarks, function(jmark){
  
  print(jmark)
  inf <- file.path(indir, paste0("count_mat_", jmark, "_l2r_filt.2020-09-13.minl2r_-1.varfilt_", jvarfilt, ".glmpcaout.penalty_1.maxiter_", jmaxiter, ".stochastic.avagrad.tol_", jtol, ".devfilt_5000.by_plate.RData"))
  print(basename(inf))
  load(inf, v=T)
  
  dat.umap.long <- DoUmapAndLouvain(glmpcaout$factors, jsettings = jsettings)
  
  dat.umap.long$dim1 <- glmpcaout$factors[, 1]
  dat.umap.long$dim2 <- glmpcaout$factors[, 2]
  dat.umap.long$dim3 <- glmpcaout$factors[, 3]
  dat.umap.long$dim4 <- glmpcaout$factors[, 4]
  
  dat.umap.long$mark <- jmark
  
  return(dat.umap.long)
})

cells.keep.lst <- lapply(dat.umap.lst, function(jdat) jdat$cell)


cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

lapply(jmarks, function(jmark){
  m <- ggplot(dat.umap.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_manual(values = cbPalette) +
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
})


# Assign louvain and mark -------------------------------------------------

k4me1.hash <- hash(list("7" = "Granus", "4" = "Granus", "1" = "Granus", "8" = "HSCs", "9" = "HSCs", "10" = "Bcells", "2" = "Eryths"))
k4me3.hash <- hash(list("6" = "Granus", "3" = "Granus", "2" = "Eryths", "4" = "Bcells", "5" = "HSCs", "5" = "HSCs"))
k27me3.hash <- hash(list("5" = "HSCs", "2" = "Bcells", "4" = "Bcells", "8" = "Granus", "6" = "Granus", "7" = "Eryths"))

hlst <- list("H3K4me1" = k4me1.hash, "H3K4me3" = k4me3.hash, "H3K27me3" = k27me3.hash)

dat.umap.annot.lst <- lapply(jmarks, function(jmark){
  dat.tmp <- dat.umap.lst[[jmark]] %>%
    rowwise() %>%
    mutate(ctype = AssignHash(x = as.character(louvain), hlst[[jmark]], null.fill = NA))
})

ggplot(dat.umap.annot.lst$H3K4me1, aes(x = umap1, y = umap2, color = ctype)) + 
  geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Add spikein/chromo by celltype ------------------------------------------


# Load sipkeins -----------------------------------------------------------

inf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.chromo2spikeinfilt.VAN5046/spikein_info.txt"
dat.spikeins.mat <- fread(inf.spikeins)

# dat.merge <- dat.spikeins.mat %>%
dat.merge <- left_join(dat.umap.annot.lst %>% bind_rows(), subset(dat.spikeins.mat, select = c(samp, spikeincounts, chromocounts)), by = c("cell" = "samp")) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"), 
         plate = as.numeric(strsplit(experi, "-")[[1]][[6]]),
         rowcoord = AddPlateCoordinates(cell)$rowcoord,
         colcoord = AddPlateCoordinates(cell)$colcoord,
         stype = AnnotateSortFromLayout(plate, rowcoord, colcoord))

# statistically significant? 

dat.merge$mark <- factor(dat.merge$mark, levels = jmarks)

ggplot(dat.merge %>% filter(!is.na(ctype)), aes(x = ctype, y = log2(chromocounts / spikeincounts))) + 
  geom_boxplot(outlier.color = "white") + 
  geom_jitter(width = 0.1, alpha = 0.25)  + 
  facet_wrap(~mark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge %>% filter(!is.na(ctype)), aes(x = ctype, y = log2(chromocounts))) + 
  geom_boxplot(outlier.color = "white") + 
  geom_jitter(width = 0.1, alpha = 0.25)  + 
  facet_wrap(~mark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Load count tables -------------------------------------------------------

indir.counts <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM/tagged_bams/merged_bams/TSS_count_tables.NewFilters")
assertthat::assert_that(dir.exists(indir.counts))
infs <- lapply(jmarks, function(jmark){
  fname <- paste0("PZ-ChIC-mouse-BM-", jmark, "-merged.sorted.tagged.countTable.", jsuffix, ".csv")
  # fname <- paste0("PZ-ChIC-mouse-BM-", jmark, "-merged.sorted.tagged.countTable.TSS.bsize_10000.csv")
  inf <- file.path(indir.counts, fname)
  assertthat::assert_that(file.exists(inf))
  return(inf)
})


mats.lst <- lapply(infs, function(inf) ReadMatTSSFormat(inf, as.sparse = TRUE, add.coord = FALSE, sort.rnames = FALSE))
mats.filt.lst <- lapply(jmarks, function(jmark) mats.lst[[jmark]][ , cells.keep.lst[[jmark]]])

cuts.gbodies <- lapply(mats.filt.lst, function(jmat) data.frame(cell = colnames(jmat), genecounts = colSums(jmat)))

jmat.bulk.long.lst <- lapply(jmarks, function(jmark){
  jcnames.keep.lst <- lapply(split(dat.merge, f = dat.merge$ctype), function(x) x$cell)
  jmat.pbulk <- SumAcrossClusters(mats.filt.lst[[jmark]], cnames.keep.lst = jcnames.keep.lst)
  
  jmat.pbulk.long <- bind_rows(jmat.pbulk) %>%
    as.data.frame()
  jmat.pbulk.long$gene <- names(jmat.pbulk[[1]]) 
  jmat.pbulk.long <- jmat.pbulk.long %>%
    melt()
  colnames(jmat.pbulk.long) <- c("gene", "ctype", "cuts")
  return(jmat.pbulk.long)
})


dat.pbulk.long <- lapply(jmarks, function(jmark){
  dat.pbulk <- dat.merge %>%
    filter(mark == jmark) %>%
    group_by(mark, ctype) %>%
    summarise(spikeincounts = sum(spikeincounts),
              chromocounts = sum(chromocounts))
  jmat.pbulk.long.annot <- left_join(jmat.bulk.long.lst[[jmark]], dat.pbulk)  
}) %>%
  bind_rows()

dat.pbulk.long$mark <- factor(dat.pbulk.long$mark, levels = jmarks)
dat.pbulk.long$ctype <- gsub("^HSCs", "aHSCs", dat.pbulk.long$ctype)

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

ggplot(dat.pbulk.long, aes(x = log2(cuts / spikeincounts), fill = ctype)) + 
  geom_density(alpha = 0.5) + 
  theme_bw() + 
  facet_grid(ctype ~ mark) + 
  scale_fill_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pbulk.long, aes(x = log2(cuts / chromocounts), fill = ctype)) + 
  geom_density(alpha = 0.5) + 
  theme_bw() + 
  facet_grid(ctype ~ mark) + 
  scale_fill_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pbulk.long, aes(x = mark, fill = ctype, y = log2(cuts / spikeincounts))) + 
  geom_boxplot() + 
  theme_bw() + 
  # facet_grid(ctype ~ mark) + 
  scale_fill_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# 2D plot
library(reshape2)
dat.pbulk.wide.spike <- dcast(data = dat.pbulk.long %>%
                          rowwise() %>%
                          mutate(cuts2spike = log2(cuts / spikeincounts),
                                 cuts2chromo = log2(cuts / chromocounts)), 
                        formula = gene ~ ctype + mark, value.var = "cuts2spike")

dat.pbulk.wide.chromo <- dcast(data = dat.pbulk.long %>%
                          rowwise() %>%
                          mutate(cuts2spike = log2(cuts / spikeincounts),
                                 cuts2chromo = log2(cuts / chromocounts)), 
                        formula = gene ~ ctype + mark, value.var = "cuts2chromo")


m.spike <- ggplot(dat.pbulk.wide.spike, aes(x = Granus_H3K4me1 - aHSCs_H3K4me1, y = Granus_H3K27me3 - aHSCs_H3K27me3)) + 
  geom_point(alpha = 0.1) +
  theme_bw() + 
  geom_density_2d() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0, color = 'blue') + geom_hline(yintercept = 0, color = 'blue')


m.chromo <- ggplot(dat.pbulk.wide.chromo, aes(x = Granus_H3K4me1 - aHSCs_H3K4me1, y = Granus_H3K27me3 - aHSCs_H3K27me3)) + 
  geom_point(alpha = 0.1) +
  theme_bw() + 
  geom_density_2d() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0, color = 'blue') + geom_hline(yintercept = 0, color = 'blue')

JFuncs::multiplot(m.spike, m.chromo, cols = 2)
