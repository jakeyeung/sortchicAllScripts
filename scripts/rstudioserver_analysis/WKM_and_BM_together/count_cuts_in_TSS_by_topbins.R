# Jake Yeung
# Date of Creation: 2020-06-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/count_cuts_in_TSS_by_topbins.R
# 


rm(list=ls())

jstart <- Sys.time()

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(hash)
library(JFuncs)
library(scchicFuncs)



jdist <- c(10000L)
jnorm <- c("ByTotalFromBins")


# jdist <- 10000
# jnorm <- "ByHetero"
# # jnorm <- "ByTotalFromBins"

jfactor <- jdist / 10000

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/tss_hetero_totalcuts.countsByCluster"
# outpdf <- file.path(outdir, paste0("tss_hetero_and_totalcuts.pdf"))
# outcelltable <- file.path(outdir, paste0("tss_hetero_and_totalcuts.celltable.txt"))

jsuffix <- "MouseBMFromTopics.2000"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bins_tss_totalcuts.countsByCluster"
outpdf <- file.path(outdir, paste0(jsuffix, ".bins_tss_hetero_and_totalcuts.pdf"))
outcelltable <- file.path(outdir, paste0(jsuffix, ".bins_tss_hetero_and_totalcuts.celltable.txt"))
outprefix <- file.path(outdir, paste0(jsuffix, ".bins_tss_hetero_and_totalcuts"))

pdf(outpdf, useDingbats = FALSE)


fewer.k27me3 <- TRUE
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression"
dir.create(indir, showWarnings = TRUE)
jprefix <- file.path(indir, paste0("integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_", fewer.k27me3, ".forPoissonRegression.CountR1only.", "2020-06-05"))


infrdata <- paste0(jprefix, ".smaller.RData")

inf.offsets <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/offsets_hetero_and_totalcuts/MouseBM_HeteroTotalCounts_50kb_bins.2020-06-16.RData"

assertthat::assert_that(file.exists(infrdata))

load(inf.offsets, v=T)
load(infrdata, v=T)

cells.keep <- lapply(tss.mats.filt.fromref.cellfilt, function(jmat){
  colnames(jmat)
})

dat.annots.filt.forfit <- lapply(dat.annots.filt, function(jdat){
  jdat <- subset(jdat, select = c(cell, cluster.new)) %>%
    mutate(Cluster = ifelse(cluster.new == "HSPCs", "aHSPCs", as.character(cluster.new)))  # set HSPC as intercept
  return(jdat)
})

# Read in topbins ---------------------------------------------------------

jclsts <- c("HSPCs", "Granu", "Eryth", "Bcell"); names(jclsts) <- jclsts
# jclst <- "HSPCs"
jannot.main <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations")
jannot.indir <- file.path(jannot.main, jsuffix)

jannot.bed.long <- lapply(jclsts, function(jclst){
  fname <- paste0("MouseBM_TSS_FromTopics.", jclst, ".bsize_2.bed")
  jannot.inf <- file.path(jannot.indir, fname)
  jannot.bed <- fread(jannot.inf, col.names = c("chromo", "Start", "End", "gene")) %>%
    rowwise() %>%
    mutate(mdpt.from2bp = (Start + End) / 2,
           clst = jclst)
}) %>%
  bind_rows()

genes.inbins <- unique(jannot.bed.long$gene)

rnames.common <- Reduce(intersect, lapply(tss.mats.filt.fromref.cellfilt, function(x) rownames(x)))
coords.common <- sapply(rnames.common, function(x) strsplit(x, ";")[[1]][[1]])
genes.common <- sapply(rnames.common, function(x) strsplit(x, ";")[[1]][[2]])

rnames.bed <- data.frame(chromo = sapply(coords.common, GetChromo), 
                         Start = as.integer(sapply(coords.common, GetStart)), 
                         End = as.integer(sapply(coords.common, GetEnd)), 
                         gene = genes.common, 
                         rname = rnames.common, 
                         stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(mdpt.from10kb = (Start + End) / 2) %>%
  left_join(., subset(jannot.bed.long, select = c(gene, mdpt.from2bp, clst), by = "gene"))

# check
ggplot(rnames.bed %>% filter(!is.na(clst)), aes(x = mdpt.from2bp, y = mdpt.from10kb)) + geom_point() + scale_x_log10() + scale_y_log10()



# Filter out  -------------------------------------------------------------

rnames.bed.filt <- subset(rnames.bed, !is.na(clst))
rnames.keep.withdupes <- rnames.bed.filt$rname

tss.counts.binfilt.long <- lapply(jmarks, function(jmark){
  jmat <- tss.mats.filt.fromref.cellfilt[[jmark]]
  jmat.filt <- jmat[rnames.keep.withdupes, ]
  tss.counts <- data.frame(cell = colnames(jmat.filt), ncuts.inbins = colSums(jmat.filt), ncuts.alltss = colSums(jmat), stringsAsFactors = FALSE) %>%
    left_join(., dat.annots.filt.forfit[[jmark]]) %>%
    left_join(., dat.ncuts.hetero.total[[jmark]]) %>%
    rowwise() %>%
    mutate(mark = jmark)
}) %>%
  bind_rows()

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(tss.counts.binfilt.long, aes(x = ncuts.inbins, fill = cluster.new)) + geom_density(alpha = 0.25) + 
  facet_wrap(~mark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_log10() + 
  scale_fill_manual(values = cbPalette)

ggplot(tss.counts.binfilt.long, aes(x = ncuts.alltss, fill = cluster.new)) + geom_density(alpha = 0.25) + 
  facet_wrap(~mark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_log10() + 
  scale_fill_manual(values = cbPalette)

ggplot(tss.counts.binfilt.long, aes(x = ncuts.inbins / ncuts.alltss, fill = cluster.new)) + geom_density(alpha = 0.25) + 
  facet_wrap(~mark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_log10() + 
  scale_fill_manual(values = cbPalette)

ggplot(tss.counts.binfilt.long, aes(x = ncuts.inbins / ncuts.total, fill = cluster.new)) + geom_density(alpha = 0.25) + 
  facet_wrap(~mark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_log10() + 
  scale_fill_manual(values = cbPalette)

ggplot(tss.counts.binfilt.long, aes(x = ncuts.inbins / ncuts.total, fill = cluster.new)) + geom_density(alpha = 0.25) + 
  facet_wrap(~mark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = cbPalette)

ggplot(tss.counts.binfilt.long, aes(x = ncuts.alltss / ncuts.total, fill = cluster.new)) + geom_density(alpha = 0.25) + 
  facet_wrap(~mark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = cbPalette)

ggplot(tss.counts.binfilt.long, aes(x = ncuts.inbins / ncuts.alltss, fill = cluster.new)) + geom_density(alpha = 0.25) + 
  facet_wrap(~mark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = cbPalette)

ggplot(tss.counts.binfilt.long, aes(x = ncuts.inbins / ncuts.hetero, fill = cluster.new)) + geom_density(alpha = 0.25) + 
  facet_wrap(~mark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_log10() + 
  scale_fill_manual(values = cbPalette)


# output cnames:
# cluster.new     ncuts.total     ncuts.hetero    tss.cuts        ncells  mark
tss.counts.pbulk.binfilt.long <- tss.counts.binfilt.long %>%
  group_by(cluster.new, mark) %>%
  summarise(ncuts.inbins = sum(ncuts.inbins),
            ncuts.alltss = sum(ncuts.alltss),
            ncuts.total = sum(ncuts.total),
            ncuts.hetero = sum(ncuts.hetero),
            ncells = length(cell)) %>%
  dplyr::select(c(cluster.new, ncuts.total, ncuts.hetero, ncuts.alltss, ncuts.inbins, ncells, mark))  # add mark at end as padding for last column 

ggplot(tss.counts.pbulk.binfilt.long, aes(x = cluster.new, y = ncuts.inbins, fill = cluster.new)) + geom_col() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~mark)

ggplot(tss.counts.pbulk.binfilt.long, aes(x = cluster.new, y = ncuts.inbins / ncuts.alltss, fill = cluster.new)) + geom_col() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~mark) 

ggplot(tss.counts.pbulk.binfilt.long, aes(x = cluster.new, y = log2(ncuts.inbins / ncuts.alltss), fill = cluster.new)) + geom_col() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~mark) 


dev.off()

# Write to output ---------------------------------------------------------

jclsts.vec <- unique(tss.counts.binfilt.long$cluster.new); names(jclsts.vec) <- jclsts.vec
jmarks.vec <- unique(tss.counts.binfilt.long$mark); names(jmarks.vec) <- jmarks.vec

lapply(jclsts.vec, function(jclst){
  lapply(jmarks.vec, function(jmark){
    outf.tmp1 <- paste0(outprefix, ".", jclst, ".", jmark, ".txt")
    outf.tmp2 <- paste0(outprefix, ".", jclst, ".", jmark, ".NoCname.txt")
    jout <- subset(tss.counts.pbulk.binfilt.long, cluster.new == jclst & mark == jmark)
    fwrite(jout, file = outf.tmp1, quote = FALSE, sep = "\t", col.names = TRUE)
    fwrite(jout, file = outf.tmp2, quote = FALSE, sep = "\t", col.names = FALSE)
  })
})


# # write cell table
fwrite(tss.counts.binfilt.long, file = outcelltable, quote = FALSE, sep = "\t", col.names = TRUE)

