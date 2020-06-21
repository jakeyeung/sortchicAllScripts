# Jake Yeung
# Date of Creation: 2020-06-20
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/count_cuts_in_TSS.R
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

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/tss_hetero_totalcuts.countsByCluster"

outpdf <- file.path(outdir, paste0("tss_hetero_and_totalcuts.pdf"))
outcelltable <- file.path(outdir, paste0("tss_hetero_and_totalcuts.celltable.txt"))


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

# take transcripts: can be different window siezs
refmark <- "H3K4me3"
rnames <- rownames(tss.mats.filt.fromref.cellfilt[[refmark]])
tx.keep <- sapply(rnames, function(rname) paste(strsplit(rname, split = ";")[[1]][2:3], collapse = ";"), USE.NAMES = FALSE)


rm(tss.mats.filt.fromref.cellfilt)  # load a new one later


dat.annots.filt.forfit <- lapply(dat.annots.filt, function(jdat){
  jdat <- subset(jdat, select = c(cell, cluster.new)) %>%
    mutate(Cluster = ifelse(cluster.new == "HSPCs", "aHSPCs", as.character(cluster.new)))  # set HSPC as intercept
  return(jdat)
})

cells.clstrfilt <- lapply(dat.annots.filt.forfit, function(jdat){
  return(jdat$cell)
})


if (jnorm == "ByHetero"){
  ncuts.for.fit <- lapply(dat.ncuts.hetero.total, function(jdat){
    # expects ncuts.total
    subset(jdat, select = c(cell, ncuts.hetero)) %>%
      dplyr::rename(ncuts.total = ncuts.hetero)
  })
} else if (jnorm == "ByTotalFromBins"){
  ncuts.for.fit <- lapply(dat.ncuts.hetero.total, function(jdat){
    # expects ncuts.total
    subset(jdat, select = c(cell, ncuts.total))
  })
} else {
  warning("jnorm must be ByHetero or ByTotalFromBins", jnorm)
}

# adjust by jfactor
ncuts.for.fit <- lapply(ncuts.for.fit, function(jdat){
  jdat$ncuts.total * jfactor  # multiply by factor so different winsizes are comparable
  return(jdat)
})



# Load the new tables ------------------------------------------------

# must be TSS tables
indir.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/debug/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS.all_tx.noR2.Buys"
tss.mats.filt.fromref.cellfilt <- lapply(jmarks, function(jmark){
  print(jmark)
  # fname <- paste0(jmark, ".imputevarfilt.lessstringent.mapq_40.remerged.countTable.TSS.csv")
  fname <- paste0(jmark, ".countTableTSS.mapq_40.TSS_", jdist, ".blfiltered.csv")
  inf <- file.path(indir.tss, fname)
  mat <- ReadMatTSSFormat(inf = inf, as.sparse = TRUE, sort.rnames = TRUE)
  # filter rownames to match tx.keep
  tx.all <- sapply(rownames(mat), function(rname) paste(strsplit(rname, ";")[[1]][2:3], collapse = ";"))
  tx.filt <- tx.all %in% tx.keep
  return(mat[tx.filt, cells.clstrfilt[[jmark]]])
})


# Calculate total TSS counts ----------------------------------------------

tss.sum.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jmat <- tss.mats.filt.fromref.cellfilt[[jmark]]
  tss.dat <- data.frame(cell = colnames(jmat), tss.cuts = colSums(jmat), stringsAsFactors = FALSE) %>%
    left_join(dat.ncuts.hetero.total[[jmark]]) %>%
    left_join(dat.annots.filt.forfit[[jmark]]) %>%
    rowwise() %>%
    mutate(mark = jmark)
}) %>%
  bind_rows()


pdf(outpdf, useDingbats = FALSE)

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(tss.sum.lst, aes(x = tss.cuts / ncuts.total, fill = Cluster)) + geom_density(alpha = 0.25) + facet_wrap(~mark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = cbPalette)

ggplot(tss.sum.lst, aes(x = tss.cuts / ncuts.hetero, fill = Cluster)) + geom_density(alpha = 0.25) + facet_wrap(~mark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = cbPalette)

ggplot(tss.sum.lst, aes(x = tss.cuts / ncuts.hetero, fill = Cluster)) + geom_density(alpha = 0.25) + facet_wrap(~mark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = cbPalette)

ggplot(tss.sum.lst, aes(x = tss.cuts, fill = Cluster)) + geom_density(alpha = 0.25) + facet_wrap(~mark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = cbPalette) + scale_x_log10() 

# Calculate frac of TSS  --------------------------------------------------

# write table TSS, TSS/total, and TSS/hetero H3K4me1, H3K4me3, and H3K27me3 

tss.sum.pbulk.lst <- tss.sum.lst %>%
  group_by(cluster.new, mark) %>%
  summarise(tss.cuts = sum(tss.cuts),
            ncuts.total = sum(ncuts.total), 
            ncuts.hetero = sum(ncuts.hetero),
            ncells = length(cell)) %>%
  dplyr::select(c(cluster.new, ncuts.total, ncuts.hetero, tss.cuts, ncells, mark))
  
ggplot(tss.sum.pbulk.lst, aes(x = cluster.new, y = tss.cuts)) + geom_col() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~mark)

ggplot(tss.sum.pbulk.lst, aes(x = cluster.new, y = tss.cuts / ncells)) + geom_col() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~mark)
  
ggplot(tss.sum.pbulk.lst, aes(x = cluster.new, y = ncuts.total)) + geom_col() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~mark)

ggplot(tss.sum.pbulk.lst, aes(x = cluster.new, y = ncuts.total / ncells)) + geom_col() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~mark)

ggplot(tss.sum.pbulk.lst, aes(x = cluster.new, y = ncuts.hetero / ncells)) + geom_col() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~mark)

dev.off()

# write to output

jclsts <- unique(tss.sum.pbulk.lst$cluster.new)
jmarks <- unique(tss.sum.pbulk.lst$mark)


for (jclst in jclsts){
  for (jmark in jmarks){
    print(paste(jclst, jmark))
    fname1 <- paste0("tss_hetero_and_totalcuts.", jclst, ".", jmark, ".txt")
    fname2 <- paste0("tss_hetero_and_totalcuts.", jclst, ".", jmark, ".NoCname.txt")
    outf1 <- file.path(outdir, fname1)
    outf2 <- file.path(outdir, fname2)
    fwrite(subset(tss.sum.pbulk.lst, cluster.new == jclst & mark == jmark), file = outf1, quote = FALSE, sep = "\t", na = "NA", col.names = TRUE)
    fwrite(subset(tss.sum.pbulk.lst, cluster.new == jclst & mark == jmark), file = outf2, quote = FALSE, sep = "\t", na = "NA", col.names = FALSE)
  }
}

fwrite(tss.sum.lst, file = outcelltable, quote = FALSE, sep = "\t", na = "NA", col.names = TRUE)


