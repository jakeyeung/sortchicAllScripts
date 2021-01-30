# Jake Yeung
# Date of Creation: 2021-01-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/DE_analysis/2-downstream_DE_analysis_TSS.R
# 


rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs) 

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned"
outname <- paste0("TSS_peaks_TEs_analysis_overlaps.", Sys.Date(), ".more.pdf")
outpdf <- file.path(outdir, outname)

make.plots <- TRUE

if (make.plots){
  pdf(file = outpdf, useDingbats = FALSE)
}

# Load .metasmetas  -------------------------------------------------------------

indir.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned"
dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  fname.tmp <- paste0("cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
  inf <- file.path(indir.meta, fname.tmp)
  jdat <- fread(inf)
  if (jmark == "H3K9me3"){
    jdat$mark <- jmark
  }
  return(jdat)
})


# Check DEs ---------------------------------------------------------------

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again"
# fname <- 

# load outputs
outs.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  fname.tmp <- paste0("poisson_fit_TSS.", jmark, ".2020-12-12.newannot2.witherrors.TSS.RData")
  inf <- file.path(indir, fname.tmp)
  load(inf, v=T)
  return(list(jfits.lst = jfits.lst, jmat.mark = jmat.mark))
})



# Check fraction of cuts at TSS vs total  ---------------------------------

# check overlapping??
rnames <- rownames(outs.lst$H3K4me3$jmat.mark)
jgenes <- sapply(rnames, function(x) strsplit(x, "\\.")[[1]][[4]])
rnames.dedup <- rnames[!duplicated(jgenes)]
coords.dedup <- paste("chr", sapply(rnames.dedup, function(x) strsplit(x, split = ";")[[1]][[1]]), sep = "")

cuts.at.tss <- lapply(jmarks, function(jmark){
  print(jmark)
  jout <- outs.lst[[jmark]]
  jmat.mark <- jout$jmat.mark
  # rnames <- rownames(jmat.mark)
  # jgenes <- sapply(rnames, function(x) strsplit(x, "\\.")[[1]][[4]])
  # rnames.dedup <- rnames[!duplicated(jgenes)]
  if (jmark != "H3K27me3"){
    rnames.keep <- rownames(jmat.mark) %in% rnames.dedup
  } else {
    rnames.keep <- rownames(jmat.mark) %in% coords.dedup
  }
  print(dim(jmat.mark))
  jmat.mark <- jmat.mark[rnames.keep, ]
  print(dim(jmat.mark))
  jdat.out <- data.frame(cell = colnames(jmat.mark), ncuts.tss = colSums(jmat.mark), stringsAsFactors = FALSE)
  return(jdat.out)
})

dat.metas.merge <- lapply(jmarks, function(jmark){
  left_join(dat.metas[[jmark]], cuts.at.tss[[jmark]])
}) %>%
  bind_rows() %>%
  mutate(frac = ncuts.tss / cuts_total,
         fracpeaks = ncuts.tss / cuts_in_peak)

ggplot(dat.metas.merge, aes(x = cluster, y = log2(ncuts.tss / cuts_total))) + 
  geom_boxplot() + 
  facet_wrap(~mark, scales = "free_x", nrow = 1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.metas.merge, aes(x = cluster, y = ncuts.tss / cuts_total)) + 
  geom_boxplot() + 
  scale_y_log10() + 
  facet_wrap(~mark, scales = "free_x", nrow = 1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.metas.merge, aes(x = cluster, y = ncuts.tss / cuts_total)) + 
  geom_boxplot() + 
  facet_wrap(~mark, scales = "free_x", nrow = 1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# plot mark by mark, order decreasing to increasing

jylim2 <- c(-8, 0)
jylim <- c(0, 1)

for (jmark in jmarks){
  m <- ggplot(dat.metas.merge %>% filter(mark == jmark), 
              aes(x = forcats::fct_reorder(.f = cluster, .x = ncuts.tss / cuts_total, .fun = median, .desc = FALSE), 
                  # y = log2(ncuts.tss / cuts_total))) + 
                  y = (ncuts.tss) / cuts_total)) + 
    scale_y_log10() + 
    geom_boxplot() + 
    xlab("") + 
    ggtitle(jmark) +
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
  print(m + coord_cartesian(ylim = jylim + 0.001))
}

for (jmark in jmarks){
  m <- ggplot(dat.metas.merge %>% filter(mark == jmark), 
              aes(x = forcats::fct_reorder(.f = cluster, .x = ncuts.tss / cuts_total, .fun = median, .desc = FALSE), 
                  # y = log2(ncuts.tss / cuts_total))) + 
                  y = log2(ncuts.tss / cuts_total))) + 
    geom_boxplot() + 
    xlab("") + 
    ggtitle(jmark) +
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
  print(m + coord_cartesian(ylim = jylim2))
}

for (jmark in jmarks){
  m <- ggplot(dat.metas.merge %>% filter(mark == jmark), 
              aes(x = forcats::fct_reorder(.f = cluster, .x = ncuts.tss / cuts_total, .fun = median, .desc = FALSE), 
                  # y = log2(ncuts.tss / cuts_total))) + 
                  y = log2(ncuts.tss / cuts_total))) + 
    geom_boxplot() + 
    xlab("") + 
    ggtitle(jmark) +
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
  print(m + coord_cartesian(ylim = jylim2))
}


ggplot(dat.metas.merge, aes(x = cluster, y = log2(ncuts.tss / cuts_in_peak))) + 
  geom_boxplot() + 
  facet_wrap(~mark, scales = "free_x", nrow = 1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.metas.merge, aes(x = cluster, y = ncuts.tss / cuts_in_peak)) + 
  geom_boxplot() + 
  facet_wrap(~mark, scales = "free_x", nrow = 1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

for (jmark in jmarks){
  m <- ggplot(dat.metas.merge %>% filter(mark == jmark), 
              aes(x = forcats::fct_reorder(.f = cluster, .x = ncuts.tss / cuts_in_peak, .fun = median, .desc = FALSE), 
                  # y = log2(ncuts.tss / cuts_total))) + 
                  y = ncuts.tss / cuts_in_peak)) + 
    scale_y_log10() + 
    geom_boxplot() + 
    xlab("") + 
    ggtitle(jmark) +
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
  print(m + coord_cartesian(ylim = jylim + 0.001))
}

for (jmark in jmarks){
  m <- ggplot(dat.metas.merge %>% filter(mark == jmark), 
              aes(x = forcats::fct_reorder(.f = cluster, .x = ncuts.tss / cuts_in_peak, .fun = median, .desc = FALSE), 
                  # y = log2(ncuts.tss / cuts_total))) + 
                  y = ncuts.tss / cuts_in_peak)) + 
    geom_boxplot() + 
    xlab("") + 
    ggtitle(jmark) +
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
  print(m + coord_cartesian(ylim = jylim))
}



for (jmark in jmarks){
  m <- ggplot(dat.metas.merge %>% filter(mark == jmark), 
              aes(x = forcats::fct_reorder(.f = cluster, .x = ncuts.tss / cuts_in_peak, .fun = median, .desc = FALSE), 
                  y = log2(ncuts.tss / cuts_in_peak))) + 
    geom_boxplot() + 
    xlab("") + 
    ggtitle(jmark) +
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
  print(m + coord_cartesian(ylim = jylim2))
}

  
ggplot(dat.metas.merge, aes(x = cluster, y = log2(ncuts.tss / spikein_cuts))) + 
  geom_boxplot() + 
  facet_wrap(~mark, scales = "free_x", nrow = 1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
 

# Load K27me3 peaks  ------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

indir.peaks <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/hiddendomains_outputs/hiddendomains_outputs_minlength_500.mincount_-10.FromR.maxcount_10_40_60")
dir.exists(indir.peaks)


ctypes <- unique(dat.metas$H3K27me3$cluster)
names(ctypes) <- ctypes

ctype <- ctypes[[1]]


dat.peaks <- lapply(ctypes, function(ctype){
  jdir <- list.files(indir.peaks, pattern = paste0("-", ctype, "-"), full.names = TRUE)
  jfile <- list.files(jdir, pattern = "cutoff_analysis.bed", full.names = TRUE)
  print(jfile)
  jdat.bed <- fread(jfile)
  colnames(jdat.bed) <- c("Chr", "Start", "End", "Name")
  jdat.bed <- jdat.bed %>%
    rowwise() %>%
    mutate(Length = End - Start,
           Mid = Start  + Length / 2)
  jdat.bed$ctype <- ctype
  return(jdat.bed)
})


# Compare total lengths ---------------------------------------------------

dat.peaks.sum <- dat.peaks %>%
  bind_rows() %>%
  group_by(ctype) %>%
  summarise(LengthTotal = sum(Length))

ggplot(dat.peaks.sum, aes(x = ctype, y = log2(LengthTotal))) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.peaks.sum, aes(x = ctype, y = LengthTotal)) + 
  geom_col() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.peaks %>% bind_rows() %>% filter(Length > 1000), aes(x = Length, fill = ctype)) + 
  geom_density() + 
  scale_x_log10() + 
  theme_bw() + 
  # coord_cartesian(xlim = c(1000, 10^6)) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~ctype)


# Assign peaks to nearest gene  -------------------------------------------

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
jinf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.100000.bed"

jctype <- "HSPCs"
jctype <- "Eryths"
jctype <- "Granulocytes"

regions.annotated.lst <- lapply(ctypes, function(jctype){
  print(jctype)
  jpeaks <- paste(paste("chr", dat.peaks[[jctype]]$Chr, sep = ""), paste(dat.peaks[[jctype]]$Mid - 1, dat.peaks[[jctype]]$Mid + 1, sep = "-"), sep = ":")
  print(head(jpeaks))
  bins.annot <- AnnotateCoordsFromList(coords.vec = jpeaks, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)
  return(bins.annot$regions.annotated)
})

regions.annotated.lst2 <- lapply(ctypes, function(ctype){
  regions.annotated.lst[[ctype]]$ctype <- ctype
  return(regions.annotated.lst[[ctype]])
})

ggplot(regions.annotated.lst2 %>% bind_rows(), aes(x = abs(distanceToTSS), fill = ctype)) + 
  geom_density(alpha = 0.5) + 
  scale_x_log10() + 
  theme_bw() + 
  facet_wrap(~ctype, ncol = 1) + 
  ggtitle(jctype) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(regions.annotated.lst2 %>% bind_rows(), aes(x = ctype, y = abs(distanceToTSS))) + 
  geom_boxplot(alpha = 0.5) + 
  scale_y_log10() + 
  theme_bw() + 
  ggtitle("Dist to TSS center of peak") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Overlap with TEs --------------------------------------------------------

jval.cutoff <- 2000
# read TEs
inf.tes <- file.path(hubprefix, "jyeung/data/databases/repeatmaskers/mm10/mm10_repeatmasker_track_viz.gz")
dat.tes <- fread(inf.tes, header = FALSE, col.names = c("Chr", "Start", "End", "Name", "Value", "Strand"))

dat.tes$Class <- sapply(dat.tes$Name, function(x) strsplit(x, split = "#")[[1]][[2]])

dat.tes.filt <- subset(dat.tes, !grepl("\\?", Class)) %>%
  rowwise() %>%
  mutate(Class2 = strsplit(Class, "\\/")[[1]][[1]]) %>%
  ungroup() %>%
  filter(Class2 %in% c("LINE", "LTR", "SINE", "DNA"),
         Value > jval.cutoff)

ggplot(dat.tes.filt, aes(x = Value, fill = Class2)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + 
  geom_vline(xintercept = jval.cutoff) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

tes.gr <- makeGRangesFromDataFrame(dat.tes.filt, keep.extra.columns = TRUE)

jctype <- "HSPCs"
dat.overlap.lst <- lapply(ctypes, function(jctype){
  print(jctype)
  jpeaks <- paste(paste("chr", dat.peaks[[jctype]]$Chr, sep = ""), paste(dat.peaks[[jctype]]$Mid - 1, dat.peaks[[jctype]]$Mid + 1, sep = "-"), sep = ":")
  peaks.gr <- makeGRangesFromDataFrame(dat.peaks[[jctype]] %>% mutate(Chr = paste("chr", Chr, sep = "")), keep.extra.columns = TRUE)
  gr.out <- findOverlaps(peaks.gr, tes.gr, type = "within")
  dat.overlap <- data.frame(peaks.gr[queryHits(gr.out), ], tes.gr[subjectHits(gr.out), ], stringsAsFactors = FALSE)
})

dat.overlap.sum.lst <- lapply(dat.overlap.lst, function(dat.overlap){
  dat.overlap %>%
    group_by(seqnames, start, end, Length, Name, Class2) %>%
    summarise(ntes = length(Name.1)) %>%
    mutate(ctype = jctype)
})

dat.overlap.sum2.long <- lapply(ctypes, function(jctype){
  jdat <- dat.overlap.sum.lst[[jctype]]
  dat.overlap.sum2 <- jdat %>%
    group_by(Class2, ctype) %>%
    summarise(ntes = sum(ntes), 
              Length = sum(Length))
}) %>%
  bind_rows()

ggplot(dat.overlap.sum2.long, aes(x = ctype, y = ntes / Length, fill = Class2)) + 
  geom_col(position = "dodge") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.overlap.sum2.long %>% filter(Class2 != "LINE"), aes(x = ctype, y = ntes, fill = Class2)) + 
  geom_col(position = "dodge") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.overlap.sum2.long, aes(x = ctype, y = Length, fill = Class2)) + 
  geom_col() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


if (make.plots){
  dev.off()
}