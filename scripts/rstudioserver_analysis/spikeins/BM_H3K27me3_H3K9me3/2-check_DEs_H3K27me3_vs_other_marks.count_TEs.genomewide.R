# Jake Yeung
# Date of Creation: 2021-01-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/2-check_DEs_H3K27me3_vs_other_marks.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(JFuncs)
library(scchicFuncs)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)


jstart <- Sys.time()

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

fits.out <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.fits <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData")
  # inf.fits <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.allmarks_spikeins/poisson_fit_bins.", jmark, ".2020-12-25.newannot2.witherrors.MoreBins.newestmeta.RData")
  load(inf.fits, v=T)
  params.long <- SummarizeParamsPvalues(jfits.lst, jmark = jmark, paramname = "Cluster") %>%
    mutate(log2fc = estimate / log(2))
  params.long$padj <- p.adjust(params.long$pval.param)
  means.long <- SummarizeMeanValue(jfits.lst, jmark = jmark)
  
  jnames <- names(jfits.lst); names(jnames) <- jnames
  pvals.long <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  
  
  return(list(params.long = params.long, means.long = means.long, pvals.long = pvals.long))
})

# Compare H3K27me3 vs other marks -----------------------------------------


jmark.ref <- "H3K27me3"
jmark.compare <- "H3K4me1"

jmark.ref <- "H3K27me3"
jmark.compare <- "H3K9me3"
params.long.ref <- fits.out[[jmark.ref]]$params.long
params.long2 <- fits.out[[jmark.compare]]$params.long

params.common <- intersect(unique(params.long2$param), unique(params.long.ref$param))

params.long.merge <- left_join(params.long2, params.long.ref %>% filter(param %in% params.common), by = c("bin", "param"))

ggplot(params.long.merge, aes(x = log2fc.x, y = log2fc.y)) + 
  geom_point(alpha = 0.25) + 
  geom_density_2d() + 
  facet_wrap(~param) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# ggplot(params.long.merge %>% filter(padj.x < 10^-1), aes(x = log2fc.x, y = log2fc.y)) +
ggplot(params.long.merge %>% filter(pval.param.x < 10^-10), aes(x = log2fc.x, y = log2fc.y)) +
# ggplot(params.long.merge %>% filter(padj.y < 10^-50), aes(x = log2fc.x, y = log2fc.y)) + 
  geom_point(alpha = 0.25) + 
  geom_density_2d() + 
  facet_wrap(~param) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Claculate TEs in bins  --------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

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


jbins <- params.long.merge$bin

dat.bins <- data.frame(Chr = sapply(jbins, JFuncs::GetChromo), 
                       Start = sapply(jbins, JFuncs::GetStart, returnAsInt = TRUE),
                       End = sapply(jbins, JFuncs::GetEnd, returnAsInt = TRUE), 
                       bin = jbins,
                       stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(Length = End - Start,
         Mid = Start  + Length / 2)


bins.gr <- makeGRangesFromDataFrame(dat.bins, keep.extra.columns = TRUE)
# gr.out <- findOverlaps(bins.gr, tes.gr, type = "within")
# dat.overlap <- data.frame(bins.gr[queryHits(gr.out), ], tes.gr[subjectHits(gr.out), ], stringsAsFactors = FALSE)
gr.out <- findOverlaps(query = tes.gr, subject = bins.gr, type = "within")
dat.overlap <- data.frame(bins.gr[subjectHits(gr.out), ], tes.gr[queryHits(gr.out), ], stringsAsFactors = FALSE)

outf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_gc_analysis/te_counts_50kb_bins.rds")

saveRDS(dat.overlap, file = outf)



