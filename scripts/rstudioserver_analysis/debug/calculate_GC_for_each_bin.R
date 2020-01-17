# Jake Yeung
# Date of Creation: 2020-01-14
# File: ~/projects/scchic/scripts/rstudioserver_analysis/debug/calculate_GC_for_each_bin.R
# Get GC for eac hbin genome-wide

rm(list=ls())

library(genoset)
library(BSgenome.Mmusculus.UCSC.mm10)

# load LDA 

inf <- "/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_Unenriched_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_Unenriched_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj"
load(inf, v=T)

# coords <- unique(mat.sub.merge$coord)
coords <- unique(rownames(count.mat))
gr.dat <- data.frame(seqnames = sapply(coords, GetChromo), start = sapply(coords, GetStart), end = sapply(coords, GetEnd))
gr <- GenomicRanges::makeGRangesFromDataFrame(gr.dat)
gr.gc <- calcGC(object = gr, bsgenome = BSgenome.Mmusculus.UCSC.mm10)

gr.gc.dat <- data.frame(chromo = sapply(names(gr.gc), GetChromo), 
                        start = as.numeric(sapply(names(gr.gc), GetStart)),
                        end = as.numeric(sapply(names(gr.gc), GetEnd)),
                        gc = gr.gc) %>%
  mutate(midpt = start + (end - start) / 2)

# ggplot(gr.gc.dat %>% filter(midpt > 40e6 & midpt < 60e6), aes(x = midpt / 10^6, y = gc)) + geom_line() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Compare with chromosome  ------------------------------------------------

# save object

save(gr.gc.dat, file = "/home/jyeung/hpc/scChiC/from_rstudioserver/rdata_robjs/gr_gc_dat.RData")
# save(gr.gc.dat, file = "/Users/yeung/data/scchic/robjs/gr_gc_dat.RData")

