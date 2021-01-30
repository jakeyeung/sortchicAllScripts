# Jake Yeung
# Date of Creation: 2021-01-26
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/6-overlap_bin_with_H3K4me1_peaks_make_counts.R
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




library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)



# Load fits ---------------------------------------------------------------

jmark <- "H3K9me3"
inf.fits <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData")

load(inf.fits, v=T)

params.long <- SummarizeParamsPvalues(jfits.lst, jmark = jmark, paramname = "Cluster") %>%
  mutate(log2fc = estimate / log(2))
params.long$padj <- p.adjust(params.long$pval.param)
jnames <- names(jfits.lst); names(jnames) <- jnames
pvals.long <- lapply(jnames, function(jname){
  x <- jfits.lst[[jname]]
  xvec <- x$pval
  data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
}) %>%
  bind_rows()

ggplot(params.long %>% filter(abs(estimate) < 5), aes(x = estimate, fill = param)) + 
  geom_density() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  facet_wrap(~param) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Get celltype specific genes  --------------------------------------------

pval.cutoff <- 1e-10

# # qval.cutoff <- 1e-10
# # params.filt <- subset(params.long, padj < qval.cutoff) %>%
# params.filt <- subset(params.long, pval.param < pval.cutoff) %>%
#   arrange(padj)
# print(dim(params.filt))


bins.filt <- subset(pvals.long, pval < pval.cutoff)$bin

ggplot(params.long %>% filter(bin %in% bins.filt), aes(x = estimate, fill = param)) + 
  geom_density() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  facet_wrap(~param) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# print(params.long)
# jbin <- "chr1:12800000-12850000"
# jbin <- bins.filt[[2]]
jbin <- bins.filt[[1]]
print(jbin)

jsub <- subset(params.long, bin == jbin)

print(jsub)



# Load LDA outputs --------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

# for K27me3 
inf.lda.lst <- lapply(jmarks, function(jmark){
  if (jmark != "H3K9me3"){
    # inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/tm_results_from_projs/tm_result_old_to_new.", jmark, ".2020-12-28.RData"))
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins.from_sitecount_mat.from_same_annot_file/lda_outputs.count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.K-30.binarize.FALSE/ldaOut.count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.K-30.Robj"))
  } else {
    # inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/tm_results_from_projs/tm_result_new_to_old.", jmark, ".2020-12-28.RData"))
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins/lda_outputs.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.binarize.FALSE/ldaOut.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.Robj"))
  }
  assertthat::assert_that(file.exists(inf.lda.tmp))
  
  return(inf.lda.tmp)
})

out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.lda <- inf.lda.lst[[jmark]]
  load(inf.lda, v=T)  # out.lda, count.mat
  tm.result <- posterior(out.lda)
  return(list(tm.result = tm.result, count.mat = count.mat))
})

count.mat.lst <- lapply(out.lst, function(out){
  out$count.mat
})

count.mat.renamed.lst <- lapply(count.mat.lst, function(jmat.tmp){
  rownames(jmat.tmp) <- sapply(rownames(jmat.tmp), function(x) strsplit(x, ";")[[1]][[2]]) 
  return(jmat.tmp)
})

# rnames.orig.k4me1 <- rownames(count.mat.lst$H3K4me1)
# rnames.new.k4me1 <- sapply(rnames.orig.k4me1, function(x) strsplit(x, ";")[[1]][[2]])


# Intersect bins with peaks  ----------------------------------------------

dat.bins <- data.frame(Chr = sapply(bins.filt, JFuncs::GetChromo), 
                       Start = sapply(bins.filt, JFuncs::GetStart, returnAsInt = TRUE),
                       End = sapply(bins.filt, JFuncs::GetEnd, returnAsInt = TRUE),
                       Name = bins.filt,
                       stringsAsFactors = FALSE)
gr.bins <- makeGRangesFromDataFrame(dat.bins, keep.extra.columns = TRUE)



jmarks.filt <- c("H3K4me1", "H3K4me3", "H3K27me3")
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.from_peaks.overlap_H3K9me3_50kb_bins"
for (jmark in jmarks.filt){
  # dat.peaks <- data.frame(Chr = sapply(rnames.new.k4me1, JFuncs::GetChromo),
  #                         Start = sapply(rnames.new.k4me1, JFuncs::GetStart, returnAsInt = TRUE),
  #                         End = sapply(rnames.new.k4me1, JFuncs::GetEnd, returnAsInt = TRUE),
  #                         Name = rnames.new.k4me1,
  #                         stringsAsFactors = FALSE)
  rnames.new <- rownames(count.mat.renamed.lst[[jmark]])
  dat.peaks <- data.frame(Chr = sapply(rnames.new, JFuncs::GetChromo),
                          Start = sapply(rnames.new, JFuncs::GetStart, returnAsInt = TRUE),
                          End = sapply(rnames.new, JFuncs::GetEnd, returnAsInt = TRUE),
                          Name = rnames.new,
                          stringsAsFactors = FALSE)
  
  # bins.filt
  gr.peaks <- makeGRangesFromDataFrame(dat.peaks, keep.extra.columns = TRUE)
  gr.out <- findOverlaps(query = gr.bins, subject = gr.peaks, type = "any")
  
  
  # Load peaks, filter out K9me3 bins only  ---------------------------------
  
  dat.overlap <- data.frame(gr.peaks[subjectHits(gr.out), ], gr.bins[queryHits(gr.out), ], stringsAsFactors = FALSE)
 
  
  
  
  # Write new count mat for LDA  --------------------------------------------
  
  jmat <- count.mat.renamed.lst[[jmark]]
  rnames.keep <- rownames(jmat)[rownames(jmat) %in% dat.overlap$Name]
  jmat.filt <- jmat[rnames.keep, ]
  print(jmark)
  print(dim(jmat.filt))
  
  
  # Write output ------------------------------------------------------------
  
  fname <- paste0("count_mat_by_peaks.", jmark, ".overlap_dynamic_k9_bins.rds")
  outrds <- file.path(outdir, fname)
  saveRDS(jmat.filt, file = outrds)
}


