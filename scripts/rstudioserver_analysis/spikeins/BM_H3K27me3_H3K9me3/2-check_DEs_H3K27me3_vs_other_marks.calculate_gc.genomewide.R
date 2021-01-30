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



# Check AT/GC content?  ---------------------------------------------------


library(genoset)
library(BSgenome.Mmusculus.UCSC.mm10)


# bins.keep.k9me3.fg <- subset(params.long.merge, pval.param.x < 10^-10)$bin
# bins.keep.k9me3.bg <- sample(subset(params.long.merge, pval.param.x > 10^-10)$bin, size = length(bins.keep.k9me3.fg), replace = FALSE)
# bins.keep.k27me3.fg <- subset(params.long.merge, pval.param.y < 10^-30)$bin
# bins.keep.k27me3.bg <- sample(subset(params.long.merge, pval.param.y > 10^-30)$bin, size = length(bins.keep.k27me3.fg), replace = FALSE)
# 
# bins.lst <- list(bins.keep.k9me3.fg, bins.keep.k9me3.bg, bins.keep.k27me3.fg, bins.keep.k27me3.bg)
# names(bins.lst) <- c("H3K9me3.fg", "H3K9me3.bg", "H3K27me3.fg", "H3K27me3.bg")



# # system.time(
#   gr.gc.dat.lst <- lapply(bins.lst, function(bins.keep){
#   })
# # )
  
coords <- params.long.merge$bin
print(length(coords))
gr.dat <- data.frame(seqnames = sapply(coords, GetChromo), start = sapply(coords, GetStart), end = sapply(coords, GetEnd), bname = coords)
print(head(gr.dat))
gr <- GenomicRanges::makeGRangesFromDataFrame(gr.dat)
names(gr) <- coords
gr.gc <- calcGC(object = gr, bsgenome = BSgenome.Mmusculus.UCSC.mm10)
gr.gc.dat <- data.frame(bname = names(gr.gc), gc = gr.gc, stringsAsFactors = FALSE)

save(gr.gc.dat, file = "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_gc_analysis/gcs_genomewide.RData")

print(Sys.time() - jstart)

# # Check mutual exclusive --------------------------------------------------
# 
# 
# params.mean.ref <- fits.out[[jmark.ref]]$means.long %>%
#   filter(estimate > -200)
# params.mean2 <- fits.out[[jmark.compare]]$means.long %>%
#   filter(estimate > -200)
# 
# params.mean.merge <- left_join(params.mean.ref, params.mean2, by = "bin")
# 
# ggplot(params.mean.merge, aes(x = exp(estimate.x), y  = exp(estimate.y))) + 
#   geom_point(alpha = 0.5) + 
#   xlab(jmark.ref) + 
#   ylab(jmark.compare) + 
#   geom_density_2d() + 
#   theme_bw() + 
#   ggtitle(paste(jmark.ref, jmark.compare)) + 
#   # coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())







