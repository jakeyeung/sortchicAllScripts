# Jake Yeung
# Date of Creation: 2019-05-17
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/11-plot_cell_counts_on_umap_all_marks.R
# Plot umap all marks 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/AuxLDA.R")


# Constants ---------------------------------------------------------------

# outdir <- "/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory"
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 1){
  outdir <- args[[1]]
} else {
  outdir <- "/tmp"
}
assertthat::assert_that(dir.exists(outdir))
indir <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/bamlist_for_merging_build95_B6"
assertthat::assert_that(dir.exists(indir))
inf.traj <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/B6_objs/traj_objs_all_marks.Rdata"
assertthat::assert_that(file.exists(inf.traj))

# Load LDA objects --------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmark <- jmarks[["H3K4me1   ]]
infs <- lapply(jmarks, function(jmark) list.files(indir, pattern = paste0(jmark, ".RData"), full.names = TRUE))

tm.result.lst <- lapply(infs, LoadGetTmResult)


# Load trajectories -------------------------------------------------------
# 
load(inf.traj, v=T)


# Load count mats nonbin --------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jbin <- FALSE
kstr <- "30_40_50"
infs.lst <- lapply(jmarks, function(jmark) paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/ldaAnalysisBins_BinCellFilt2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.", jbin, ".no_filt/lda_out_meanfilt.B6_", jmark, "_pcutoff_0.CountThres0.K-", kstr, ".Robj"))
lapply(infs.lst, function(inf) assertthat::assert_that(file.exists(inf)))

# print('/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/ldaAnalysisBins_BinCellFilt2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.FALSE.no_filt/lda_out_meanfilt.B6_H3K27me3_pcutoff_0.CountThres0.K-30_40_50.Robj')
# print(infs.lst[["H3K27me3"]])

kchoose <- "auto"
count.mat.lst <- lapply(jmarks, function(jmark){
  inf <- infs.lst[[jmark]]
  out.objs <- LoadLDABins(jmark, jbin = NA, top.thres = 0.995, inf = inf, convert.chr20.21.to.X.Y = FALSE, add.chr.prefix = TRUE, choose.k = kchoose)
  return(out.objs$count.mat)
})

cell.sums <- lapply(jmarks, function(jmark){
  count.mat <- count.mat.lst[[jmark]]
  cell.sum <- Matrix::colSums(count.mat)
  dat <- data.frame(mark = jmark,
                    cell = names(cell.sum),
                    cell.sum = cell.sum,
                    cell.sum.log10 = log10(cell.sum))
  return(dat)
}) %>%
  bind_rows()

cell.sums <- left_join(cell.sums, dat.umap.long.trajs %>% bind_rows() %>% dplyr::select(cell, umap1, umap2, mark))

jsize <- 4

# plot the PDF 

# x axis is trajectory, y axis is cell sums in log and linear 
# for 3 trajectories: 
jtrajs <- c("granu", "lymphoid", "eryth")

traj.info.all <- lapply(jmarks, function(jmark){
  traj.info <- lapply(jtrajs, function(jtraj){
    return(trajs[[jmark]][[jtraj]] %>% mutate(traj = jtraj, mark = jmark) %>% dplyr::select(c(-umap1, -umap2)))
  }) %>%
    bind_rows()
}) %>% bind_rows()

cell.sums.merged.all <- left_join(cell.sums, traj.info.all)

m.all.linear <- ggplot(cell.sums.merged.all %>% filter(!is.na(traj)), aes(x = lambda, y = cell.sum, color = traj)) + geom_point(alpha = 0.3) + facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
m.all <- ggplot(cell.sums.merged.all %>% filter(!is.na(traj)), aes(x = lambda, y = cell.sum, color = traj)) + geom_point(alpha = 0.3) + facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_log10() 
m.all.log <- ggplot(cell.sums.merged.all %>% filter(!is.na(traj)), aes(x = lambda, y = cell.sum.log10, color = traj)) + geom_point(alpha = 0.3) + facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

pdf(paste0("/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory/cell_counts_over_umap.", Sys.Date(), ".pdf"))
# plot cell sums in linear
for (jmark in jmarks){
  m <- PlotXYWithColor(cell.sums %>% filter(mark == jmark), xvar = "umap1", yvar = "umap2", cname = "cell.sum", jcol = scales::muted("darkred"), jtitle = jmark, jsize = jsize)
  print(m)
}
# plot cell sums in log 
for (jmark in jmarks){
  m.log <- PlotXYWithColor(cell.sums %>% filter(mark == jmark), xvar = "umap1", yvar = "umap2", cname = "cell.sum.log10", jcol = scales::muted("darkred"), jtitle = jmark, jsize = jsize)
  print(m.log)
}
print(m.all.linear)
print(m.all)
print(m.all.log)
for (jmark in jmarks){
  traj.info <- lapply(jtrajs, function(jtraj){
    return(trajs[[jmark]][[jtraj]] %>% mutate(traj = jtraj, mark = jmark) %>% dplyr::select(c(-umap1, -umap2)))
  }) %>%
    bind_rows()
  cell.sums.merged <- left_join(cell.sums, traj.info)
  m.traj <- ggplot(cell.sums.merged %>% filter(!is.na(traj)), aes(x = lambda, y = cell.sum.log10, color = traj)) + facet_wrap(~traj) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark)
  m.traj.lin <- ggplot(cell.sums.merged %>% filter(!is.na(traj)), aes(x = lambda, y = cell.sum, color = traj)) + facet_wrap(~traj) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark)
  print(m.traj)
  print(m.traj.lin)
}
# put marks in one graphs
dev.off()
