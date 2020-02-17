# Jake Yeung
# Date of Creation: 2019-11-29
# File: ~/projects/scchic/scripts/scripts_analysis/debugging/check_LDA_plot_raw_data.R
# H3K4me3 looks weird???

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

# Load LDA  ---------------------------------------------------------------

# load old LDA
inf.dat.stringent <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"
assertthat::assert_that(file.exists(inf.dat.stringent))
load(inf.dat.stringent, v=T)
dat.umap.long.stringent <- LoadUmap(inf.dat.stringent)

# load new LDA 
inf.new <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All.H3K4me3_only_first_try/lda_outputs.PZ-Bl6-BM-All_Unenriched.H3K4me3.2019-11-22.K-50.binarize.FALSE/ldaOut.PZ-Bl6-BM-All_Unenriched.H3K4me3.2019-11-22.K-50.Robj"
load(inf.new, v=T)

# get umap
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(posterior(out.lda)$topics, config = jsettings)
dat.umap.long.new <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])


# Load raw  ---------------------------------------------------------------

# load old mat
inf.old <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA/stringent_filter/B6_H3K4me3_pcutoff_0.95_binfilt_cellfilt.2019-06-03.stringent_filter.RData"
load(inf.old, v=T)
mat.old <- count.dat$counts
rownames(mat.old) <- paste("chr", rownames(mat.old), sep = "")

# load new mat
mat.new <- readRDS("/Users/yeung/data/scchic/quality_control_PZ_Bl6-BM-Linneg-Stemcells_AllMerged/B6BM_Unenriched_H3K4me3_pcutoff_0.5_binfilt_cellfilt.2019-11-27.rds")
rows.com <- intersect(rownames(mat.new), rownames(mat.old))
cols.com <- intersect(colnames(mat.new), colnames(mat.old))

print(nnzero(mat.new) / length(mat.new))
print(nnzero(mat.old) / length(mat.old))

mat.old.com <- mat.old[, cols.com]
mat.new.com <- mat.new[, cols.com]

# Plot side by side  ------------------------------------------------------


# add raw data to the umap

# S100a7
jgene <- "S100"
jgene <- "^Sox6$"
jgene <- "^Ccl5$"
jgene <- "^Ccl2$"
(jterm <- subset(out.objs$regions.annotated, grepl(jgene, SYMBOL) & distanceToTSS == 0)$region_coord[[1]])

jvec.old <- data.frame(cell = colnames(mat.old.com), exprs = mat.old.com[jterm, ] / colSums(mat.old.com[, ]), stringsAsFactors = FALSE)
jvec.new <- data.frame(cell = colnames(mat.new.com), exprs = mat.new.com[jterm, ] / colSums(mat.new.com[, ]), stringsAsFactors = FALSE)

oldumap.oldexprs <- left_join(dat.umap.long.stringent %>% filter(cell %in% cols.com), jvec.old)
oldumap.newexprs <- left_join(dat.umap.long.stringent %>% filter(cell %in% cols.com), jvec.new)

newumap.oldexprs <- left_join(dat.umap.long.new %>% filter(cell %in% cols.com), jvec.old)
newumap.newexprs <- left_join(dat.umap.long.new %>% filter(cell %in% cols.com), jvec.new)

m.old.exprs.old <- PlotXYWithColor(oldumap.oldexprs, xvar = "umap1", yvar = "umap2", cname = "exprs", jsize = 3)
print(m.old.exprs.old)

m.old.exprs.new <- PlotXYWithColor(oldumap.newexprs, xvar = "umap1", yvar = "umap2", cname = "exprs", jsize = 3)
print(m.old.exprs.new)

m.old.exprs.old <- PlotXYWithColor(newumap.oldexprs, xvar = "umap1", yvar = "umap2", cname = "exprs", jsize = 3, cont.color = TRUE)
print(m.old.exprs.old)

# m.new.exprs.old <- ggplot(newumap.oldexprs, aes(x = umap1, y = umap2, color = log10(exprs))) + 
#   geom_point(size = 3) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   scale_color_viridis_c()
# print(m.new.exprs.old)

m.old.exprs.new <- PlotXYWithColor(newumap.newexprs, xvar = "umap1", yvar = "umap2", cname = "exprs", jsize = 3)
print(m.old.exprs.new)

