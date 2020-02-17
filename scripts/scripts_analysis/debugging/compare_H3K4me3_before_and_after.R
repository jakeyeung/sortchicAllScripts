# Jake Yeung
# Date of Creation: 2019-11-28
# File: ~/projects/scchic/scripts/scripts_analysis/debugging/compare_H3K4me3_before_and_after.R
# Compare raw counts before and after


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

# load old mat
# inf.old <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA/stringent_filter/B6_H3K4me3_pcutoff_0.95_binfilt_cellfilt.2019-06-03.stringent_filter.RData"
inf.old <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA/B6_H3K4me1_pcutoff_0.95_binfilt_cellfilt.2019-05-11.RData"
load(inf.old, v=T)
mat.old <- count.dat$counts
rownames(mat.old) <- paste("chr", rownames(mat.old), sep = "")

# load new mat
# mat.new <- readRDS("/Users/yeung/data/scchic/quality_control_PZ_Bl6-BM-Linneg-Stemcells_AllMerged/B6BM_Unenriched_H3K4me3_pcutoff_0.5_binfilt_cellfilt.2019-11-27.rds")
mat.new <- readRDS("/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-All_allmarks/PZ-Bl6-BM-All_Unenriched.H3K4me1.2019-11-23.rds")  # old antibody
rows.com <- intersect(rownames(mat.new), rownames(mat.old))
cols.com <- intersect(colnames(mat.new), colnames(mat.old))

print(nnzero(mat.new) / length(mat.new))
print(nnzero(mat.old) / length(mat.old))

# compare sample by sampe??
# jcell <- "B6-13W1-BM-H3K4me3-3_175"
jcell <- cols.com[[1]]
jcell <- cols.com[[10]]
jdat <- as.data.frame(cbind(old = mat.old[rows.com, jcell], new = mat.new[rows.com, jcell]))
ggplot(jdat, aes(x = old, y = new)) + geom_jitter() + ggtitle(jcell)

# hierarchical clusterinig???

# check H3K4me1

# try with common rows?
lsi.out.old <- RunLSI(as.matrix(mat.old[rows.com, ]))
lsi.out.new <- RunLSI(as.matrix(mat.new[rows.com, ]))

jsettings <- umap.defaults
jsettings$n_neighbors <- 50
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

umap.out.lsi.old <- umap(lsi.out.old$u, config = jsettings)
dat.umap.long.lsi.old <- data.frame(cell = rownames(umap.out.lsi.old$layout), umap1 = umap.out.lsi.old$layout[, 1], umap2 = umap.out.lsi.old$layout[, 2], stringsAsFactors = FALSE)
ggplot(dat.umap.long.lsi.old, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Old count matrix")

umap.out.lsi.new <- umap(lsi.out.new$u, config = jsettings)
dat.umap.long.lsi.new <- data.frame(cell = rownames(umap.out.lsi.new$layout), umap1 = umap.out.lsi.new$layout[, 1], umap2 = umap.out.lsi.new$layout[, 2], stringsAsFactors = FALSE)
ggplot(dat.umap.long.lsi.new, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("New count matrix")

# merge the two together

# colnames(mat.new) <- paste(colnames(mat.new), "new", sep = ".")
# colnames(mat.old) <- paste(colnames(mat.old), "old", sep = ".")

plot(density(log2(unlist(as.matrix(mat.old[rows.com, ])) + 1)))
plot(density(log2(unlist(as.matrix(mat.new[rows.com, ])) + 1)))

# do global comparison??
x <- unlist(mat.old[rows.com, cols.com])
y <- unlist(mat.new[rows.com, cols.com])

plot(x, y)     
