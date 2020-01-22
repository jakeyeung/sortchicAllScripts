# Jake Yeung
# Date of Creation: 2019-11-29
# File: ~/projects/scchic/scripts/scripts_analysis/debugging/compare_H3K4me3_before_and_after_retag.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(irlba)
library(hash)
library(igraph)
library(umap)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmarks <- jmarks[[2]]

jsettings <- umap.defaults
jsettings$n_neighbors <- 30

jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# load old mat
set.seed(0)

# find this bin
jbin.viewer <- "chr7:115,927,442-115,930,186"
jbin.bin <- "chr7:115920000-116020000"
jcell.tmp <- "B6-13W1-BM-H3K4me3-2_155"

# rnames <- rownames(count.dat$counts)
# which(sapply(rnames, function(x) startsWith(x, "7:115[0-9][0-9][0-9][0-9][0-9][0-9]")))
# which(sapply(rnames, function(x) startsWith(x, "7:1159")))

for (jmark in jmarks)({
  
  print(jmark)
  
  if (jmark == "H3K4me3"){
    inf.old <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA/stringent_filter/B6_H3K4me3_pcutoff_0.95_binfilt_cellfilt.2019-06-03.stringent_filter.RData"
  } else {
    inf.old <- paste0("/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA/B6_", jmark, "_pcutoff_0.95_binfilt_cellfilt.2019-05-11.RData")
  }
  load(inf.old, v=T)
  mat.old <- count.dat$counts
  rownames(mat.old) <- paste("chr", rownames(mat.old), sep = "")
  
  inf.new <- paste0("/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-All_allmarks/PZ-Bl6-BM-All_Unenriched.", jmark, ".2019-11-23.rds")
  mat.new <- readRDS(inf.new)
  # inf.new.raw <- paste0("/Users/yeung/data/scchic/from_cluster/count_mats_B6BM_and_LinNeg_StemCells.retag/countTables/", jmark, "-BM_SC-merged.tagged.retagged.countTable.demuxbugfixed.csv")
  # mat.new <- ReadMatSlideWinFormat(inf.new.raw)
  
  rows.com <- intersect(rownames(mat.new), rownames(mat.old))
  cols.com <- intersect(colnames(mat.new), colnames(mat.old))
  
  cells.all <- cols.com
  names(cells.all) <- cols.com
  
  # check a cell in a bin
  # jtmp <- mat.old[jbin.bin, jcell.tmp]  # 3 reads in chr7:115920000-116020000 in B6-13W1-BM-H3K4me3-2_155, IGV verified
  # jtmp.new <- mat.new[jbin.bin, jcell.tmp]  # 4 reads here 
  
  print(nnzero(mat.old) / length(mat.old))
  print(nnzero(mat.new) / length(mat.new))
  
  choose.i <- sample(seq(length(cols.com)), size = 1)
  jcell <- cols.com[[choose.i]]
  jdat <- as.data.frame(cbind(old = mat.old[rows.com, jcell], new = mat.new[rows.com, jcell]))
  m <- ggplot(jdat, aes(x = old, y = new)) + geom_jitter() + ggtitle(paste(jmark, jcell)) 
  print(m)
  
  # calculate correlation for each cell
  jtmp.old <- mat.old[rows.com, cells.all]
  jtmp.new <- mat.new[rows.com, cells.all]
  jcells.all.i <- seq(length(cells.all))
  names(jcells.all.i) <- cells.all
  jcorr.lst <- lapply(jcells.all.i, function(i){
    jcorr <- cor(jtmp.old[, i], jtmp.new[, i])
  })
  # hierarchical clusterinig???
  
  mat.norm.old <- NormalizeMatrix(as.matrix(mat.old[rows.com, ]))
  mat.norm.new <- NormalizeMatrix(as.matrix(mat.new[rows.com, ]))
  jcorr.norm.lst <- lapply(jcells.all.i, function(i){
    jcorr <- cor(mat.norm.old[, i], mat.norm.new[, i])
  })
  hist(unlist(jcorr.norm.lst))
  
  # check sums
  sums.old <- colSums(jtmp.old) / 5
  sums.new <- colSums(jtmp.new) / 5
  
  # check H3K4me1
  # try with common rows?
  lsi.out.old <- RunLSI(as.matrix(jtmp.old))
  lsi.out.new <- RunLSI(as.matrix(jtmp.new))
  
  # which sells are bad?
  # bad.cells <- colnames(mat.new)[!colnames(mat.new) %in% cols.com]
  # jtmp.new.bad <- mat.new[rows.com, bad.cells]
  # 
  # plot(density(sums.new[names(sums.new) %in% bad.cells]))
  # plot(density(sums.new))
  
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 50
  jsettings$min_dist <- 0.1
  jsettings$random_state <- 123
  
  umap.out.lsi.old <- umap(lsi.out.old$u, config = jsettings)
  dat.umap.long.lsi.old <- data.frame(cell = rownames(umap.out.lsi.old$layout), umap1 = umap.out.lsi.old$layout[, 1], umap2 = umap.out.lsi.old$layout[, 2], stringsAsFactors = FALSE)
  m.umap.old <- ggplot(dat.umap.long.lsi.old, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle("LSI on Old count matrix", jmark)
  
  umap.out.lsi.new <- umap(lsi.out.new$u, config = jsettings)
  dat.umap.long.lsi.new <- data.frame(cell = rownames(umap.out.lsi.new$layout), umap1 = umap.out.lsi.new$layout[, 1], umap2 = umap.out.lsi.new$layout[, 2], stringsAsFactors = FALSE)
  m.umap.new <- ggplot(dat.umap.long.lsi.new, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle("LSI on New count matrix", jmark)
  
  
  pdf(paste0("/Users/yeung/data/scchic/pdfs/compare_old_new_pipeline/compare_pipelines.retag.", jmark, ".samecols.nodedup.pdf"))
    plot(density(unlist(jcorr.lst)))
    hist(unlist(jcorr.lst))
    plot(sums.old, sums.new, log = "xy")
    print(m)
    print(m.umap.old)
    print(m.umap.new)
  dev.off()
})



# 
# # inf.old <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA/stringent_filter/B6_H3K4me3_pcutoff_0.95_binfilt_cellfilt.2019-06-03.stringent_filter.RData"
# inf.old <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA/B6_H3K4me1_pcutoff_0.95_binfilt_cellfilt.2019-05-11.RData"
# load(inf.old, v=T)
# mat.old <- count.dat$counts
# rownames(mat.old) <- paste("chr", rownames(mat.old), sep = "")
# 
# # load new mat
# # mat.new <- readRDS("/Users/yeung/data/scchic/quality_control_PZ_Bl6-BM-Linneg-Stemcells_AllMerged/B6BM_Unenriched_H3K4me3_pcutoff_0.5_binfilt_cellfilt.2019-11-27.rds")
# mat.new <- readRDS("/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-All_allmarks/PZ-Bl6-BM-All_Unenriched.H3K4me1.2019-11-23.rds")  # old antibody
# rows.com <- intersect(rownames(mat.new), rownames(mat.old))
# cols.com <- intersect(colnames(mat.new), colnames(mat.old))
# 
# print(nnzero(mat.new) / length(mat.new))
# print(nnzero(mat.old) / length(mat.old))
# 
# # compare sample by sampe??
# # jcell <- "B6-13W1-BM-H3K4me3-3_175"
# jcell <- cols.com[[1]]
# jcell <- cols.com[[10]]
# jdat <- as.data.frame(cbind(old = mat.old[rows.com, jcell], new = mat.new[rows.com, jcell]))
# ggplot(jdat, aes(x = old, y = new)) + geom_jitter() + ggtitle(jcell)
# 
# # hierarchical clusterinig???
# 
# # check H3K4me1
# 
# # try with common rows?
# lsi.out.old <- RunLSI(as.matrix(mat.old[rows.com, ]))
# lsi.out.new <- RunLSI(as.matrix(mat.new[rows.com, ]))
# 
# jsettings <- umap.defaults
# jsettings$n_neighbors <- 50
# jsettings$min_dist <- 0.1
# jsettings$random_state <- 123
# 
# umap.out.lsi.old <- umap(lsi.out.old$u, config = jsettings)
# dat.umap.long.lsi.old <- data.frame(cell = rownames(umap.out.lsi.old$layout), umap1 = umap.out.lsi.old$layout[, 1], umap2 = umap.out.lsi.old$layout[, 2], stringsAsFactors = FALSE)
# ggplot(dat.umap.long.lsi.old, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   ggtitle("Old count matrix")
# 
# umap.out.lsi.new <- umap(lsi.out.new$u, config = jsettings)
# dat.umap.long.lsi.new <- data.frame(cell = rownames(umap.out.lsi.new$layout), umap1 = umap.out.lsi.new$layout[, 1], umap2 = umap.out.lsi.new$layout[, 2], stringsAsFactors = FALSE)
# ggplot(dat.umap.long.lsi.new, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   ggtitle("New count matrix")
# 
# # merge the two together
# 
# # colnames(mat.new) <- paste(colnames(mat.new), "new", sep = ".")
# # colnames(mat.old) <- paste(colnames(mat.old), "old", sep = ".")
# 
# plot(density(log2(unlist(as.matrix(mat.old[rows.com, ])) + 1)))
# plot(density(log2(unlist(as.matrix(mat.new[rows.com, ])) + 1)))
# 
# # do global comparison??
# x <- unlist(mat.old[rows.com, cols.com])
# y <- unlist(mat.new[rows.com, cols.com])
     