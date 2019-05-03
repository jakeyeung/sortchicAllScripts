# Jake Yeung
# Date of Creation: 2019-04-30
# File: ~/projects/scchic/scripts/scripts_analysis/pseudobulk_analysis/pseudobulk_vs_sorted_again.R
# Look at H3K27me3 more closely

library(data.table)
library(ggplot2)

inf1 <- "/Users/yeung/data/scchic/from_cluster/pseudobulk_comparisons/merged_softlinks_textfile_log1p/H3K27me3_MatBcell_rep1_comparison.txt"
inf2 <- "/Users/yeung/data/scchic/from_cluster/pseudobulk_comparisons/merged_softlinks_textfile_log1p/H3K27me3_MatBcell_rep2_comparison.txt"

dat1 <- fread(inf1)
dat2 <- fread(inf2)

# merge the two replicates together???

bcell.vec <- apply(cbind(dat1$H3K27me3_MatBcell_rep1_reltoinput.bw, dat2$H3K27me3_MatBcell_rep2_reltoinput.bw), 1, mean)

# compare with clusters 2 and 5

compare.dat <- data.frame(bcell = bcell.vec, clstr2 = dat1$H3K27me3_cluster_2.log1p.bw, clstr5 = dat1$H3K27me3_cluster_5.log1p.bw)

cor(compare.dat$bcell, compare.dat$clstr2, use = "complete", method = "spearman")
cor(compare.dat$bcell, compare.dat$clstr5, use = "complete", method = "spearman")

jsub <- subset(compare.dat, clstr2 < 0.25 & clstr5 < 0.25 & clstr2 > 0.025 & clstr5 > 0.025)

cor(jsub$bcell, jsub$clstr2, use = "complete", method = "pearson")
cor(jsub$bcell, jsub$clstr5, use = "complete", method = "pearson")

plot(jsub$bcell, jsub$clstr2, pch = ".")
plot(jsub$bcell, jsub$clstr5, pch = ".")

# try to winsorize??
