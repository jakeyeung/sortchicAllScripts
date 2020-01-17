# Jake Yeung
# Date of Creation: 2020-01-14
# File: ~/projects/scchic/scripts/rstudioserver_analysis/debug/fit_GLM_test.R
# Try fit GLM

jstart <- Sys.Date()

library(scchicFuncs)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(glmpca)
library(topicmodels)

# Load and fit ------------------------------------------------------------

outdir <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs"
assertthat::assert_that(dir.exists(outdir))


inf <- "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_50000_50000.NoSliding/lda_outputs.count_mat.H3K4me1.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me1.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj"
load(inf, v=T)

tm.result <- posterior(out.lda)
dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

experis <- as.factor(sapply(colnames(count.mat), function(x) ClipLast(x, jsep = "_")))
experis.int <- as.integer(experis)

print(unique(experis))

Y <- as.matrix(count.mat)

terms.all <- colnames(tm.result$terms)
terms.keep <- apply(tm.result$terms, 1, function(x) terms.all[order(x, decreasing = TRUE)][1:1000])
terms.keep <- unique(as.vector(terms.keep))

Y <- Y[terms.keep, ]
# which(rowSums(Y) == 0)
# print(range(rowSums(Y)))
# print(range(colSums(Y)))

# clean up zeros
# which(colSums(Y) == 0)

X.dat <- data.frame(cell = colnames(Y))
X.dat <- left_join(X.dat, dat.var)
X <- matrix(X.dat$cell.var.within.sum.norm, ncol = 1, byrow = TRUE)  # intrachromo var

assertthat::assert_that(all(rowSums(Y) > 0))

print("Running GLM....")
system.time(
  glmout.regressvar <- glmpca(Y = Y, L = 30, fam = "mult", X = X)
)
outf.regressvar <- file.path(outdir, "fit_GLM_test_regressvar_glmout.RData")
print("Done running GLM...")
saveRDS(object = glmout.regressvar, file = outf.regressvar)

print("Running GLM.... no regress var")
system.time(
  glmout <- glmpca(Y = Y, L = 30, fam = "mult")
)
outf <- file.path(outdir, "fit_GLM_test_noregressvar_glmout.RData")
print("Done running GLM...")
saveRDS(object = glmout, file = outf)


print(Sys.Date() - jstart)