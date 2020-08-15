# Jake Yeung
# Date of Creation: 2020-08-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/dinuc_analysis/init_dinuc_analysis.R
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Granu -------------------------------------------------------------------

inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/dinuc_freq/H3K4me3-BM_AllMerged.Granulocytes.sorted.cleaned.txt"
dat <- fread(inf, sep = "\t", header = FALSE)

# jlength <- nchar(dat$V2[[1]])
jlength <- 200
dinuc.mat <- matrix(data = 0, nrow = nrow(dat), ncol = jlength)

# pos strands only 
pos.i <- sapply(dat$V1, function(x) endsWith(x, suffix = "+"))

dat.sub <- dat[pos.i, ]

system.time(
  dinuc.mat <- sapply(dat.sub$V2, function(seq){
    svec <- strsplit(seq, split = "")[[1]]
    countvec <- sapply(svec, function(s){
      if (s %in% c("A", "T")){
        return(1)
      } else {
        return(0)
      }
    }, USE.NAMES = FALSE)
    return(countvec[1:jlength])
  }, USE.NAMES = FALSE)
)

dinuc.mat2 <- apply(dinuc.mat, 2, FUN = function(jcol){
  rollapply(jcol, 2, FUN = prod)
})

dinuc.sum <- apply(dinuc.mat2, MARGIN = 2, mean)

plot(dinuc.sum[1:100])
