# Jake Yeung
# Date of Creation: 2018-12-22
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/downstream_LDA_merged_peaks.R
# Analyze LDA output on merged peaks from 25kb distance
# If this works it serves as input to metacell so we can compare

library(dplyr)
library(ggplot2)

source("scripts/Rfunctions/Aux.R")

# Download data from server -----------------------------------------------

outdir <- "/tmp/lda_outputs.meanfilt.merge_25000"
if (!dir.exists(outdir)){
  # passwords not handled automatically, so this fails at the moment
  jcmd <- paste0("scp -r t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_outputs.meanfilt.merge_25000 ", outdir)
  system(jcmd)
}

# download count matrix
count.inf <- file.path(outdir, "PZ-BM-H3K4me1.merged.NoCountThres.merge_25000.Robj")
if (!file.exists(count.inf)){
  print("Downloading count matrix")
  jcmd2 <- paste0("scp -r t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats.merge_25000/PZ-BM-H3K4me1.merged.NoCountThres.merge_25000.Robj ", count.inf)
  system(jcmd2)
}

# Load stuff --------------------------------------------------------------

load(count.inf, v=T)
count.mat <- count.dat$counts

# Plot output -------------------------------------------------------------

dat.meanvar <- data.frame(Sum = Matrix::rowSums(count.mat), 
                          Mean = Matrix::rowMeans(count.mat),
                          Var = apply(count.mat, 1, var),
                          peak = rownames(count.mat), 
                          stringsAsFactors = FALSE)
dat.meanvar <- dat.meanvar %>%
  rowwise() %>%
  mutate(CV = sqrt(Var) / Mean,
         peaksize = GetPeakSize(peak)) 

# maybe a few suspicious peaks?
ggplot(dat.meanvar, aes(x = peaksize, y = Sum)) + geom_point(alpha = 0.1) +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_log10()

ggplot(dat.meanvar, aes(x = peaksize)) + geom_density() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meanvar, aes(x = log2(Mean), y = log2(CV), size = peaksize)) + geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(slope = -0.5) + ggtitle("Includes some suspicious peaks")

# plot cell sizes 
plot(density(colSums(count.mat)))
plot(hist(colSums(count.mat), breaks = 75))

# Just go for it LDA ------------------------------------------------------

# clean up bad cells
cellmin.thres <- 2000
cellmax.thres <- 27500
smallcells.i <- Matrix::colSums(count.mat) < cellmin.thres
bigcells.i <- Matrix::colSums(count.mat) > cellmax.thres

print(paste("N smallcells: ", length(which(smallcells.i))))
print(paste("N bigcells: ", length(which(bigcells.i))))

plot(hist(colSums(count.mat), breaks = 75), col = 'lightblue')
abline(v = c(cellmin.thres, cellmax.thres))

# remove bad cells
print(dim(count.mat))
count.mat <- count.mat[, !as.logical(smallcells.i + bigcells.i)]
print(dim(count.mat))

# remove bad regions
nclst <- 15
out.lda <- LDA(x = t(as.matrix(count.mat)), k = nclst, method = "Gibbs")
save(out.lda, count.mat, file="outputs_R/lda_output/lda_outputs.meanfilt.merge_25000.Robj")
