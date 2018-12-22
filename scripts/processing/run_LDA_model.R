# Jake Yeung
# run_LDA_model.R
# Run LDA model 
# 2018-12-19

jstart <- Sys.time() 

library(topicmodels)
library(dplyr)
library(ggplot2)
library(ldatuning)


print(paste("Work directory: ", getwd()))

source("scripts/Rfunctions/ParseStrings.R")

args <- commandArgs(trailingOnly=TRUE)

print("input args:")
print(args)

inpath <- args[[1]]
outdir <- args[[2]]
nclst <- StrToNumeric(args[[3]])
topic.vec <- StrToVector(args[[4]], delim = ",")
tunemodels <- StrToBool(args[[5]])
meanmax <- StrToNumeric(args[[6]])  # remove suspicious peaks 
cellmin <- StrToNumeric(args[[7]])  # remove cells with low counts
cellmax <- StrToNumeric(args[[8]])  # remove suspiciious cells 

if (is.na(nclst)){
  stop(paste("nclst must be numeric, found", nclst.str))
}
if (is.na(meanmax)){
  stop(paste("meanmax must be numeric, found", meanmax.str))
}
if (is.na(tunemodels)){
  stop(paste("tunemodels must be TRUE or FALSE, found", tunemodels.str))
}
print(paste("Will iterate through", length(topic.vec), "Ks"))
print(topic.vec)

plotpath <- file.path(outdir, "plots.meanfilt.pdf")
outpath <- file.path(outdir, paste0("lda_out.meanfilt.K-", nclst, ".Robj"))
tunepath <- file.path(outdir, paste0("lda_tuning.meanfilt.K-", nclst.str, ".Robj"))

GetPeakSize <- function(coord){
  # chr1:3005258-3006803 -> 1545
  coord <- as.character(coord)
  jstart <- as.numeric(strsplit(strsplit(coord, ":")[[1]][[2]], "-")[[1]][[1]])
  jend <- as.numeric(strsplit(strsplit(coord, ":")[[1]][[2]], "-")[[1]][[2]])
  return(jend - jstart)
}

# constants
# countmax <- 100  # peaks with more than these counts are filtered out for suspicious 
# this works for small peaks but probably needs to be increased for larger peaks?
# meanmax <- 1  # peaks with more than these counts are filtered out for suspicious 
cellmin <- 

# Load counts -------------------------------------------------------------

load(inpath, v=T)

count.mat <- count.dat$counts

# Plot mean and variance --------------------------------------------------

dat.meanvar <- data.frame(Sum = Matrix::rowSums(count.mat), 
                          Mean = Matrix::rowMeans(count.mat),
                          Var = apply(count.mat, 1, var),
                          peak = rownames(count.mat),
                          stringsAsFactors=FALSE)
dat.meanvar <- dat.meanvar %>%
  rowwise() %>%
  mutate(CV = sqrt(Var) / Mean,
         peaksize = GetPeakSize(peak))

p1 <- ggplot(dat.meanvar, aes(x = peaksize)) + geom_density() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


p2 <- ggplot(dat.meanvar, aes(x = log10(Mean), y = log10(CV), size = peaksize)) + geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(slope = -0.5)

p3 <- ggplot(dat.meanvar, aes(x = peaksize, y = Sum)) + geom_point(alpha = 0.1) +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_log10()

# save plots
pdf(plotpath, useDingbats = FALSE)
  print(p1)
  print(p2)
  print(p3)
dev.off()

# suspicious peaks
bad.peaks <- dat.meanvar %>%
  filter(Mean > meanmax)
print(head(bad.peaks$peak))
print(paste("There are", nrow(bad.peaks), "peaks with more than", "counts. Removing them..."))

# filter them out before running count.mat
print("Dimensions before filtering...")
print(dim(count.mat))

count.mat <- count.mat[which(!rownames(count.mat) %in% bad.peaks$peak), ]
print("Dimensions after filtering...")
print(dim(count.mat))

# Remove cells with zero entries
count.mat <- count.mat[, which(Matrix::colSums(count.mat) > 0)]

# Run LDA on count matrix -------------------------------------------------

# nclst <- 10
print("Running LDA")
out.lda <- LDA(x = t(as.matrix(count.mat)), k = nclst, method = "Gibbs", control=list(seed=0))

# save output
print("Saving LDA")
save(out.lda, file = outpath)
print("Time elapsed after LDA")
print(Sys.time() - jstart)

# tune LDA 

if (tunemodels){
    # topic.vec <- c(4, 9, 11, 14, 16, 18)
    print("Running tuning")
    optimal.topics <- FindTopicsNumber(t(as.matrix(count.mat)), topics=topic.vec, mc.cores = length(topic.vec), method="Gibbs", metrics=c("Arun2010", "CaoJuan2009", "Griffiths2004", "Deveaud2014"), control = list(seed=0))
    FindTopicsNumber_plot(optimal.topics)

    save(optimal.topics, file = tunepath)
    print("Saving optimal topics")
}

print("Time elapsed after tuning")
print(Sys.time() - jstart)

