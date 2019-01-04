# Jake Yeung
# run_cisTopic.R
# Run cisTopic 
# 2019-01-03

jstart <- Sys.time() 

# library(topicmodels)
library(dplyr)
library(ggplot2)
# library(ldatuning)
library(cisTopic)


print(paste("Work directory: ", getwd()))

source("scripts/Rfunctions/ParseStrings.R")
source("scripts/Rfunctions/Aux.R")

args <- commandArgs(trailingOnly=TRUE)

print("input args:")
print(args)

inpath <- args[[1]]
outdir <- args[[2]]
nclst <- StrToNumeric(args[[3]])
topic.vec <- as.numeric(StrToVector(args[[4]], delim = ","))
meanmax <- StrToNumeric(args[[5]])  # remove suspicious peaks 
cellmin <- StrToNumeric(args[[6]])  # remove cells with low counts
cellmax <- StrToNumeric(args[[7]])  # remove suspiciious cells
projname <- args[[8]]

if (is.na(nclst)){
  stop(paste("nclst must be numeric, found", nclst))
}
if (is.na(meanmax)){
  stop(paste("meanmax must be numeric, found", meanmax))
}
print(paste("Will iterate through", length(topic.vec), "Ks"))
print(topic.vec)

plotpath <- file.path(outdir, "plots.meanfilt.pdf")
outpath <- file.path(outdir, paste0("lda_out.meanfilt.K-", nclst, ".Robj"))
tunepath <- file.path(outdir, paste0("lda_tuning.meanfilt.K-", nclst, ".Robj"))


# constants
# countmax <- 100  # peaks with more than these counts are filtered out for suspicious 
# this works for small peaks but probably needs to be increased for larger peaks?
# meanmax <- 1  # peaks with more than these counts are filtered out for suspicious 
# cellmin <- 

# print args
print("Args:")
print(paste(inpath, outdir, nclst, topic.vec, cellmin, cellmax))

print(topic.vec)

# stop('debugging')

# Load counts -------------------------------------------------------------

load(inpath, v=T)

count.mat <- count.dat$counts

print(dim(count.mat))

# # Plot mean and variance --------------------------------------------------
# 
dat.meanvar <- data.frame(Sum = Matrix::rowSums(count.mat), 
                          Mean = Matrix::rowMeans(count.mat),
                          Var = apply(count.mat, 1, var),
                          peak = rownames(count.mat),
                          stringsAsFactors=FALSE)
dat.meanvar <- dat.meanvar %>%
  rowwise() %>%
  mutate(CV = sqrt(Var) / Mean,
         peaksize = GetPeakSize(peak))
# 
# p1 <- ggplot(dat.meanvar, aes(x = peaksize)) + geom_density() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# p2 <- ggplot(dat.meanvar, aes(x = log10(Mean), y = log10(CV), size = peaksize)) + geom_point(alpha = 0.25) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   geom_abline(slope = -0.5)
# 
# p3 <- ggplot(dat.meanvar, aes(x = peaksize, y = Sum)) + geom_point(alpha = 0.1) +
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_y_log10()
# 
# # save plots
# pdf(plotpath, useDingbats = FALSE)
#   print(p1)
#   print(p2)
#   print(p3)
# dev.off()

# suspicious peaks
bad.peaks <- dat.meanvar %>%
  filter(Mean > meanmax)
print(head(bad.peaks$peak))
print(paste("There are", nrow(bad.peaks), "peaks with more than", "counts. Removing them..."))

# filter them out before running count.mat
print("Dimensions before filtering peaks...")
print(dim(count.mat))

count.mat <- count.mat[which(!rownames(count.mat) %in% bad.peaks$peak), ]
print("Dimensions after filtering peaks...")
print(dim(count.mat))

# remove peaks in chrM
# M.peaks <- grep("chrM|chrY|chrX", dat.meanvar$peak, value=TRUE)
M.peaks <- grep("chrM", dat.meanvar$peak, value=TRUE)

count.mat <- count.mat[which(!rownames(count.mat) %in% M.peaks), ]
print("Dimensions after filtering peaks X, Y, M chromos")
print(dim(count.mat))


# Remove cells with zero entries
print("Dimensions before filtering cells")
print(dim(count.mat))
count.mat <- count.mat[, which(Matrix::colSums(count.mat) > cellmin)]
count.mat <- count.mat[, which(Matrix::colSums(count.mat) < cellmax)]
print("Dimensions after filtering cells...")
print(dim(count.mat))

# Run LDA on count matrix -------------------------------------------------

# nclst <- 10
print("Running cisTopic")
cisTopicObject <- createcisTopicObject(count.mat, project.name=projname)
cisTopicObject <- runModels(cisTopicObject, topic=topic.vec, seed=0, nCores=length(topic.vec), burnin = 120, iterations = 175, addModels=FALSE)
save(cisTopicObject, file = file.path(outdir, paste0(projname, "_cisTopicOutput.Robj")))

print(Sys.time() - jstart)
