# Jake Yeung
# Date of Creation: 2020-07-18
# File: ~/projects/scchic/scripts/rstudioserver_analysis/cuts_analysis/K562_temperature_differences_on_cut_distances.R
# Explore K562 cut distances


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(zoo)
library(Matrix)


# Plot oscillations in different temperatures -----------------------------

xstart <- 100
jtemp <- "20C"
# jtemp <- "4C"
jmeth <- "sap"
# jmeth <- "EtOH"
# inf.counts <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/temperature_exp/raw_data/tagged_bams.scmo3/CutDistances/PZ-K562-", jmeth, "-H3K4me3-", jtemp, "-1.scmo3.again_contigfixed.tagged/counts.csv")
inf.counts <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed/CutDistances/PZ-K562-H3K4me3-spikein-75U.scmo3.again_contigfixed.tagged/counts.csv"
assertthat::assert_that(file.exists(inf.counts))

jcounts <- fread(inf.counts)
jcounts$V1 <- NULL
jcell <- which(colSums(jcounts) == max(colSums(jcounts)))
x <- jcounts[xstart:nrow(jcounts), ..jcell]

xsmooth <- zoo::rollapply(x, width = 50, FUN = mean, align = "left")
plot(xsmooth, type = "l", main = paste(jmeth, jtemp))




# points(x, pch = 20)

# plot(log(x + 1), pch=20)



maxcounts <- 3
cellmaxcounts <- 100
jcounts.filt <- as.matrix(jcounts)
# remove some cells
print(dim(jcounts.filt))
jcounts.filt <- jcounts.filt[, which(colSums(jcounts.filt) > cellmaxcounts)]
print(dim(jcounts.filt))




cmeans <- colMeans(jcounts.filt[20:nrow(jcounts.filt), ])
# cmeans <- colMeans(jcounts.filt[1:nrow(jcounts.filt), ])
jcounts.filt <- sweep(jcounts.filt, MARGIN = 2, STATS = cmeans, FUN = "/")

# plot(unlist(jcounts[20:nrow(jcounts), ]))

jcounts.filt[which(jcounts.filt > maxcounts, arr.ind = TRUE)] <- maxcounts

jcounts.filt <- jcounts.filt[, which(!is.nan(colSums(jcounts.filt)))]

heatmap3::heatmap3(t(jcounts.filt), Rowv = NA, Colv = NA, scale = "none", RowSideLabs = "TotalCuts", main = "test", margins = c(5, 20), col = viridis::viridis(50, end = 1))

plot(unlist(jcounts.filt[, 1]))
plot(unlist(jcounts.filt[, 2]))



