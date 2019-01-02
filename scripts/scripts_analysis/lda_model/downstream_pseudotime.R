# Jake Yeung
# Date of Creation: 2019-01-01
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/downstream_pseudotime.R
# Downstream of poisson regression

library(dplyr)
library(ggplot2)


# Functions ---------------------------------------------------------------

source("scripts/Rfunctions/FitFunctions.R")


# Load --------------------------------------------------------------------



# load("~/projects/scChiC/outputs_R/fit_output/glm_fits_100kb_windows.Robj", v=T)  # glm.fits
# load("~/projects/scChiC/outputs_R/fit_output/glm_fits_100kb_windows_withOffset.Robj", v=T)  # glm.fits
# load("/private/tmp/lda_output/BM-H3K27me3.AvO_filt.Robj", v=T)
# load("/private/tmp/lda_output/lda_out.meanfilt.K-12.Robj", v=T)


load("~/projects/scChiC/outputs_R/fit_output/glm_fits_peak_calling_windows_withOffset.Robj", v=T)  # glm.fits
# load("/private/tmp/lda_output/lda_outputs.meanfilt_1.merge_1000_NoXYM.cellmin_1000.cellmax_50000/lda_out.meanfilt.K-12.Robj", v=T)
load("/private/tmp/lda_output/lda_outputs.meanfilt_1.merge_1000.cellmin_1000.cellmax_50000/lda_out.meanfilt.K-12.Robj", v=T)
load("/private/tmp/lda_output/PZ-BM-H3K27me3.merged.NoCountThres.Robj", v=T)


count.mat <- count.dat$counts

# needed sometimes
rownames(out.lda@gamma) <- out.lda@documents
cells.keep <- out.lda@documents
peaks.keep <- out.lda@terms
cells.keep.i <- which(colnames(count.mat) %in% cells.keep)
peaks.keep.i <- which(rownames(count.mat) %in% peaks.keep)
count.mat <- count.mat[peaks.keep.i, cells.keep.i]

# Summarize it ------------------------------------------------------------

fits.sum <- lapply(glm.fits, function(x) return(data.frame(int = x$int,
                                                           pseudo = x$pseudo,
                                                           pval = x$pval)))
jnames <- names(fits.sum)
fits.sum <- bind_rows(fits.sum)
fits.sum$region <- jnames
fits.sum <- fits.sum %>%
  arrange(pval)


# Plot outputs ------------------------------------------------------------

jregions <- grep("chr15:1029", fits.sum$region, value = TRUE)

fits.sum %>% 
  filter(region %in% jregions)

print(head(fits.sum %>% arrange(pseudo)))

(hoxc.peaks <- grep("chr15:1029", rownames(count.mat), value = TRUE))
(hoxc.peaks.i <- grep("chr15:1029", rownames(count.mat), value = FALSE))

jpeak <- "chr2:74660000-74760000"

jpeak <- "chrY:90811097-90813994"
jpeak <- "chr5:28459482-28460573"

jpeak <- "chr16:34825504-34827428"

jpeak <- "chrY:90811097-90813994"

fits.again <- FitGlmRow(count.mat[jpeak, ], jpseudo, jsize, returnobj = TRUE)

# pred.x <- runif(100, min(jpseudo), max(jpseudo))
# pred.y <- exp(predict(fits.again, newdata = data.frame(pseudo = pred.x, size = mean(jsize))))

pred.x <- jpseudo
pred.y <- fits.again$fitted.values

jdat.real <- data.frame(x = jpseudo, y = count.mat[jpeak, ], type = "real")
# jdat.pred <- data.frame(x = pred.x, y = pred.y, type = "pred")
jdat.pred <- data.frame(x = pred.x, y = pred.y, type = "pred")

ggplot() + geom_point(data = jdat.real, aes(x = x, y = y), alpha=0.15) + 
  geom_line(data = jdat.pred, aes(x = x, y = y)) + theme_bw() + ggtitle(jpeak)

# dat$phat <- m1.pois$fitted.values
# ggplot(dat) + geom_point(aes(x = pseudo, y = counts)) + geom_line(aes(x = pseudo, y = phat))

pred.x <- runif(100, 0, 0.25)
jpeak <- "chr12:114660000-114760000"
jpeak <- "chr15:102900000-103000000"
jpeak <- "chr2:74660000-74760000"
jpeak <- "chr11:96420000-96520000"
jpeak <- "chr1:39800000-39900000"
jpeak <- "chr3:57540000-57640000"
jpeak <- "chr20:170700000-170800000"
jsub <- subset(fits.sum, region == jpeak)
pred.y <- exp(jsub$int + jsub$pseudo* pred.x)

jdat.real <- data.frame(x = jpseudo, y = count.mat[jpeak, ], type = "real")
jdat.pred <- data.frame(x = pred.x, y = pred.y, type = "pred")

ggplot() + geom_point(data = jdat.real, aes(x = x, y = y), alpha=0.15) + geom_line(data = jdat.pred, aes(x = x, y = y))
