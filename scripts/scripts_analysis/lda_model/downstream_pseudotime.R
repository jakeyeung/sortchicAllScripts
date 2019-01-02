# Jake Yeung
# Date of Creation: 2019-01-01
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/downstream_pseudotime.R
# Downstream of poisson regression

library(dplyr)
library(ggplot2)

# load("~/projects/scChiC/outputs_R/fit_output/glm_fits_100kb_windows.Robj", v=T)  # glm.fits
load("~/projects/scChiC/outputs_R/fit_output/glm_fits_100kb_windows_withOffset.Robj", v=T)  # glm.fits

load("/private/tmp/lda_output/BM-H3K27me3.AvO_filt.Robj", v=T)
load("/private/tmp/lda_output/lda_out.meanfilt.K-12.Robj", v=T)

count.mat <- count.dat$counts

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


(hoxc.peaks <- grep("chr15:1029", rownames(count.mat), value = TRUE))
(hoxc.peaks.i <- grep("chr15:1029", rownames(count.mat), value = FALSE))

jpeak <- "chr2:74660000-74760000"
fits.again <- FitGlmRow(count.mat[jpeak, ], jpseudo, jsize, returnobj = TRUE)

pred.x <- runif(100, -5, 2)
pred.y <- exp(predict(fits.again, newdata = data.frame(pseudo = pred.x, size = mean(jsize))))

pred.x <- jpseudo
pred.y <- fits.again$fitted.values

jdat.real <- data.frame(x = jpseudo, y = count.mat[jpeak, ], type = "real")
# jdat.pred <- data.frame(x = pred.x, y = pred.y, type = "pred")
jdat.pred <- data.frame(x = pred.x, y = pred.y, type = "pred")

ggplot() + geom_point(data = jdat.real, aes(x = x, y = y), alpha=0.15) + 
  geom_line(data = jdat.pred, aes(x = x, y = y))


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
