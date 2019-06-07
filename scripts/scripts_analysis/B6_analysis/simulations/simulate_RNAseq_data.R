# Jake Yeung
# Date of Creation: 2019-06-01
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/simulations/simulate_RNAseq_data.R
# Simulate RNA-seq data and do mean-variance relationships 

library(splatter)
simsim <- lun2Simulate()
dat <- sim@assays$data$counts

# do standard

jmean <- rowMeans(dat)
jvar <- apply(dat, 1, var)

cv2 <- jvar / jmean ^ 2

plot(log(jmean), log(cv2), main = "Summarize per gene across cells")
abline(a = 0, b = -1)
plot(log(jmean), log(jvar), main = "Summarize per gene across cells")


# now do the reverse?

jmean.cell <- colMeans(dat)
jvar.cell <- apply(dat, 2, var)
cv2.cell <- jvar.cell / jmean.cell ^ 2


plot(log(jmean.cell), log(jvar.cell), main = "Summarize per cell across genes")
plot(log(jmean.cell), log(cv2.cell), main = "Summarize per cell across genes")
plot(jmean.cell, sqrt(jvar.cell), main = "Summarize per cell across genes")
