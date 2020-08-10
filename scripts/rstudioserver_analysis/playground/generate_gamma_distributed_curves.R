# Jake Yeung
# Date of Creation: 2020-06-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/playground/generate_gamma_distributed_curves.R
# Gamma

rm(list=ls())

jshape <- 1.5
jscale <- 5
N <- 1000

x <- rgamma(N, shape = jshape, scale = jscale)

plot(density(x))

rpois(n = N, lambda = jscale)

x <- rgeom(N, 0.3)
y <- rgeom(N, 0.7)

plot(hist(x))
plot(hist(y))
plot(hist(c(x, y)))


jmean <- mean(x)
jvar <- var(x)
jskew <- moments::skewness(x)
jkurt <- moments::kurtosis(x)

jmean / jvar ^ 0.375
jmean / jskew ^ 0.2
jmean / jkurt ^ 0.125


