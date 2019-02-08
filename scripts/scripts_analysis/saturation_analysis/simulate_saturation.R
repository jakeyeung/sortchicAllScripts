# Jake Yeung
# Date of Creation: 2019-02-04
# File: ~/projects/scchic/scripts/scripts_analysis/saturation_analysis/simulate_saturation.R
# Saturdation by simulation

set.seed(0)

N <- 1000  # total number of molecules

jsd <- 1  # unevenness of distribution, set to 0 for exponential saturation

# distribution of counts (strength of magnet to your hand)
weight.vec <- rlnorm(N, 0, sd = jsd)

# weight.vec <- rep(1, N)  # equal weights recapitulates the exponential curve

# distribution of weights 
plot(density(weight.vec))

# Sample ------------------------------------------------------------------

nreads <- seq(1, 10*N, by = 1)
jsamps <- lapply(nreads, function(nread) sample(seq(N), size = nread, replace = TRUE, prob = weight.vec))
umis <- unlist(lapply(jsamps, function(x) length(unique(x))))

# plot number of counts per molecule for read depth of X*N
plot(hist(table(jsamps[[3*N]])), xlab = "counts")

# median(table(jsamps[[3*N]]))

# Plot sat curve ----------------------------------------------------------

# fit the exponential
jdat <- data.frame(x = nreads, y = umis)
jfit <- nls(formula = y ~ A * (1 - exp(-x / k)), 
            data = jdat, 
            start = list(A = N, k = 1000))
pred.fit <- predict(jfit)
# fit spline
spl <- smooth.spline(nreads, y=umis)
pred.spline <- predict(spl)

plot(nreads, umis, pch = 20, cex = 0.5, main = paste("Fit from exponential, sig=", jsd))
lines(pred.fit, col='red')
abline(h = N)

plot(nreads, umis, pch = 20, cex = 0.5, main = paste("Fit from spline, sig=", jsd))
lines(pred.spline, col='orange')
abline(h = N)

pred.fit.prime <- diff(pred.fit) / diff(jdat$x)
plot(jdat$x[2:nrow(jdat)], pred.fit.prime, type = "l", main = paste("dy/dx from exponential. sig=", jsd))

pred.spline.prime <- predict(spl, deriv = 1)
plot(pred.spline.prime$x, log(pred.spline.prime$y), type = "l", main = paste("dy/dx from spline. sig=", jsd))
