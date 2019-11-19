# Jake Yeung
# Date of Creation: 2019-06-04
# File: ~/projects/scchic/scripts/scripts_analysis/simulations/multinomial_likelihood.R
# Multinomial likelihood space

# Generate data -----------------------------------------------------------

ngenes <- 10000
ncells <- 50
ncounts <- 30
ntranscripts <- 10^6

p1 <- exp(rnorm(n = ngenes, mean = 0, sd = 1))
p2 <- exp(rnorm(n = ngenes, mean = 0, sd = 1))
p3 <- mapply(function(x, y) sqrt(x * y), p1, p2)

#1 generate single-cell transcriptomes 

clstr1 <- rmultinom(ncells, ncounts, prob = p1)
clstr2 <- rmultinom(ncells, ncounts, prob = p2)
clstr3 <- rmultinom(ncells, ncounts, prob = p3)

dat.merged <- cbind(clstr1, clstr2, clstr3)

dat.pca <- prcomp(dat.merged)

colvec <- c(rep("blue", ncells), rep("red", ncells), rep("black", ncells))
plot(dat.pca$rotation[, 1], dat.pca$rotation[, 2], col = colvec, pch = 20, xlab = "PC1 Loadings", ylab = "PC2 Loadings", cex = 5, cex.lab = 1.5, cex.main = 1.5, main = "PCA on Observed Space")

# Get likelihoods ---------------------------------------------------------

ll1 <- apply(dat.merged, 2, function(x) dmultinom(x = x, prob = p1, log = TRUE))
ll2 <- apply(dat.merged, 2, function(x) dmultinom(x = x, prob = p2, log = TRUE))

plot(ll1, ll2, col = colvec, pch = 20, xlab = "Loglikelihood 1", ylab = "Loglikelihood 2", cex = 5, cex.lab = 1.5, cex.main = 1.5, main = "Cells in Latent Space")

