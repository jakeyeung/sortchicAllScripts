# Jake Yeung
# Date of Creation: 2019-04-25
# File: ~/projects/scchic/scripts/scripts_analysis/sorted_bulk_analysis/simulate_cell_mixings.R
# Simulate
# reef: 
# https://rdrr.io/cran/flexmix/man/FLXmclust.html 
# https://cran.r-project.org/web/packages/flexmix/vignettes/flexmix-intro.pdf 
# # https://cran.r-project.org/web/packages/flexmix/flexmix.pdf 


rm(list=ls())


# Simulate? ---------------------------------------------------------------


library(dplyr)
library(ggplot2)
library(flexmix)
library(reshape2)
library(tidyr)
library(rstan)
library(brms)
# simulate data: 10 genes?
# 3 cell types

set.seed(0)
jsig <- 0.9
ngenes <- 10
nctypes <- 3
nsamps <- 1500
# ctypes.frac1 <- runif(n = nsamps, min = 0, max = 1)
# ctypes.frac2 <- runif(n = nsamps, min = 0, max = sapply(ctypes.frac1, function(x) 1 - x))
# a third is frac1, frac2, frac3

mhigh <- 10
mlow <- 1

ctypes.frac1.x <- rnorm(n = nsamps / 3, mean = mhigh, sd = 0.1)
ctypes.frac1.y <- rnorm(n = nsamps - length(ctypes.frac1.x), mean = mlow, sd = 0.1)
ctypes.frac1 <- c(ctypes.frac1.x, ctypes.frac1.y)
# ctypes.frac1 <- sample(x = c(0, 0.5, 1))
# ctypes.frac2 <- runif(n = nsamps, min = 0, max = 0.1)
ctypes.frac2.x <- rnorm(n = nsamps / 3, mean = mlow, sd = 0.1)
ctypes.frac2.y <- rnorm(n = nsamps / 3, mean = mhigh, sd = 0.1)
ctypes.frac2 <- c(ctypes.frac2.x, ctypes.frac2.y, ctypes.frac2.x)

ctypes.frac3.x <- rnorm(n = nsamps * 4 / 5, mean = mlow, sd = 0.1)
ctypes.frac3.y <- rnorm(n = nsamps - length(ctypes.frac3.x), mean = mhigh, sd = 0.1)
ctypes.frac3 <- c(ctypes.frac3.x, ctypes.frac3.y)

ctypes.frac.mat <- cbind(ctypes.frac1, ctypes.frac2, ctypes.frac3)
# normalize
ctypes.frac.mat <- sweep(exp(ctypes.frac.mat), MARGIN = 1, STATS = rowSums(exp(ctypes.frac.mat)), FUN = "/")

plot(ctypes.frac.mat[, 1], type = "l")
lines(ctypes.frac.mat[, 2], type = "l", col = 'blue')
lines(ctypes.frac.mat[, 3], type = "l", col = 'red')

clstrs.real <- apply(ctypes.frac.mat, MARGIN = 1, which.max)
names(clstrs.real) <- paste("samp", seq(nsamps), sep = "")

#check
assertthat::assert_that(all(rowSums(ctypes.frac.mat)) == 1)

markers <- matrix(FALSE, nrow = nctypes, ncol = ngenes)
markers[1, 1:4] <- TRUE
markers[2, 4:7] <- TRUE
markers[3, 7:10] <- TRUE

gene.betas <- rnorm(n = 10, mean = 5, sd = 3)  # celltype specific???

# set up archetypes
exprs.mat <- sweep(markers, MARGIN = 2, STATS = gene.betas, FUN = "*")

# assign clusters
clstrs <- sample(x = seq(3), size = nsamps, replace = TRUE)

# simulate exprs
exprs.sim <- matrix(data = NA, nrow = nsamps, ncol = ngenes)
for (i in seq(nsamps)){
  # mixings
  gvec <- colSums(sweep(exprs.mat, MARGIN = 1, STATS = as.vector(as.numeric(ctypes.frac.mat[i, ])), FUN = "*"))
  exprs.sim[i, ] <- gvec
}

# add noise
exprs.sim.withnoise <- exprs.sim + rnorm(length(exprs.sim), mean = 0, sd = jsig)
sizefacs <- rep(1.0, nrow(exprs.sim.withnoise))

# try to recover this
dat.pca <- prcomp(t(exprs.sim.withnoise), center = TRUE, scale. = TRUE)

plot(dat.pca$rotation[, 1], dat.pca$rotation[, 2])  # continuum


rownames(exprs.sim.withnoise) <- paste("samp", seq(nrow(exprs.sim.withnoise)), sep = "")
colnames(exprs.sim.withnoise) <- paste("gene", seq(ncol(exprs.sim.withnoise)), sep = "")

# add celltype marker information
colnames(markers) <- paste("gene", seq(ncol(markers)), sep = "")


exprs.long <- melt(exprs.sim.withnoise, varnames = c("samp", "gene"), value.name = "exprs") %>%
  group_by(samp) %>%
  mutate(exprs.cent = scale(exprs, center = TRUE, scale = FALSE))

for (i in seq(nctypes)){
  exprs.long[[paste0("exprs", i)]] <- exprs.long$exprs.cent * markers[i, ][exprs.long$gene]
  exprs.long[[paste0("des", i)]] <- as.numeric(markers[i, ][exprs.long$gene])
  # center it
  # exprs.long[[paste0("exprs", i)]] <- scale(exprs.long[[paste0("exprs", i)]], center = TRUE, scale = FALSE)
}

# frmla.lst <- lapply(seq(nctypes), function(i) return(as.formula(paste0(". ~ 1"))))
# mdl.lst <- lapply(frmla.lst, function(frmla) FLXMRglmfix(formula = frmla, family = "gaussian", fixed = ~ 1 + gene))
# mixout.ind <- flexmix(exprs.cent ~ 0 + gene | samp, exprs.long, k = 3, model = mdl.lst)
# plot(posterior(mixout.ind)[, 1])


# add marker gene
frmla.lst <- lapply(seq(nctypes), function(i) return(as.formula(paste0(". ~ 1 + gene:des", i))))
mdl.lst <- lapply(frmla.lst, function(frmla) FLXMRglmfix(formula = frmla, family = "gaussian", fixed = ~ 1 + gene))
mixout.ind <- flexmix(exprs.cent ~ 1 | samp, data = exprs.long, k = 3, model = mdl.lst)
(posterior(mixout.ind)[, 1])

# 
# # do some custom designs?
# # try mixtools 
# frmla.lst <- lapply(seq(nctypes), function(i) return(as.formula(paste0("exprs.cent ~ 0 + gene + gene:des", i))))
# mdl.mat.lst <- lapply(frmla.lst, function(x){
#   desmat <- model.matrix(x, data = exprs.long)
#   # remove empty columns
#   desmat <- desmat[, which(colSums(desmat) != 0)]
# })


# try regmixEM
  

library(mixtools)
k <- 3
Ylst <- lapply(split(exprs.long, exprs.long$samp), function(x) x$exprs.cent)
Xlst <- lapply(split(exprs.long, exprs.long$samp), function(x){
  model.matrix(exprs.cent ~ 1, x)
})
jlambda <- rep(1, k); jlambda <- jlambda / sum(jlambda)

sigma.hyp <<- 2  # mixtools 

system.time(
  em.out <- regmixEM.mixed(Ylst, Xlst, sigma = 2, arb.sigma = FALSE, arb.R = FALSE,
                           lambda = jlambda, k = 3, 
                           addintercept.random = FALSE,
                           epsilon = 1e-02, verb = TRUE)
)

em.out.bckup <- em.out
# outputs?
print(em.out)

plot(em.out)

# 
# 
# # Try brms ----------------------------------------------------------------
# 
# library(brms)
# 
# # fit1 <- brm(count ~ zAge + zBase * Trt + (1|patient), 
# #             data = epilepsy, family = poisson())
# 
# mix <- mixture(gaussian, gaussian, gaussian)
# 
# prior <- c(
#   prior(normal(0, 7), Intercept, dpar = mu1),
#   prior(normal(5, 7), Intercept, dpar = mu2),
#   prior(normal(10, 7), Intercept, dpar = mu3)
# )
# 
# fit1 <- brm(bf(exprs.cent ~ (1 | samp) + gene), exprs.long, family = mix,
#             prior = prior, chains = 2) 
# summary(fit1)
# pp_check(fit1)
# 
# 
# set.seed(1234)
# dat <- data.frame(
#   y = c(rnorm(200), rnorm(100, 6)), 
#   x = rnorm(300),
#   z = sample(0:1, 300, TRUE)
# )
# 
# ## fit a simple normal mixture model
# mix <- mixture(gaussian, gaussian, gaussian)
# jprior <- c(
#   brms::prior(normal(0, 7), Intercept, dpar = mu1),
#   brms::prior(normal(5, 7), Intercept, dpar = mu2),
#   brms::prior(normal(15, 7), Intercept, dpar = mu3)
# )
# fit1 <- brm(bf(exprs.cent ~ (1 | samp) + gene), exprs.long, family = ,
#             prior = jprior, chains = 2)
# 
# 
# fit1 <- brm(bf(y ~ x + z), dat, family = mix,
#             prior = prior, chains = 2) 
# summary(fit1)
# pp_check(fit1)
# 
# ## fit a simple normal mixture model
# mix <- mixture(gaussian, gaussian)
# jprior <- c(
#   brms::prior(normal(0, 7), Intercept, dpar = mu1),
#   brms::prior(normal(5, 7), Intercept, dpar = mu2)
# )
# fit1 <- brm(bf(y ~ x + z), dat, family = mix,
#             prior = jprior, chains = 2) 
# summary(fit1)
# pp_check(fit1)

