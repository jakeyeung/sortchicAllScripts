# Jake Yeung
# Date of Creation: 2019-04-25
# File: ~/projects/scchic/scripts/scripts_analysis/sorted_bulk_analysis/simulate_cell_mixings.R
# Simulate
# reef: 
# https://rdrr.io/cran/flexmix/man/FLXmclust.html 
# https://cran.r-project.org/web/packages/flexmix/vignettes/flexmix-intro.pdf 
# # https://cran.r-project.org/web/packages/flexmix/flexmix.pdf 


rm(list=ls())


library(dplyr)
library(ggplot2)
library(flexmix)

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
plot(dat.pca$rotation[, 2], dat.pca$rotation[, 3])  # continuum
# plot(dat.pca$x[, 1], dat.pca$x[, 2])

# Infer parameters --------------------------------------------------------
# 
# ex2 <- flexmix(yn ~ x, data = NPreg, k = 2,
#                model = list(FLXMRglm(yn ~ . + I(x^2)),
#                             FLXMRglm(yp ~ ., family = "poisson")))
# plot(ex2)
# 
# ex2.test <- flexmix(yn ~ 1, data = NPreg, k = 2,
#                model = list(FLXMRglm(yn ~ x + I(x^2)),
#                             FLXMRglm(yp ~ x, family = "poisson")))
# plot(ex2.test)
# 
# # observations are grouped
# colnames(exprs.sim.withnoise) <- paste("y", seq(ncol(exprs.sim.withnoise)), sep = "")
# 
# # fit 10 genes???
# ex2 <- flexmix(yn ~ x, data = NPreg, k = 2,
#                model = list(FLXMRglm(yn ~ . + I(x^2)),
#                             FLXMRglm(yp ~ ., family = "poisson")))
# 

# markers
library(devtools)
dev_mode(T)
# detach(name = "package:cellassign", unload = TRUE)
# remove.packages("cellassign")
# install_local("/Users/yeung/projects/cellassign", force=TRUE, upgrade = "never")
# install_github("Irrationone/cellassign")
# detach(name = "package:cellassign", unload = TRUE)
library(cellassign)
# cas <- cellassign(exprs_obj = scale(t(scale(t(exprs.sim.withnoise), center = TRUE, scale = TRUE)), center = TRUE, scale = TRUE), 
cas <- cellassign(exprs_obj = exprs.sim.withnoise, 
                  marker_gene_info = t(markers),
                  s = sizefacs)
dev_mode(F)

print(cas$mle_params$phi)


# check errors 
rss <- rowSums((ctypes.frac.mat - cas$mle_params$gamma) ^ 2)

plot(density(rss))

ctypes.frac.mat[which(rss > 0.3), ] 
cas$mle_params$gamma[which(rss > 0.3), ] 


# plot gammas 
plot(cas$mle_params$gamma[, 1], type = "l")
lines(cas$mle_params$gamma[, 2], type = "l", col = 'blue')
lines(cas$mle_params$gamma[, 3], type = "l", col = 'red')


# plot original
plot(ctypes.frac.mat[, 1], type = "l")
lines(ctypes.frac.mat[, 2], type = "l", col = 'blue')
lines(ctypes.frac.mat[, 3], type = "l", col = 'red')

plot(density(cas$mle_params$gamma))



# lets try FlexMix --------------------------------------------------------

library(flexmix)

# 
# data("Nclus", package = "flexmix")
# 
# require("MASS")
# eqscplot(Nclus)
# ex1 <- flexmix(Nclus ~ 1, k = 4, model = FLXMCmvnorm())
# print(ex1)
# plotEll(ex1, Nclus)
# 
# ex2 <- flexmix(exprs.sim.withnoise ~ 1, k = 3, model = FLXMCmvnorm())
# 
# library(tidyr)
# # https://rsangole.netlify.com/post/finite-mixture-modeling-using-flexmix/ 
# wh_best <- ex2
# wh_best.param <- parameters(wh_best)
# # remove cov
# wh_best.param <- wh_best.param[grepl("center", rownames(wh_best.param)), ]
# wh_best.param <- data.frame(Brand=rownames(wh_best.param),
#                             wh_best.param,row.names = NULL) 
# wh_best.param <- wh_best.param %>% gather(Components,Value,Comp.1:Comp.3)
# # wh_best.param <- wh_best.param %>% left_join(y = whiskey_brands,by = 'Brand')
# ggplot(wh_best.param,aes(y=Value,x=Brand))+
#   geom_bar(stat='identity')+
#   coord_flip()+
#   facet_grid(.~Components)
# 

# do it on our data

exprs.sim.withnoise.rowcent <- sweep(exprs.sim.withnoise, MARGIN = 1, STATS = rowMeans(exprs.sim.withnoise), FUN = "-")
# mixout <- flexmix(exprs.sim.withnoise.rowcent ~ 1, k = 3, model = FLXMCmvnorm(), control = list(classify = "weighted"))


mixout <- flexmix(exprs.sim.withnoise.rowcent ~ 1, k = 3, model = FLXMCmvnorm(), control = list(classify = "weighted"))
mixweights <- posterior(mixout)

plot(mixweights)
# where are my component means?
# plot gammas 
plot(mixweights[, 1], type = "l")
lines(mixweights[, 2], type = "l", col = 'blue')
lines(mixweights[, 3], type = "l", col = 'red')


plot(ctypes.frac.mat[, 1], type = "l")
lines(ctypes.frac.mat[, 2], type = "l", col = 'blue')
lines(ctypes.frac.mat[, 3], type = "l", col = 'red')


# add a design matrix
# 3 models depending on cluster 


# mo1 <- FLXMRglm(formula = . ~ . + gene , family = "gaussian")
# mo2 <- FLXMRglm(formula = . ~ . + gene , family = "gaussian")
# flexfit <- flexmix(x ~ 1, data = data, k = 2, model = list(mo1, mo2))
# 
# mixout.ind <- flexmix(exprs.sim.withnoise.rowcent ~ 1, k = 3, model = FLXMCmvnorm(formula = ~ . ~ + "gene1" + "gene2"))
# # mixout.des <- flexmix(exprs.sim.withnoise.rowcent ~ 1, k = 3, model = FLXMCmvnorm(formula = ), control = list(classify = "weighted"))
# 

# maybe use the fixed and varying effects?
rownames(exprs.sim.withnoise.rowcent) <- paste("samp", seq(nrow(exprs.sim.withnoise.rowcent)), sep = "")
colnames(exprs.sim.withnoise.rowcent) <- paste("gene", seq(ncol(exprs.sim.withnoise.rowcent)), sep = "")

rownames(exprs.sim.withnoise) <- paste("samp", seq(nrow(exprs.sim.withnoise)), sep = "")
colnames(exprs.sim.withnoise) <- paste("gene", seq(ncol(exprs.sim.withnoise)), sep = "")

library(reshape2)
exprs.long <- melt(exprs.sim.withnoise, varnames = c("samp", "gene"), value.name = "exprs") %>%
  group_by(samp) %>%
  mutate(exprs.cent = scale(exprs, center = TRUE, scale = FALSE))

system.time(
  mixed.long <- flexmix(formula = exprs.cent ~ gene | samp, data = exprs.long, k = 3, model = FLXMRglmfix(fixed = ~ 1))
)

print(summary(mixed.long))
print(parameters(mixed.long))
mixed.long.post <- data.frame(exprs.long$samp, posterior(mixed.long))
mixed.long.post <- mixed.long.post[!duplicated(mixed.long.post), ]
plot(mixed.long.post$X1, type = "l", ylim = c(0, 1))
lines(mixed.long.post$X2, type = "l", col = "blue")
lines(mixed.long.post$X3, type = "l", col = "red")


# add celltype marker information
colnames(markers) <- paste("gene", seq(ncol(markers)), sep = "")

for (i in seq(nctypes)){
  exprs.long[[paste0("exprs", i)]] <- exprs.long$exprs.cent * markers[i, ][exprs.long$gene]
  exprs.long[[paste0("des", i)]] <- as.numeric(markers[i, ][exprs.long$gene])
  # center it
  # exprs.long[[paste0("exprs", i)]] <- scale(exprs.long[[paste0("exprs", i)]], center = TRUE, scale = FALSE)
}

# change y
frmla.lst <- lapply(seq(nctypes), function(i) return(as.formula(paste0("exprs", i, " ~ gene"))))
mdl.lst <- lapply(frmla.lst, function(frmla) FLXMRglmnet(formula = frmla, family = "gaussian"))

system.time(
  mixed.long.filt <- flexmix(formula = . ~ 1 | samp, data = exprs.long, k = 3, model = mdl.lst)
)
print(summary(mixed.long.filt))
print(parameters(mixed.long.filt))
mixed.long.filt.post <- data.frame(exprs.long$samp, posterior(mixed.long.filt))
mixed.long.filt.post <- mixed.long.filt.post[!duplicated(mixed.long.filt.post), ]

plot(mixed.long.filt.post$X1, type = 'l')
lines(mixed.long.filt.post$X2, type = 'l', col = 'blue')
lines(mixed.long.filt.post$X3, type = 'l', col = 'red')


# try another way: change x
frmla.lst.des <- lapply(seq(nctypes), function(i) return(as.formula(paste0(". ~ gene:des", i))))
mdl.lst.des <- lapply(frmla.lst.des, function(frmla) FLXMRglmfix(formula = frmla, family = "gaussian", fixed = ~ gene))

system.time(
  mixed.long.filt.des <- flexmix(formula = exprs.cent ~ 1 | samp, data = exprs.long, k = 3, model = mdl.lst.des, cluster = clstrs.real[as.character(exprs.long$samp)])
)
print(summary(mixed.long.filt.des))
print(parameters(mixed.long.filt.des))
mixed.long.filt.des.post <- data.frame(exprs.long$samp, posterior(mixed.long.filt.des))
mixed.long.filt.des.post <- mixed.long.filt.post[!duplicated(mixed.long.filt.des), ]




mixed.long.filt.des2 <- flexmix(formula = exprs.cent ~ 1 | samp, data = exprs.long, k = 3,
                                model = FLXMRglmnet(. ~ gene*des1, fixed = ~ gene),
                                control = list(verb = 10, iter = 5, nrep = 5, classify="hard"))
                                # cluster = clstrs.real[as.character(exprs.long$samp)])


fo <- sample(rep(seq(10), length = nrow(exprs.long)))


system.time(
mixed.long.filt.des2 <- flexmix(formula = exprs.cent ~ 1 + gene | samp, data = exprs.long, k = 3,
                                model = FLXMRglm(family = "gaussian"),
                                cluster = clstrs.real[as.character(exprs.long$samp)])
)

print(summary(mixed.long.filt.des2))
print(parameters(mixed.long.filt.des2))
mixed.long.filt.des2.post <- data.frame(exprs.long$samp, posterior(mixed.long.filt.des2))
mixed.long.filt.des2.post <- mixed.long.filt.des2.post[!duplicated(mixed.long.filt.des2.post), ]

plot(mixed.long.filt.des2.post$X1, type = "l")
lines(mixed.long.filt.des2.post$X2, type = "l", col = 'blue')
lines(mixed.long.filt.des2.post$X3, type = "l", col = 'red')


plot(ctypes.frac.mat[, 1], type = "l")
lines(ctypes.frac.mat[, 2], type = "l", col = 'blue')
lines(ctypes.frac.mat[, 3], type = "l", col = 'red')



# 
# # try keeping same y, adding different design
# 
# exprs.long$des1 <- 
# 
# 
# 
# m1 <- 0
# m2 <- 50
# sd1 <- sd2 <- 5
# N1 <- 100
# N2 <- 10
# 
# a <- rnorm(n=N1, mean=m1, sd=sd1)
# b <- rnorm(n=N2, mean=m2, sd=sd2)
# x <- c(a,b)
# class <- c(rep('a', N1), rep('b', N2))
# data <- data.frame(cbind(x=as.numeric(x), class=as.factor(class)))
# 
# library("ggplot2")
# p <- ggplot(data, aes(x = x)) + 
#   geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", fill = "white") +
#   geom_vline(xintercept = m1, col = "red", size = 2) + 
#   geom_vline(xintercept = m2, col = "blue", size = 2)
# p
# 
# 
# mo1 <- FLXMRglm(family = "gaussian")
# mo2 <- FLXMRglm(family = "gaussian")
# flexfit <- flexmix(x ~ 1, data = data, k = 2, model = list(mo1, mo2))
# 
# 
# print(table(clusters(flexfit), data$class))
# 
# c1 <- parameters(flexfit, component=1)[[1]]
# c2 <- parameters(flexfit, component=2)[[1]]
# plot_mix_comps <- function(x, mu, sigma, lam) {
#   lam * dnorm(x, mu, sigma)
# }
# 
# lam <- table(clusters(flexfit))
# 
# ggplot(data) +
#   geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", fill = "white") +
#   stat_function(geom = "line", fun = plot_mix_comps,
#                 args = list(c1[1], c1[2], lam[1]/sum(lam)),
#                 colour = "red", lwd = 1.5) +
#   stat_function(geom = "line", fun = plot_mix_comps,
#                 args = list(c2[1], c2[2], lam[2]/sum(lam)),
#                 colour = "blue", lwd = 1.5) +
#   ylab("Density")

