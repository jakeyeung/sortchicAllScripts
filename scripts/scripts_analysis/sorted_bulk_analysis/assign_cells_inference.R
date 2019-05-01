# Jake Yeung
# Date of Creation: 2019-04-25
# File: ~/projects/scchic/scripts/scripts_analysis/sorted_bulk_analysis/assign_cells_inference.R
# Assign cells probabilisticaly kA


rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(umap)
library(tidytext)

# library(devtools)
# dev_mode(T)
# install_local("/Users/yeung/projects/cellassign", force=TRUE, upgrade = "never")
library(cellassign)


library(nnet)
library(msgl)
library(doParallel); library(foreach)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/TrajFunctions.R")


jmark <- "H3K4me1"

jfac <- 10^6
jpseudo <- 0

# Get trajs mixed ---------------------------------------------------------

trajs.mixed.out <- GetTrajMixed()
trajs.mixed <- trajs.mixed.out$trajs.mixed
dat.umap.mixed <- trajs.mixed.out$dat.umap.mixed


# Load bulk ---------------------------------------------------------------

inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"
dat <- fread(inf.bulkdat, sep = "\t")

colnames(dat) <- gsub(" ", "_", colnames(dat))

dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))

# Load UMAP ---------------------------------------------------------------

inf.umap <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
assertthat::assert_that(file.exists(inf.umap))
load(inf.umap, v=T)

# Load exprs  -------------------------------------------------------------

tssdist <- 50000
jdate <- "2019-04-22"
Kvec <- "50"
inf.lda <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_peaks_exprs_merge_clusters_explore/lda_out_meanfilt.PZ-BM-H3K4me1.CountThres0.K-", Kvec, "_GeneTSS.Dedup.", jdate, ".", tssdist, ".Robj")
# inf.lda <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_peaks_exprs_merge_clusters_explore/lda_out_meanfilt.PZ-BM-H3K4me1.CountThres0.K-50_GeneTSS.Dedup.2019-04-21.20000.Robj"

assertthat::assert_that(file.exists(inf.lda))
load(inf.lda, v=T)

out.lda <- out.lda[[length(out.lda)]]
out.lda@documents <- SwitchColnames(unname(out.lda@documents), jsplit = "-")

tm.result <- posterior(out.lda)

# get UMAP coords
umap.out <- umap(tm.result$topics)

umap.long <- data.frame(cell = unname(rownames(umap.out$layout)), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)

top.cells <- tidytext::tidy(out.lda, matrix = "gamma") %>%
  group_by(topic) %>%
  arrange(desc(gamma)) %>%
  mutate(rnk = seq(length(gamma))) %>%
  mutate(gamma.zscore = scale(gamma, center = TRUE, scale = TRUE)) %>%
  dplyr::rename(cell = document)

top.cells.sum <- top.cells %>%
  group_by(topic) %>% # do entropy on 1 to 99% of cells
  filter(gamma.zscore < quantile(gamma.zscore, 0.95)) %>%
  mutate(zscore.prob = exp(gamma.zscore) / sum(exp(gamma.zscore))) %>%
  summarise(entropy = -sum(zscore.prob * log(zscore.prob))) %>%
  arrange(entropy)

top.peaks <- tidytext::tidy(out.lda, matrix = "beta", log = FALSE) %>%
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  mutate(rnk = seq(length(beta))) %>%
  mutate(beta.zscore = scale(beta, center = TRUE, scale = TRUE)) %>%
  rowwise() %>%
  mutate(gene = strsplit(term, ";")[[1]][[2]])

mat.impute <- t(tm.result$topics %*% tm.result$terms)

# filter top terms
m.celldens <- ggplot(top.cells, aes(x = gamma.zscore)) + geom_density() + facet_wrap(~topic) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m.celldens)

m.umap <- ggplot(umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.umap)

m.umap.orig <- ggplot(dat.umap.long.trajs[[jmark]], aes(x = umap1, y = umap2)) + geom_point()

print(m.umap.orig) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# get gene list
rnames <- rownames(mat.impute)
rnames.keep <- grepl(";", rnames)

mat.impute.sub <- mat.impute[rnames.keep, ]

genes <- sapply(rownames(mat.impute.sub), function(x){
  g <- tryCatch({
    return(strsplit(x, ";")[[1]][[2]])
  }, error = function(e){
    return("Peak")
  })
}, USE.NAMES = FALSE)

exprs.long <- data.frame(peak = rownames(mat.impute.sub), gene = genes, as.data.frame(mat.impute.sub)) %>%
  tidyr::gather(key = "cell", value = "exprs", c(-peak, -gene))

exprs.long <- left_join(exprs.long, dat.umap.long.trajs[[jmark]] %>% dplyr::select(cell, umap1, umap2, louvain))


# Get gene expression vector for each cell? -------------------------------

# can we correlate cells?
# library(multinom)
# jfit <- readRDS("/Users/yeung/data/scchic/robjs/multinom_fit.0.999.rds")

exprs.long <- exprs.long %>%
  mutate(exprs.log = log2(exprs * jfac + jpseudo)) %>%
  group_by(cell) %>%
  mutate(zscore = scale(exprs.log, center = TRUE, scale = TRUE))

# add trajectory info

m1 <- PlotXYWithColor(exprs.long %>% filter(gene == "S100a8"), xvar = "umap1", yvar = "umap2", cname = "zscore")
print(m1)


# Select marker genes -----------------------------------------------------
ntops <- 8
topics.keep <- top.cells.sum$topic[1:ntops]
keepn <- 50
genes.keep <- subset(top.peaks, rnk < keepn & topic %in% topics.keep)$gene
genes.top1 <- subset(top.peaks, rnk == 1 & topic %in% topics.keep)$gene

# take only intersection of genes
genes.keep.intersect <- intersect(genes.keep, dat.long$Gene_Name)

# # plot topics on umap??
# for (g in genes.top1){
#   print(g)
#   m1 <- PlotXYWithColor(exprs.long %>% filter(gene == g), xvar = "umap1", yvar = "umap2", cname = "zscore", jtitle = g)
#   print(m1)
# }


# Create design matrix ----------------------------------------------------

ctypes.all <- as.character(unique(dat.long$CellType))
ctypes <- c("Kit_and_Sca1-positive_hematopoietic_stem_cell", "T_cell", "granulocyte", "nucleate_erythrocyte", "natural_killer_cell", "lymphocyte_of_B_lineage", "fetal_liver_hematopoietic_progenitor_cell", "megakaryocyte")

# get marker genes for each cell type
dat.sub <- subset(dat.long, CellType %in% ctypes & Gene_Name %in% genes.keep.intersect) %>%
  group_by(Gene_Name) %>%
  mutate(is.marker = ifelse(zscore > 1.5, TRUE, FALSE))


# Use cell assign ---------------------------------------------------------


# try on real data???
dat <- subset(exprs.long, gene %in% genes.keep.intersect)

# cell by gene
dat.mat <- tidyr::spread(dat %>% dplyr::select(gene, cell, zscore) %>% mutate(gene = paste("gene", gene, sep = "_")), key = gene, value = zscore) %>%
  as.data.frame()
rownames(dat.mat) <- dat.mat$cell; dat.mat$cell <- NULL

marker.genes <- tidyr::spread(dat.sub %>% dplyr::select(Gene_Name, CellType, is.marker) %>% mutate(gene = paste("gene", Gene_Name, sep = "_")), key = CellType, value = is.marker) %>%
  ungroup() %>%
  dplyr::select(-Gene_Name) %>%
  as.data.frame()
rownames(marker.genes) <- marker.genes$gene; marker.genes$gene <- NULL
assertthat::assert_that(all(rownames(marker.genes) == colnames(dat.mat)))


# make into matrices
dat.mat.mat <- as.matrix(dat.mat)
marker.genes.mat <- as.matrix(marker.genes)
# matrix(as.integer(marker.genes), nrow = nrow(marker.genes), ncol = ncol(marker.genes))

sizefacs <- rep(1.0, nrow(dat.mat))

system.time(
  cas <- cellassign(exprs_obj = dat.mat.mat,
                    marker_gene_info = marker.genes.mat,
                    s = sizefacs)  
)

saveRDS(cas, file = paste0("~/data/scchic/robjs/cellassign_output.ntop.", ntops, ".keepn.", keepn, ".rds"))
dev_mode(F)


exprs.umap <- subset(exprs.long, gene == "Tal1")

cas.long <- data.frame(cell = rownames(dat.mat.mat), as.data.frame(cas$mle_params$gamma), stringsAsFactors = FALSE)
cas.long <- left_join(cas.long, exprs.umap %>% dplyr::select(cell, umap1, umap2))

pdf(paste0("~/data/scchic/pdfs/assign_cells_output.ntops.", ntops, ".keepn.", keepn,  ".pdf"), useDingbats = FALSE)
PlotXYWithColor(cas.long, xvar = "umap1", yvar = "umap2", cname = "nucleate_erythrocyte")
PlotXYWithColor(cas.long, xvar = "umap1", yvar = "umap2", cname = "T_cell")
PlotXYWithColor(cas.long, xvar = "umap1", yvar = "umap2", cname = "natural_killer_cell")
PlotXYWithColor(cas.long, xvar = "umap1", yvar = "umap2", cname = "lymphocyte_of_B_lineage")
PlotXYWithColor(cas.long, xvar = "umap1", yvar = "umap2", cname = "granulocyte")
PlotXYWithColor(cas.long, xvar = "umap1", yvar = "umap2", cname = "fetal_liver_hematopoietic_progenitor_cell")
PlotXYWithColor(cas.long, xvar = "umap1", yvar = "umap2", cname = "Kit_and_Sca1.positive_hematopoietic_stem_cell")
PlotXYWithColor(cas.long, xvar = "umap1", yvar = "umap2", cname = "megakaryocyte")
dev.off()


# Try unsupervised way  ---------------------------------------------------

library(flexmix)

# row center
dat.mat.rowcent <- sweep(dat.mat.mat, MARGIN = 1, STATS = rowMeans(dat.mat.mat), FUN = "-")
mixout <- flexmix(dat.mat.rowcent ~ 1, k = 8, model = FLXMCmvnorm(diagonal = TRUE), control = list(classify = "weighted"))
mixweights <- posterior(mixout)

# plot the outputs
mixout.long <- data.frame(cell = rownames(dat.mat.rowcent), as.data.frame(mixweights), stringsAsFactors = FALSE) %>%
  left_join(., exprs.umap %>% dplyr::select(cell, umap1, umap2))

pdf("~/data/scchic/pdfs/flexmix_mvnorm.pdf", useDingbats = FALSE)
cnames <- grep("^V[0-9]", colnames(mixout.long), value = TRUE)
for (cname in cnames){
  print(cname)
  print(PlotXYWithColor(mixout.long, xvar = "umap1", yvar = "umap2", cname = cname))
}
dev.off()

# can we try desmat?
dat.mat.long <- tidyr::gather(as.data.frame(dat.mat.rowcent) %>% mutate(cell = rownames(dat.mat.rowcent)), key = gene, value = exprs, -cell)

system.time(
  mixout2 <- flexmix(formula = exprs ~ 1 + gene | cell, data = dat.mat.long, k = 7, 
                     cluster = apply(posterior(mixout), 1, which.max), 
                     model = FLXMRglm(family = "gaussian"))
)

# check my design matrix is sound
dat.mat.long.markers <- left_join(dat.mat.long, marker.genes %>% mutate(gene = rownames(marker.genes))) %>%
  mutate(cell = as.factor(cell),
         gene = as.factor(gene))

# add marker genes
des.mat <- model.matrix(exprs ~ 1 + gene : granulocyte, dat.mat.long.markers)
des.mat <- model.matrix(exprs ~ 1 + I(gene * granulocyte), dat.mat.long.markers)






# Simulate? ---------------------------------------------------------------



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


# maybe use the fixed and varying effects?
rownames(exprs.sim.withnoise.rowcent) <- paste("samp", seq(nrow(exprs.sim.withnoise.rowcent)), sep = "")
colnames(exprs.sim.withnoise.rowcent) <- paste("gene", seq(ncol(exprs.sim.withnoise.rowcent)), sep = "")

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

frmla.lst <- lapply(seq(nctypes), function(i) return(as.formula(paste0(". ~ 1"))))
mdl.lst <- lapply(frmla.lst, function(frmla) FLXMRglmfix(formula = frmla, family = "gaussian", fixed = ~ 1 + gene))
mixout.ind <- flexmix(exprs.cent ~ 0 + gene | samp, exprs.long, k = 3, model = mdl.lst)
plot(posterior(mixout.ind)[, 1])

# add marker gene
frmla.lst <- lapply(seq(nctypes), function(i) return(as.formula(paste0(". ~ 1"))))
mdl.lst <- lapply(frmla.lst, function(frmla) FLXMRglmfix(formula = frmla, family = "gaussian", fixed = ~ 1 + gene:des1))
mixout.ind <- flexmix(exprs.cent ~ 1 | samp, exprs.long, k = 3, model = mdl.lst)
plot(posterior(mixout.ind)[, 1])


# do some custom designs?
# try mixtools 
frmla.lst <- lapply(seq(nctypes), function(i) return(as.formula(paste0("exprs.cent ~ 0 + gene + gene:des", i))))
mdl.mat.lst <- lapply(frmla.lst, function(x){
  desmat <- model.matrix(x, data = exprs.long)
  # remove empty columns
  desmat <- desmat[, which(colSums(desmat) != 0)]
})



# create 

MakeMatRemoveEmptyCols <- function(dat, frmla){
  desmat <- model.matrix(frmla, data = dat)
  return(desmat[, which(colSums(desmat) != 0)])
}
k <- 3
Y <- lapply(split(exprs.long, exprs.long$samp), function(x) return(x$exprs.cent))
# X <- lapply(split(exprs.long, exprs.long$samp), function(x) MakeMatRemoveEmptyCols(x, as.formula("exprs.cent ~ 0 + gene + gene:des1 + gene:des2 + gene:des3")))
X <- lapply(split(exprs.long, exprs.long$samp), function(x) model.matrix(exprs.cent ~ 1, x))
lambda.start <- rep(1, k); lambda.start <- lambda.start / sum(lambda.start)
R.set <- lapply(seq(Y), function(i) diag(1, k))
# mixem.out <- regmixEM.mixed(Y, X, lambda = lambda.start, R = R.set, addintercept.fixed = TRUE)
mu.mat <- matrix(runif(n = 30), nrow = nrow(X[[1]]), ncol = k)
# mixem.out <- regmixEM.mixed(Y[1:2], X[1:2], lambda = lambda.start, addintercept.fixed = FALSE, k=k, addintercept.random = TRUE, verb = 1, mu = mu.mat, sigma = 2, arb.sigma = FALSE, arb.R = FALSE)
sigma.hyp <<- 2
mixem.out <- regmixEM.mixed(Y, X, lambda = lambda.start, addintercept.fixed = TRUE, k=k, arb.sigma = FALSE, arb.R = FALSE)

library(mixtools)  
data(RanEffdata)
set.seed(100)
x <- lapply(1:length(RanEffdata), function(i) 
+             matrix(RanEffdata[[i]][, 2:3], ncol = 2))
x <- x[1:20]
y <- lapply(1:length(RanEffdata), function(i) 
+             matrix(RanEffdata[[i]][, 1], ncol = 1))
y <- y[1:20]
lambda <- c(0.45, 0.55)
mu <- matrix(c(0, 4, 100, 12), 2, 2)
sigma <- 2
R <- list(diag(1, 2), diag(1, 2))

# each list if a sample
em.out <- regmixEM.mixed(y, x, sigma = sigma, arb.sigma = FALSE,
                                                     lambda = lambda, mu = mu, R = R, k = 2,
                                                     addintercept.random = FALSE,
                                                     epsilon = 1e-02, verb = TRUE)

# # plot the hits
# head(cas$mle_params$gamma)
# 
# cas.long <- 

# 
# # 
# library(flexmix)
# 
# # fit two classes 
# 
# data("NPreg", package = "flexmix")
# 
# Model_n <- FLXMRglm(yn ~ . + I(x^2))
# Model_p <- FLXMRglm(yn ~ .)
# m1 <- flexmix(. ~ x, data = NPreg, k = 2, model = list(Model_n, Model_p), control = list(verbose = 10))
# 
# # try on real data???
# dat <- subset(exprs.long, gene %in% genes.keep.intersect)
# 
# # cell by gene
# dat.mat <- tidyr::spread(dat %>% dplyr::select(gene, cell, zscore) %>% mutate(gene = paste("gene", gene, sep = "_")), key = gene, value = zscore)
# rownames(dat.mat) <- dat.mat$cell; dat.mat$cell <- NULL
# 
# marker.genes <- tidyr::spread(dat.sub %>% dplyr::select(Gene_Name, CellType, is.marker) %>% mutate(gene = paste("gene", Gene_Name, sep = "_")), key = CellType, value = is.marker) %>%
#   ungroup() %>%
#   dplyr::select(-Gene_Name) %>%
#   as.data.frame()
# rownames(marker.genes) <- marker.genes$gene; marker.genes$gene <- NULL
# 
# 
# 
# assertthat::assert_that(all(rownames(marker.genes) == colnames(dat.mat)))
# 
# # make into matrices
# dat.mat.mat <- as.matrix(dat.mat)
# marker.genes.mat <- matrix(as.integer(marker.genes), nrow = nrow(marker.genes), ncol = ncol(marker.genes))
# 
# sizefacs <- rep(1.0, nrow(dat.mat))
# library(devtools)
# dev_mode(T)
# install_local("/Users/yeung/projects/cellassign", force=TRUE)
# library(cellassign)
# cas <- cellassign(exprs_obj = dat.mat.mat,
#                   marker_gene_info = marker.genes.mat,
#                   s = sizefacs)
# dev_mode(F)
# 
# 
# 
# # Y <- dat$zscore
# # # do celltype-specific design matrix
# 
# # simulate data
# 
# 
# 
# # simulate
# # 
# # library(brms)
# # 
# # set.seed(1234)
# # dat <- data.frame(
# #   y = c(rnorm(200), rnorm(100, 6)), 
# #   x = rnorm(300),
# #   z = sample(0:1, 300, TRUE)
# # )
# # 
# # ## fit a simple normal mixture model
# # mix <- mixture(gaussian, gaussian)
# # prior <- c(
# #   prior(normal(0, 7), Intercept, dpar = mu1),
# #   prior(normal(5, 7), Intercept, dpar = mu2)
# # )
# # fit1 <- brm(bf(y ~ x + z), dat, family = mix,
# #             prior = prior, chains = 2) 
# # summary(fit1)
# # pp_check(fit1)
# # 
# # ## use different predictors for the components
# # fit2 <- brm(bf(y ~ 1, mu1 ~ x, mu2 ~ z), dat, family = mix,
# #             prior = prior, chains = 2) 
# # summary(fit2)
# 
# # 
# # 
# # 
# # # check example
# # library(mixtools)
# # data("CO2data")
# # attach(CO2data)
# # CO2reg <- regmixEM(CO2, GNP, lambda = c(1, 3) / 4, beta = matrix(c(8, -1, 1, 1), 2, 2), sigma = c(2, 1), addintercept = TRUE)
# 
