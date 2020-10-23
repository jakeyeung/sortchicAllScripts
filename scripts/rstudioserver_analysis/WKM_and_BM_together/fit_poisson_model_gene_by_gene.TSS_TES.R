# Jake Yeung
# Date of Creation: 2020-06-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/fit_poisson_model_gene_by_gene.R
# 


rm(list=ls())

jstart <- Sys.time()

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(hash)
library(JFuncs)
library(scchicFuncs)


FitGlmRowClusters <- function(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, jbin = NULL, returnobj=FALSE){
  # use Offset by size of library
  # https://stats.stackexchange.com/questions/66791/where-does-the-offset-go-in-poisson-negative-binomial-regression
  # fit GLM for a row of a sparse matrix, should save some space?
  
  # pvalue by deviance goodness of fit: https://thestatsgeek.com/2014/04/26/deviance-goodness-of-fit-test-for-poisson-regression/
  # offset is in log because the model says the log counts is equal to RHS
  
  if (!is.null(nrow(jrow))){
    # probably a matrix of many rows, sum them up
    print(paste("Merging", nrow(jrow), "rows"))
    row <- Matrix::colSums(jrow)
  }
  dat <- data.frame(cell = cnames, ncuts = jrow, stringsAsFactors = FALSE) %>%
    left_join(., dat.annots.filt.mark, by = "cell") %>%
    left_join(., ncuts.cells.mark, by = "cell")
  
  # m1.pois <- glm(ncuts ~ 1 + Cluster + offset(ncuts.total), data = dat, family = "poisson")
  m1.pois <- glm(ncuts ~ 1 + Cluster + offset(log(ncuts.total)), data = dat, family = "poisson")
  mnull.pois <- glm(ncuts ~ 1 + offset(log(ncuts.total)), data = dat, family = "poisson")
  
  if (!returnobj){
    jsum <- anova(mnull.pois, m1.pois)
    pval <- pchisq(jsum$Deviance[[2]], df = jsum$Df[[2]], lower.tail = FALSE)
    out.dat <- data.frame(pval = pval, 
                          dev.diff = jsum$Deviance[[2]],
                          df.diff = jsum$Df[[2]],
                          t(as.data.frame(coefficients(m1.pois))), 
                          stringsAsFactors = FALSE)
    if (!is.null(jbin)){
      out.dat$bin <- jbin
      rownames(out.dat) <- jbin
    }
    return(out.dat)
  } else {
    return(list(fit.full = m1.pois, fit.null = mnull.pois, dat.input = dat))
  }
}



fewer.k27me3 <- TRUE

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression"
dir.create(indir, showWarnings = TRUE)
jprefix <- file.path(indir, paste0("integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_", fewer.k27me3, ".forPoissonRegression.CountR1only.", "2020-06-05"))

outdir <- indir
outfits <- file.path(outdir, paste0("fit_poisson_model_on_TSS.RData"))

infrdata <- paste0(jprefix, ".smaller.RData")

assertthat::assert_that(file.exists(infrdata))

load(infrdata, v=T)

dat.annots.filt.forfit <- lapply(dat.annots.filt, function(jdat){
  jdat <- subset(jdat, select = c(cell, cluster.new)) %>%
    mutate(Cluster = ifelse(cluster.new == "HSPCs", "aHSPCs", cluster.new))  # set HSPC as intercept
  return(jdat)
})

# system.time(
#   jlong.cellfilt.counts <- lapply(jmarks, function(jmark){
#     jmat <- tss.mats.filt.fromref.cellfilt[[jmark]]
#     # set up the data frame for melting
#     jdat <- data.frame(bin = rownames(jmat), as.data.frame(as.matrix(jmat)), stringsAsFactors = FALSE)
#     colnames(jdat) <- c("bin", colnames(jmat))
#     # annotate bins
#     jdat$mark <- jmark
#     jdat$gene <- sapply(as.character(jdat$bin), function(b) strsplit(b, ";")[[1]][[2]])
#     jdat$ens <- sapply(jdat$gene, function(g) AssignHash(g, g2e, null.fill = g))
#     jlong <- reshape2::melt(jdat, id.vars = c("bin", "mark", "gene", "ens"), variable.name = "cell", value.name = "ncuts")
#     jlong <- left_join(jlong, ncuts.cells[[jmark]])  # add total counts
#     jlong <- left_join(jlong, subset(dat.annots.filt[[jmark]], select = c(cell, cluster.new)))
#     return(jlong)
#   })
# )

# lapply(jlong.cellfilt.counts, dim)


# 
# # Run poisson model  ------------------------------------------------------
# 
# jmark <- "H3K27me3"
# # jmark <- "H3K4me3"
# jgene <- "Hlf"
# jgene <- "Pax5"
# jgene <- "Sox6"
# jgene <- "S100a7a"
# # jbin.i <- sample(rownames(tss.mats.filt.fromref.cellfilt[[jmark]]), size = 1)
# 
# jbin <- rownames(tss.mats.filt.fromref.cellfilt[[jmark]])[grepl(jgene, rownames(tss.mats.filt.fromref.cellfilt[[jmark]]))]
# 
# jrow <- tss.mats.filt.fromref.cellfilt[[jmark]][jbin, ]
# cnames <- colnames(tss.mats.filt.fromref.cellfilt[[jmark]])
# dat.annots.filt.mark <- subset(dat.annots.filt[[jmark]], select = c(cell, cluster.new)) %>%
#   mutate(Cluster = ifelse(cluster.new == "HSPCs", "aHSPCs", cluster.new))  # set HSPC as intercept
# ncuts.cells.mark <- ncuts.cells[[jmark]]
# 
# jfit <- FitGlmRowClusters(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, returnobj = TRUE)
# jout <- FitGlmRowClusters(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, jbin = jbin, returnobj = FALSE)
# 
# print(jfit$fit.full)
# 
# jsum <- anova(jfit$fit.null, jfit$fit.full)
# print(jsum)
# 
# pchisq(jfit$fit.full$deviance, jfit$fit.full$df.residual, lower.tail = TRUE)
# pval <- pchisq(jsum$Deviance[[2]], df = jsum$Df[[2]], lower.tail = FALSE)
# 
# params.dat <- 
# 
# out.dat <- data.frame(pval = pval, 
#                       dev.diff = jsum$Deviance[[2]],
#                       df.diff = jsum$Df[[2]],
#                       t(as.data.frame(coefficients(jfit$fit.full))), 
#                       stringsAsFactors = FALSE)
# rownames(out.dat) <- jbin
# 
# 
# 
# # Plot fits ---------------------------------------------------------------
# 
# pred.cell <- jfit$dat.input$cell
# pred.y <- jfit$fit.full$fitted.values
# 
# pred.y.linear <- predict(jfit$fit.full, newdata = jfit$dat.input, type = "response")
# pred.y.log <- predict(jfit$fit.full, newdata = jfit$dat.input)
# 
# pred.dat <- data.frame(ypred = pred.y.log, jfit$dat.input, stringsAsFactors = FALSE)
# 
# ggplot(pred.dat, aes(x = Cluster, y = ypred, color = log(ncuts.total))) + 
#   geom_jitter(width = 0.2, height = 0)  + 
#   scale_color_viridis_c() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   ggtitle(paste(jgene, jbin, jmark)) + 
#   ylab("logLambda")


# Fit genome-wide ---------------------------------------------------------

# jmark <- "H3K4me3"
# jmat.mark <- tss.mats.filt.fromref.cellfilt[[jmark]]
# dat.annots.filt.mark <- dat.annots.filt.forfit[[jmark]]
# ncuts.cells.mark <- ncuts.cells[[jmark]]
# cnames <- colnames(jmat.mark)

print("Fitting... ")

system.time(
  jfits.lst.bymark <- lapply(jmarks, function(jmark){
    print(jmark)
    jmat.mark <- tss.mats.filt.fromref.cellfilt[[jmark]]
    dat.annots.filt.mark <- dat.annots.filt.forfit[[jmark]]
    ncuts.cells.mark <- ncuts.cells[[jmark]]
    cnames <- colnames(jmat.mark)
    jfits.lst <- apply(jmat.mark, MARGIN = 1, FUN = function(jrow){
      jout <- FitGlmRowClusters(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, jbin = NULL, returnobj = FALSE)
      return(jout)
    })
    return(jfits.lst)
  })
)

save(tss.mats.filt.fromref.cellfilt, dat.annots.filt.forfit, ncuts.cells, jfits.lst.bymark, file = outfits)

print(Sys.time() - jstart)



