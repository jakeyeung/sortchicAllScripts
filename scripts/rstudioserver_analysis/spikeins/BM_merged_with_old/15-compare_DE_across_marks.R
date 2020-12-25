# Jake Yeung
# Date of Creation: 2020-12-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/15-compare_DE_across_marks.R
# Check DE across marks

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)



# Load DE outputs ---------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

de.out <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.de <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3/poisson_fit_TES.", jmark, ".2020-11-14.newannot2.RData")
  load(inf.de, v=T)
  # create fold changes lreatiev to HSPCs
  effects.lst <- lapply(jfits.lst, function(x){
    indx <- grepl(pattern = "Cluster", names(x))
    # intercept <- c("ClusterHSPCs" = 0)
    xvec <- data.frame("ClusterHSPCs" = 0, x[indx])
    return(xvec)
  }) 
  effects.mat <- as.data.frame(effects.lst %>% bind_rows())
  rownames(effects.mat) <- names(jfits.lst)
  
  colnames(effects.mat) <- paste(jmark, colnames(effects.mat), sep = "_")
  
  pval.lst <- lapply(jfits.lst, function(x){
    return(x$pval)
  }) %>%
    unlist()
  
  pval.dat <- data.frame(rname = names(pval.lst), pval = pval.lst, stringsAsFactors = FALSE)
  colnames(pval.dat) <- c("rname", paste0("pval_", jmark))
  return(list(effects.mat = effects.mat, pval.dat = pval.dat, jmat.mark = jmat.mark))
})

# combine mats
rownames.common <- Reduce(intersect, lapply(de.out, function(x) rownames(x$effects.mat)))

pvalthres <- 1e-10
# rows.keep <- lapply(de.out, function(x) subset(x$pval.dat, rname %in% rownames.common))
pval.mat.merge <- Reduce(left_join, lapply(de.out, function(x) subset(x$pval.dat, rname %in% rownames.common))) %>%
  as.data.frame()

rownames(pval.mat.merge) <- pval.mat.merge$rname

pval.mat.merge$rname <- NULL

rows.keep <- rownames(pval.mat.merge)[which(apply(pval.mat.merge, 1, function(jrow) min(jrow) < pvalthres))]

# effects.mat.merge <- cbind(de.out$H3K4me1$effects.mat[rows.keep, ], de.out$H3K4me3$effects.mat[rows.keep, ])
effects.mat.merge <- do.call(cbind, lapply(de.out, function(x) x$effects.mat[rows.keep, ]))

# remove rows with effect sizes beyond +/- 10
good.rows <- rownames(effects.mat.merge)[apply(effects.mat.merge, 1, function(jrow) max(abs(jrow)) < 5)]

plot(density(unlist(effects.mat.merge[good.rows, ])))

rows.keep2 <- intersect(good.rows, rows.keep)
pca.out <- prcomp(t(effects.mat.merge[rows.keep2, ]), center = TRUE, scale. = TRUE)
plot(pca.out$x[, 1], pca.out$x[, 2], pch = 20)
text(pca.out$x[, 1], pca.out$x[, 2], labels = rownames(pca.out$x))


# normalize mark by mar 
effects.mat.merge.renorm <- cbind(t(scale(t(de.out$H3K4me1$effects.mat[rows.keep2, ]), center = TRUE, scale = TRUE)), t(scale(t(de.out$H3K4me3$effects.mat[rows.keep2, ]), center = TRUE, scale = TRUE)))
# effects.mat.merge.renorm <- cbind(t(scale(t(de.out$H3K4me1$effects.mat[rows.keep2, ]), center = TRUE, scale = TRUE)), t(scale(t(de.out$H3K4me3$effects.mat[rows.keep2, ]), center = TRUE, scale = TRUE)), -t(scale(t(de.out$H3K27me3$effects.mat[rows.keep2, ]), center = TRUE, scale = TRUE)))

pca.out2 <- prcomp(t(effects.mat.merge.renorm), center = FALSE, scale. = FALSE)
plot(pca.out2$x[, 1], pca.out2$x[, 2], pch = 20)
text(pca.out2$x[, 1], pca.out2$x[, 2], labels = rownames(pca.out2$x))

# Check some well known genes ---------------------------------------------

# jgene.grep <- "Hbb-y"
jgene.grep <- "Ebf1"
jgene.grep <- "Hlf"
jgene.grep <- "S100a8"
jgene.grep <- "Dach1"
jgene.grep <- "Sox6"
jgene.grep <- "Tal1"
jgene.grep <- "Cln8"
jgene.grep <- "Ltf"
jgene.grep <- "Cdh20"
jgene.grep <- "Klhl14"
jgene.grep <- "Pld4"
jgene.grep <- "Inpp4a"
jgene.grep <- "Galnt2"
jgene.grep <- "Siglech"
jgene.grep <- "Fyn"
rname <- grep(jgene.grep, rownames(pval.mat.merge), value = TRUE)[[1]]
print(pval.mat.merge[rname, ])

plot(unlist(effects.mat.merge[rname, ]), main = rname)
text(unlist(effects.mat.merge[rname, ]), labels = colnames(effects.mat.merge))

plot(unlist(effects.mat.merge.renorm[rname, ]), main = rname)
text(unlist(effects.mat.merge.renorm[rname, ]), labels = colnames(effects.mat.merge.renorm))


# Get celltype specific gene set for K4me1 and check in K27me3  -----------


# get granu-specific

plot(density(-log10(pval.mat.merge$pval_H3K4me1)))
plot(density(-log10(pval.mat.merge$pval_H3K4me3)))
plot(density(-log10(pval.mat.merge$pval_H3K27me3)))

# jmat.ref <- t(scale(t(de.out$H3K4me3$effects.mat[rows.keep, ])))

# mark.ref <- "H3K4me1"
mark.ref <- "H3K27me3"
jmat.ref <- de.out[[mark.ref]]$effects.mat[rows.keep, ]
# handle outliers
jmat.ref[which(jmat.ref > 5, arr.ind = TRUE)] <- 5
jmat.ref[which(jmat.ref < -5, arr.ind = TRUE)] <- -5

# plot(density(unlist(de.out[[mark.ref]]$effects.mat[rows.keep, ])))
plot(density(unlist(jmat.ref)))

jmat.ref <- t(scale(t(jmat.ref), center = TRUE, scale = TRUE))

jgrep <- "Ebf1"
pval.mat.merge[grep(pattern = jgrep, rownames(pval.mat.merge), value = TRUE),] 
rname.test <- grep(pattern = jgrep, rownames(jmat.ref), value = TRUE)

plot(jmat.ref[rname.test, ])
text(jmat.ref[rname.test, ], labels = colnames(jmat.ref))


# sort by granuspecific
cname <- "H3K4me1_ClusterEryths"
cname <- "H3K27me3_ClusterBcells"
cname <- "H3K4me1_ClusterBcells"

jmat.ref <- jmat.ref[order(jmat.ref[, cname], decreasing = TRUE), ]
# jmat.ref <- jmat.ref[order(jmat.ref[, cname], decreasing = FALSE), ]

ctype.rnames <- rownames(jmat.ref)[1:100]

print(head(ctype.rnames))

# check geness in k27me3
jmat.check <- t(scale(t(de.out$H3K27me3$effects.mat[ctype.rnames, ]), center = TRUE, scale = TRUE))
jmat.check2 <- t(scale(t(de.out$H3K4me3$effects.mat[ctype.rnames, ]), center = TRUE, scale = TRUE))



boxplot(jmat.check, las = 2)
boxplot(jmat.check2, las = 2)
boxplot(jmat.ref[ctype.rnames, ], las = 2)



# Load celltype specific genes from rnaseq  -------------------------------


infrdata <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression/integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_TRUE.forPoissonRegression.CountR1only.2020-06-05.smaller.RData"
assertthat::assert_that(file.exists(infrdata))
load(infrdata, v=T)

infrdata2 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData"
load(infrdata2, v=T)



library(hash)
e2g <- hash::invert(g2e)

# gset <- sapply(de.ens.sorted.stringent$Erythroblast, AssignHash, e2g)
eset <- de.ens.sorted.stringent$Erythroblast
eset <- de.ens.sorted.stringent$Bcell
eset <- de.ens.sorted.stringent$HSCs
eset <- de.ens.sorted.stringent$Neutrophil
# filter rnames
rnames.genes <- sapply(sapply(rownames(jmat.ref), function(x) strsplit(x, "\\.")[[1]][[4]]), function(x) strsplit(x, "_")[[1]][[1]])

rnames.ens <- sapply(rnames.genes, function(x) AssignHash(x, g2e, null.fill = NA))

rnames.ens.filt <- rnames.ens[rnames.ens %in% eset]

jsums <- rowSums(de.out$H3K27me3$jmat.mark)

# filter only rownames with some counts for K27me3
enough.counts <- names(jsums[which(jsums > 1000)])

rnames.filt1 <- names(rnames.ens.filt)
rnames.filt2 <- rnames.filt1[which(rnames.filt1 %in% enough.counts)]

# jmat.ref.filt <- jmat.ref[names(rnames.ens.filt), ]
# jmat.check.de <- t(scale(t(de.out$H3K4me1$effects.mat[names(rnames.ens.filt), ]), center = TRUE, scale = TRUE))
# jmat.check2.de <- t(scale(t(de.out$H3K4me3$effects.mat[names(rnames.ens.filt), ]), center = TRUE, scale = TRUE))

jmat.ref.filt <- jmat.ref[rnames.filt2, ]
jmat.check.de <- t(scale(t(de.out$H3K4me1$effects.mat[rnames.filt2, ]), center = TRUE, scale = TRUE))
jmat.check2.de <- t(scale(t(de.out$H3K4me3$effects.mat[rnames.filt2, ]), center = TRUE, scale = TRUE))

# plot(density(-log10(unlist(pval.mat.merge[names(rnames.ens.filt), 1]))))
# plot(density(-log10(unlist(pval.mat.merge[names(rnames.ens.filt), 2]))))
# plot(density(-log10(unlist(pval.mat.merge[names(rnames.ens.filt), 3]))))

boxplot(jmat.ref.filt, las = 2)
boxplot(jmat.check.de, las = 2)
boxplot(jmat.check2.de, las = 2)

jmat.ref.filt[order(jmat.ref.filt[, "H3K27me3_ClusterHSPCs"]), ][1:5, 1:5]

# check total

jrname <- "15:9506158-9529876;NM_008372.4..Il7r_-"
pval.mat.merge[jrname, ]

plot(density(de.out$H3K4me1$jmat.mark[jrname, ]))
plot(density(de.out$H3K4me3$jmat.mark[jrname, ]))
plot(density(de.out$H3K27me3$jmat.mark[jrname, ]))


