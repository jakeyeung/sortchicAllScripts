# Jake Yeung
# Date of Creation: 2019-06-25
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/explore_count_matrix/explore_TSS_count_matrix.R
# Explore TSS count matrix

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(DESeq2)

library(svs)
require(statmod)

library(fastICA)


# Load  -------------------------------------------------------------------


inf.umap <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6/dat_umap_long_with_louvain.H3K4me1.RData"
load(inf.umap, v=T)


inf <- "/Users/yeung/data/scchic/from_cluster/count_mat_TSS_B6/B6-BM-H3K4me1.merged.NoCountThres.GeneTSS.Dedup.Robj"
assertthat::assert_that(file.exists(inf))

load(inf, v=T)

# lib.size <- estimateSizeFactorsForMatrix(as.matrix(count.dat$counts))
lib.size <- Matrix::colSums(count.dat$counts)
ed <- t(t(count.dat$counts)/lib.size)

# take top genes?
genes.mean <- Matrix::rowMeans(count.dat$counts)
genes.var <- apply(count.dat$counts, 1, var)
genes.cv2 <- genes.var / genes.mean ^ 2

genes.mean <- Matrix::rowMeans(ed)
genes.var <- apply(ed, 1, var)
genes.cv2 <- genes.var / genes.mean ^ 2

par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9)
smoothScatter(log(genes.mean),log(genes.cv2))

jlong <- data.frame(genes.mean.log10 = log10(genes.mean), genes.cv2.log10 = log10(genes.cv2), cell = unname(names(genes.mean)))

ggplot(jlong, aes(x = genes.mean.log10, y = genes.cv2.log10)) + geom_point(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(slope = -1, intercept = -3.5)


# Fit curve ---------------------------------------------------------------

minMeanForFit <- unname( quantile( genes.mean[ which( genes.cv2 > .3 ) ], .95 ) )
useForFit <- genes.mean >= minMeanForFit # & spikeins
X <- cbind( a0 = 1, a1tilde = 1/genes.mean[useForFit] )
Y <- genes.cv2[useForFit]
fit <- glmgam.fit(X, Y)
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
fit$coefficients


# Plot output -------------------------------------------------------------

# repeat previous plot
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(genes.mean),log(genes.cv2));
xg <- exp(seq( min(log(genes.mean[genes.mean>0])), max(log(genes.mean)), length.out=1000 ))
vfit <- a1/xg + a0
# add fit line
lines( log(xg), log(vfit), col="black", lwd=3 )
df <- ncol(ed) - 1
# add confidence interval
lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")


# Rank genes --------------------------------------------------------------

afit <- a1/genes.mean+a0
varFitRatio <- genes.var/(afit*genes.mean^2)
varorder <- order(varFitRatio,decreasing=T)
oed <- ed[varorder,]

# repeat previous plot
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(genes.mean),log(genes.cv2)); lines( log(xg), log(vfit), col="black", lwd=3 ); lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black"); lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black");
# add top 100 genes
points(log(genes.mean[varorder[1:100]]),log(genes.cv2[varorder[1:100]]),col=2)


pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
adj.pval <- p.adjust(pval,"fdr")
sigVariedGenes <- adj.pval<1e-3;
table(sigVariedGenes)


# Cluster by top genes ----------------------------------------------------

top.genes <- names(sigVariedGenes[which(sigVariedGenes)])

ed.filt <- sqrt(as.matrix(t(ed[top.genes, ])))
# ed.filt <- log2(as.matrix(t(ed[top.genes, ])) * 10^6 + 1)
# ed.filt <- t(log2(as.matrix(t(ed.filt)) * 10^6 + 1))
# ed.filt <- t(sqrt(as.matrix(t(ed.filt))))
pca.out <- prcomp(ed.filt, center = TRUE, scale. = TRUE)

dat.proj <- as.matrix(ed.filt) %*% pca.out$rotation

plot(dat.proj[, 1], dat.proj[, 2])

# System time -------------------------------------------------------------

# mat <- as.matrix(t(count.dat$counts[top.genes, ]))
mat <- as.matrix(t(ed[top.genes, ]))
mat[which(mat > 1)] <- 1
# samps <- sample(seq(ncol(mat)), size = 700)
# mat <- mat[, samps]

# system.time(
#   out <- fast_plsa(mat, k = 50)
# )

mat <- log2(mat * 10^6 + 1)

out.pca <- prcomp(mat, center = TRUE, scale. = TRUE)
out.ica <- fastICA(mat, n.comp=10)

plot(out.pca$x[, 1], out.pca$x[, 2], pch = 20)
plot(out.pca$x[, 2], out.pca$x[, 3], pch = 20)

plot(out.ica$S[, 1], out.ica$S[, 2], pch = 20)

mat.proj <- mat %*% diag(out.pca$scale) %*% out.pca$rotation
plot(mat.proj[, 1], mat.proj[, 2])

# color by louvain?

assertthat::assert_that(all(rownmaes(out.ica$S) == rownames(mat.proj)))
jmark <- "H3K4me1"
pca.long <- data.frame(mat.proj[, 1:5], out.ica$S[, 1:5], cell.raw = names(rownames(mat.proj)), stringsAsFactors = FALSE) %>%
  rowwise()  %>%
  mutate(repcell = strsplit(cell.raw, "\\.")[[1]][[22]],
         cell = paste0("B6-13W1-BM-", jmark, "-", repcell)) %>%
  dplyr::select(-cell.raw)




pca.long.merged <- left_join(pca.long, dat.umap.long)

# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "#D3D3D3")
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "#D3D3D3")
m1 <- ggplot(pca.long.merged, aes(x = PC1, y = PC2, color = louvain)) + geom_point(size = 3) + theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)
m2 <- ggplot(pca.long.merged, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) + theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)
multiplot(m1, m2, cols = 2)



ggplot(pca.long.merged, aes(x = X1, y = X2, color = louvain)) + geom_point() + theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

ggplot(pca.long.merged, aes(x = PC2, y = PC3, color = louvain)) + geom_point() + theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

ggplot(pca.long.merged, aes(x = PC3, y = PC4, color = louvain)) + geom_point() + theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

ggplot(pca.long.merged, aes(x = PC3, y = PC4, color = louvain)) + geom_point() + theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

# Plot output -------------------------------------------------------------






