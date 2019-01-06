# Jake Yeung
# Date of Creation: 2019-01-03
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/lda_binarized.R
# Run LDA on binarized data


rm(list=ls())

library(cisTopic)


# Helper functions --------------------------------------------------------



# Load data ---------------------------------------------------------------

load("/private/tmp/lda_output/PZ-BM-H3K27me3.merged.NoCountThres.Robj", v=T)
load("/private/tmp/lda_outputs.meanfilt_1.merge_1000_NoM_cisTopic.cellmin_1000.cellmax_50000/H3K27me3_NoCountThres_1000_noM_cisTopicOutput.Robj", v=T)

# remove chromoM?
count.mat <- count.dat$counts
M.peaks <- grep("chrM", rownames(count.mat), value=TRUE)
count.mat <- count.mat[which(!rownames(count.mat) %in% M.peaks), ]

# Run analysis ------------------------------------------------------------

# cisTopicObject <- createcisTopicObject(count.mat, project.name='H3K27me_BM')

# takes some time
# cisTopicObject <- runModels(cisTopicObject, topic=c(10, 15, 20, 25), seed=0, nCores=4, burnin = 120, iterations = 150, addModels=FALSE)
# save(cisTopicObject, file = "outputs_R/lda_output/cistopic_PZ-BM-H3K27me3.merged.NoCountThres.Robj")

cisTopicObject <- selectModel(cisTopicObject)

nn <- 5
jmetric <- 'pearson'
jmindist <- 0.0001
custom.settings <- umap.defaults
custom.settings$n_neighbors <- nn
custom.settings$metric <- jmetric
custom.settings$min_dist <- jmindist

modelMat <- .modelMatSelection(cisTopicObject, target="cell", method="Probability")
dat.umap.cistopic <- umap(t(modelMat), config = custom.settings)
# color by loadings on Kvec
jcol <- modelMat[1, ]
colorPal <- grDevices::colorRampPalette(c("pink", "red", "darkred"))
jcol.rgb <- colorPal(20)[as.numeric(cut(jcol,breaks = 20))]
par(mfrow=c(1,1), mar=c(1,1,1,1))
plot(dat.umap.cistopic$layout[, 1], dat.umap.cistopic$layout[, 2], pch = 20, main = jmain, col = jcol.rgb)
# cisTopicObject <- runPCA(cisTopicObject, target='cell', method = "Z-score")
# cisTopicObject <- runtSNE(cisTopicObject, target='cell', method = "Z-score")


# cisTopicObject <- runUmap(cisTopicObject, target='cell', method = "Probability", config = custom.settings)
# par(mfrow=c(4,5), mar=c(1,1,1,1))
# plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Probability', topics = 1:20,
#              colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)

# par(mfrow=c(2,5))
# plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Z-score', topics = 11:20,
#              colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)

# plot across different 

# plotFeatures(cisTopicObject, method='PCA', target='cell', topic_contr='Z-score', 
#              colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=3, legend=TRUE)

pred.matrix <- predictiveDistribution(cisTopicObject)


# Downstream of peaks -----------------------------------------------------

cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)

system.time(
  cisTopicObject <- runUmap(cisTopicObject, target='region')  # super slow?
)
par(mfrow=c(4,5), mar=c(1,1,1,1))
plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Z-score', topics = 1:20,
             colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Z-score', topics = 1:20,
             colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)

# Link regions to genes ---------------------------------------------------

cisTopicObject <- annotateRegions(cisTopicObject, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb='org.Mm.eg.db')

# par(mfrow=c(1,1))
# signaturesHeatmap(cisTopicObject, selected.signatures = 'annotation')
plotFeatures(cisTopicObject, method='PCA', target='region', topic_contr=NULL, colorBy=c('annotation'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, intervals=20)

system.time(
  cisTopicObject <- GREAT(cisTopicObject, genome='mm10', fold_enrichment=1.5, geneHits=1, sign=0.1, request_interval=10)
)
# save(cisTopicObject, file = "outputs_R/lda_output/cistopic_GREAT_PZ-BM-H3K27me3.merged.NoCountThres.Robj")


ontos <- unique(unlist(sapply(cisTopicObject@binarized.rGREAT, names)))
jonto <- "GO Cellular Component"
jonto <- "MSigDB Immunologic Signatures"
jonto <- "GO Biological Process"
jonto <- "MSigDB Pathway"
jonto <- "MSigDB Predicted Promoter Motifs"
jonto <- "MSigDB Perturbation"
jonto <- "MSigDB Immunologic Signatures"
jonto <- "InterPro"

topics.keep <- sapply(cisTopicObject@binarized.rGREAT, function(x) is.null(x[[jonto]]))
topics.keep.i <- which(!topics.keep)
ontologyDotPlot(cisTopicObject, top=5, topics=topics.keep.i, ontology = jonto, var.y='name', order.by='Binom_Adjp_BH')


# Zscore ------------------------------------------------------------------

nn.vec <- c(5, 10, 20, 40, 60)
jmindist.vec <- c(0.0001, 0.01, 0.1, 0.25, 0.5)

interaction(nn.vec, jmindist.vec, sep = "_")

# (cisTopicObject, top=5, topics=topics.keep.i, var.y='name', order.by='Binom_Adjp_BH')

# cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.999, plot=TRUE)



# 
# par(mfrow=c(2,5))
# plotFeatures(cisTopicObject, method='PCA', target='region', topic_contr="Z-score", topics = seq(10),
#              colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
# 
# par(mfrow=c(2,5))
# plotFeatures(cisTopicObject, method='Umap', target='region', topics = seq(10),
#              colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
# par(mfrow=c(2,5))
# plotFeatures(cisTopicObject, method='Umap', target='region', topics = 11:20,
#              colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)

