# Jake Yeung
# cisTopic_downstream.R
# Downstream processing of cistopic
# 2019-01-03

library(umap)
library(cisTopic)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(parallel)

args <- commandArgs(trailingOnly=TRUE)

GetUmapSettings <- function(nn=5, jmetric='pearson', jmindist=0.0001){
  custom.settings <- umap.defaults
  custom.settings$n_neighbors <- nn
  custom.settings$metric <- jmetric
  custom.settings$min_dist <- jmindist
  return(custom.settings)
}

inf <- args[[1]]
outdir <- args[[2]]

load(inf, v=T)  # cisTopicObject

pdf(file.path(outdir, "best_model_and_binarize.pdf"))
cisTopicObject <- selectModel(cisTopicObject)
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)
dev.off()

# binarize


# plot many diffrent Umaps and plot into a file
nn.vec <- c(5, 10, 20, 40, 60)
# jmet.vec <- c('pearson', 'euclidean')
jmet.vec <- 'euclidean'
jmindist.vec <- c(0.0001, 0.01, 0.1, 0.25, 0.5)
jdistmeth <- "Z-score"
vec.lst <- levels(as.factor(nn.vec):as.factor(jmindist.vec))  # get combinations of nn.vec and jmindist.vec, length 5x5

jmet <- 'euclidean'
out.lst <- lapply(vec.lst, function(nn.jmindist){
  nn <- as.numeric(strsplit(nn.jmindist, ":")[[1]][[1]])
  jmindist <- as.numeric(strsplit(nn.jmindist, ":")[[1]][[2]])
  custom.settings <- GetUmapSettings(nn, jmet, jmindist)
  out.tmp <- runUmap(cisTopicObject, target='cell', method = jdistmeth, config = custom.settings)
  out.tmp <- runUmap(cisTopicObject, target='region', method = jdistmeth, config = custom.settings)
  return(out.tmp)
# }, mc.cores = length(vec.lst))  # 25 cores?
})

save(out.lst, file = file.path(outdir, paste0("out_list_umap_cells_regions", jdistmeth, ".", jmet, ".pdf")))

# # run regions plot
# out.lst.regions <- mclapply(vec.lst, function(nn.jmindist){
#   nn <- as.numeric(strsplit(nn.jmindist, ":")[[1]][[1]])
#   jmindist <- as.numeric(strsplit(nn.jmindist, ":")[[1]][[2]])
#   custom.settings <- GetUmapSettings(nn, jmet, jmindist)
#   out.tmp <- runUmap(cisTopicObject, target='region', method = jdistmeth, config = custom.settings)
#   return(out.tmp)
# }, mc.cores = length(vec.lst))  # 25 cores?
# save(out.lst.regions, file = file.path(outdir, paste0("out_list_umap_regions", jdistmeth, ".", jmet, ".pdf")))

# plot outputs
outplot <- file.path(outdir, paste0("umap_across_variables.", jdistmeth, ".", jmet, ".pdf"))
pdf(outplot, useDingbats=FALSE)
par(mfrow=c(1,1), mar=c(1,1,1,1))
for (i in 1:length(out.lst)){
    x <- out.lst[[i]]
    nn.jmindist <- vec.lst[[i]]
    nn <- as.numeric(strsplit(nn.jmindist, ":")[[1]][[1]])
    jmindist <- as.numeric(strsplit(nn.jmindist, ":")[[1]][[2]])
    jmain <- paste("NN:", nn, "MinDist", jmindist)
    plotFeatures(x, method='Umap', target='cell', topic_contr=jmet, topics = 1, colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
    plotFeatures(x, method='Umap', target='region', topic_contr=jmet, topics = 1, colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
    title(main = jmain, line = -1)
}
dev.off()

# run GREAT
cisTopicObject <- annotateRegions(cisTopicObject, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb='org.Mm.eg.db')
cisTopicObject.great <- GREAT(cisTopicObject, genome='mm10', fold_enrichment=1.5, geneHits=1, sign=0.1, request_interval=10)

save(cisTopicObject.great, file = file.path(outdir, paste0("cistopic_great_output.pdf")))

