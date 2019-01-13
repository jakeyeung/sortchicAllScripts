# Jake Yeung
# run_GREAT_downstream.R
# Analyze LDA downstream using GREAT 
# 2019-01-06

jstart <- Sys.time()

library(topicmodels)
library(dplyr)
library(ggplot2)
library(umap)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rGREAT)
library(JFuncs)

source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")

# ARGS

args <- commandArgs(trailingOnly=TRUE)

inf <- args[[1]]
outpath <- args[[2]]
ncores <- StrToNumeric(args[[3]])
top.thres <- StrToNumeric(args[[4]])
# top.thres <- 0.98  is sensible?

load(inf, v=T)

# Model selection

out.lda.lst <- out.lda

Kvec <- sapply(out.lda.lst, function(x) x@k)

best.K <- Kvec[which.max(sapply(out.lda.lst, function(x) x@loglikelihood))]

# plot loglikelihood: why does it decrease eventually??
# pdf(outdir, "best_K.pdf")
#   par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
#   plot(Kvec, sapply(out.lda.lst, function(x) x@loglikelihood))
# dev.off()

print(paste("Best K:", best.K))

kchoose <- best.K
out.lda <- out.lda.lst[[which(Kvec == kchoose)]]
tmResult <- posterior(out.lda)

# Run GREAT --------------------------------------------------------------

# need to assign cutoff for each peak for each topic
topic.regions <- lapply(seq(best.K), function(clst){
  return(SelectTopRegions(tmResult$terms[clst, ], colnames(tmResult$terms), method = "thres", method.val = top.thres))
})
regions <- data.frame(seqnames = sapply(colnames(tmResult$terms), GetChromo),
                      start = sapply(colnames(tmResult$terms), GetStart),
                      end = sapply(colnames(tmResult$terms), GetEnd))
regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                annoDb='org.Mm.eg.db'))
# rownames(regions.annotated) <- regions.annotated$region_coord
regions.annotated$region_coord <- names(regions.range)

print(paste("Running great multicore", ncores))
out.great.lst <- mclapply(seq(best.K), function(i){
  gr.in <- regions.range[topic.regions[[i]], ]
  out.great <- submitGreatJob(gr.in, species="mm10", request_interval = 30)
  return(out.great)
}, mc.cores = ncores)
out.tb.lst <- mclapply(out.great.lst, function(out.great){
  out.tb <- getEnrichmentTables(out.great, ontology=availableOntologies(out.great), 
                                request_interval = 30)
  return(out.tb)
}, mc.cores = ncores)

save(topic.regions, regions.annotated, out.great.lst, out.tb.lst, out.lda, file = outpath)

print(Sys.time() - jstart)
print("Done successfully")
