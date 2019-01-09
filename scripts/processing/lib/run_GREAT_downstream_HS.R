# Jake Yeung
# run_GREAT_downstream.R
# Analyze LDA downstream using GREAT 
# use Human Genome
# 2019-01-09

jstart <- Sys.time()

library(topicmodels)
library(dplyr)
library(ggplot2)
library(umap)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

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
top.thres <- 0.98

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
                                                TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                                annoDb='org.Hs.eg.db'))
rownames(regions.annotated) <- regions.annotated$region_coord

# submit fake greatjob to initialize the database?


InitGreatTest <- function(){
    great.test = tryCatch({
        i <- 1
        gr.in <- regions.range[topic.regions[[i]], ]
        out.great.tmp <- submitGreatJob(gr.in, species="hg19")
    }, error = function(e) {
        print("Error on submit Great job, dont know why")
    })
}

if (ncores == 1){
  print(paste("Running great single core", ncores))
  out.great.lst <- lapply(seq(best.K), function(i){
    tmp <- InitGreatTest()
    gr.in <- regions.range[topic.regions[[i]], ]
    out.great <- submitGreatJob(gr.in, species="hg19", request_interval = 700)
    return(out.great)
  })
  out.tb.lst <- lapply(out.great.lst, function(out.great){
    out.tb <- getEnrichmentTables(out.great, ontology=availableOntologies(out.great), 
                                  request_interval = 300)
    return(out.tb)
  })
} else {
  print(paste("Running great multicore", ncores))
  out.great.lst <- mclapply(seq(best.K), function(i){
    tmp <- InitGreatTest()
    gr.in <- regions.range[topic.regions[[i]], ]
    out.great <- submitGreatJob(gr.in, species="hg19", request_interval = 700)
    return(out.great)
  }, mc.cores = ncores)
  out.tb.lst <- mclapply(out.great.lst, function(out.great){
    out.tb <- getEnrichmentTables(out.great, ontology=availableOntologies(out.great), 
                                  request_interval = 300)
    return(out.tb)
  }, mc.cores = ncores)
}

save(topic.regions, regions.range, out.great.lst, out.tb.lst, out.lda, file = outpath)

print(Sys.time() - jstart)
print("Done successfully")
