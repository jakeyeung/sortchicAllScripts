# Jake Yeung
# integrate_datasets.R
# 2019-02-22
# DESCRIPTION
# 
#     Integrate datasets using Seurat
# 
# FOR HELP
# 
#     Rscript integrate_datasets.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-02-22
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

jstart <- Sys.time() 

setwd("~/projects/scChiC")

suppressPackageStartupMessages(library("argparse"))

library(GGally)
library(purrr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(umap)
library(data.table)
library(dplyr)
library(hash)
library(JFuncs)
library(topicmodels)
library(scales)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(cowplot)

library(Seurat)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

# parser$add_argument('infile', metavar='INFILE',
#                                             help='<+INFILE+>')
# parser$add_argument('outfile', metavar='OUTFILE',
#                                             help='<+OUTFILE+>')
# parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
#                         help="Print extra output [default]")
#                                         
# # get command line options, if help option encountered print help and exit,
# # otherwise if options not found on command line then set defaults, 
# args <- parser$parse_args()
# 
# # print some progress messages to stderr if "quietly" wasn't requested
# if ( args$verbose ) { 
#     print("Arguments:")
#     print(args)
# }
# 

# Run interactively on H3K4me1 and H3K4me3

# load LDA
# dirmain <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_100.cellmax_500000.binarize.FALSE.no_filt"
dirmain <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell"

# jmarks <- c("H3K4me1", "H3K4me3")
# jmarks <- c("H3K27me3", "H3K9me3")
# jmarks <- c("H3K4me1", "H3K27me3")
jmarks <- c("H3K4me1", "H3K9me3")
jmarks.str <- paste(jmarks, collapse = "_")

infs <- list(file.path(dirmain, "lda_outputs.meanfilt_10.cellmin_100.cellmax_500000.binarize.TRUE.no_filt", paste0("lda_out_meanfilt.BM-", jmarks[[1]], ".CountThres0.K-15_20_25_30_35.Robj")),
             file.path(dirmain, "lda_outputs.meanfilt_10.cellmin_100.cellmax_500000.binarize.FALSE.no_filt", paste0("lda_out_meanfilt.BM-", jmarks[[2]], ".CountThres0.K-15_20_25_30.Robj")))
             
print(infs)
out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks, infs, SIMPLIFY = FALSE)

names(out.objs) <- jmarks

tm.result.lst <- lapply(out.objs, function(x) posterior(x[['out.lda']]))
mat.impute.lst <- lapply(tm.result.lst, function(tm.result) t(tm.result$topic %*% tm.result$term))

print(lapply(mat.impute.lst, dim))


# Wrangle data ------------------------------------------------------------

top.thres <- 0.999
top.regions.lst <- lapply(jmarks, function(jmark){
  topic.regions <- out.objs[[jmark]]$topic.regions  # for each cluster
  top.regions <- unique(unlist(topic.regions))
  return(top.regions) 
})
top.regions.merge <- unique(unlist(top.regions.lst))

mat.merged.lst <- lapply(mat.impute.lst, function(mat.impute){
  row.i <- which(top.regions.merge %in% rownames(mat.impute))
  return(as.data.frame(mat.impute[row.i, ]))
})

rnames.lst <- lapply(mat.merged.lst, function(x) rownames(x))
top.regions.intersect <- purrr::reduce(.x = rnames.lst, .f = intersect)

mat.merged.lst <- lapply(mat.merged.lst, function(mat.merged){
  return(mat.merged[top.regions.intersect, ])
})

mat.merged <- bind_cols(mat.merged.lst[[jmarks[[1]]]], mat.merged.lst[[jmarks[[2]]]])
rownames(mat.merged) <- top.regions.intersect

coldat <- data.frame(tech = sapply(colnames(mat.merged), function(x) strsplit(x, "_")[[1]][[2]]))
rownames(coldat) <- colnames(mat.merged)

bm <- CreateSeuratObject(counts = mat.merged, meta.data = coldat)
bm.lst <- SplitObject(object = bm, split.by = "tech")

# get variable genes
for (i in 1:length(x = bm.lst)) {
  bm.lst[[i]] <- FindVariableFeatures(object = bm.lst[[i]], selection.method = "dispersion", nfeatures = 500, verbose = TRUE)
}

# sort(jmarks.all)
reference.list <- bm.lst[sort(unique(coldat$tech))]

bm.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
bm.integrated <- IntegrateData(anchorset = bm.anchors, dims = 1:30)

save(mat.merged, coldat, bm.lst, bm.anchors, bm.integrated, file = paste0("/hpc/hub_oudenaarden/jyeung/data/scChiC/integrated_datasets/integrated_", jmarks.str, "_seurat_again.RData"))

saveRDS(bm.integrated, file = paste0("/hpc/hub_oudenaarden/jyeung/data/scChiC/integrated_datasets/bm_integrated_", jmarks.str, "_seurat.rds"))

print(Sys.time() - jstart)
