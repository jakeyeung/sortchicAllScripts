# Jake Yeung
# make_count_matrix_from_bams.R
# Make count matrix from bams 
# 2018-12-19

jstart <- Sys.time() 

print(paste("In workdir:", getwd()))

library(Rsubread)
library(Matrix)
library(stringr)
library(hash)


source("scripts/Rfunctions/GetMetaData.R")

GetValueFromHash <- function(key, jhash){
  val <- jhash[[key]]
  if (is.null(val)){
    warning(paste("Key not found, defaulting to key:", key))
    val <- key
  }
  return(val)
}

args <- commandArgs(trailingOnly=TRUE)

bampaths <- args[[1]]  # file list of bam paths
regions <- args[[2]]  # file of regions to count reads
outpath <- args[[3]]  # sparse count matrix. R object?

outdir <- dirname(outpath)

jpaired <- TRUE  # probably always paired End

print(paste("Paired End mode: ", jpaired))

bamfiles <- unlist(read.csv2(bampaths, stringsAsFactors=FALSE), use.names=FALSE)

print(bamfiles)

regions.dat <- read.table(regions, header=FALSE, sep = "\t", stringsAsFactors=FALSE)


# make required columns for Rsubread
if (ncol(regions.dat) == 9){
  colnames(regions.dat) <- c("Chr", "Start", "End", "pname", "score", "strand", "signal", "mlogpval", "mlogqval")
} else if (ncol(regionsdat) == 15){
  colnames(regions.dat) <- c("Chr", "Start", "End", "pname", "score", "strand", "tStart", "tEnd", "rgb", "bCount", "bSizes", "bStarts", "signal", "mlogpval", "mlogqval")
} else {
  stop("Number of columns must be 9 or 15")
}
# GeneID as chromo:start-end
GeneID <- paste(regions.dat$Chr, ":", regions.dat$Start, "-", regions.dat$End, sep='')
regions.dat$GeneID <- GeneID

print("Regions file looks like this...")
print(head(regions.dat))

# other interesting features: countMultiMappingReads
count.dat <- Rsubread::featureCounts(bamfiles, annot.ext=regions.dat, read2pos=5, isPairedEnd=jpaired, reportReadsPath=outdir)

# make sparse
count.dat$counts <- Matrix(count.dat$counts, sparse=TRUE)  # save some space?

cnames <- sapply(bamfiles, function(fullpath){
  # make full path
  bname <- basename(fullpath)
  # 8 letter UMI
  bc <- str_match(fullpath, pattern = "[ACGT]{8}")
  tissue <- GetTissue(fullpath)
  chip <- GetChip(fullpath)
  biorep <- GetBioRep(fullpath)
  techrep <- GetTechRep(fullpath)
  experi <- GetExperiment(fullpath)
  cname.new <- paste(chip, tissue, biorep, techrep, experi, bc, sep="-")
  return(cname.new)
})

# make hash to match old targets to new names
name.hash <- hash(count.dat$targets, cnames)

count.dat$targets <- sapply(count.dat$targets, GetValueFromHash, name.hash)
colnames(count.dat$counts) <- sapply(colnames(count.dat$counts), GetValueFromHash, name.hash)
colnames(count.dat$stat) <- sapply(colnames(count.dat$stat), GetValueFromHash, name.hash)

save(count.dat, file = outpath)

print(Sys.time() - jstart)

