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

if (length(args) == 4){
  use.pname = TRUE
} else {
  use.pname = FALSE
}

assertthat::assert_that(file.exists(bampaths))
assertthat::assert_that(file.exists(regions))

outdir <- dirname(outpath)

jpaired <- TRUE  # probably always paired End

print(paste("Paired End mode: ", jpaired))

bamfiles <- unlist(read.csv2(bampaths, stringsAsFactors=FALSE), use.names=FALSE)

print(bamfiles)

regions.dat <- read.table(regions, header=FALSE, sep = "\t", stringsAsFactors=FALSE)


# make required columns for Rsubread
if (ncol(regions.dat) == 9){
  colnames(regions.dat) <- c("Chr", "Start", "End", "pname", "score", "strand", "signal", "mlogpval", "mlogqval")
} else if (ncol(regions.dat) == 15){
  colnames(regions.dat) <- c("Chr", "Start", "End", "pname", "score", "strand", "tStart", "tEnd", "rgb", "bCount", "bSizes", "bStarts", "signal", "mlogpval", "mlogqval")
} else if (ncol(regions.dat) == 4){
  # handle hidden domains output
  # add strand column, default "+"
  regions.dat$strand <- "+"
  colnames(regions.dat) <- c("Chr", "Start", "End", "pname", "strand")
} else {
  stop(paste("Number of columns must be 9 or 15, found", ncol(regions.dat)))
}
# GeneID as chromo:start-end
if (!use.pname){
  GeneID <- paste(regions.dat$Chr, ":", regions.dat$Start, "-", regions.dat$End, sep='')
} else {
  GeneID <- regions.dat$pname
}
regions.dat$GeneID <- GeneID

print("Regions file looks like this...")
print(head(regions.dat))

# other interesting features: countMultiMappingReads
count.dat <- Rsubread::featureCounts(bamfiles, annot.ext=regions.dat, read2pos=5, isPairedEnd=jpaired, reportReadsPath=outdir, allowMultiOverlap=TRUE)  # allowMultiOverlap otherwise you get zero reads sometimes if you convert TSS to counts

# make sparse
count.dat$counts <- Matrix(count.dat$counts, sparse=TRUE)  # save some space?

cnames <- sapply(bamfiles, function(fullpath){
  # make full path
  bname <- basename(fullpath)
  # 8 letter UMI
  bc <- str_match(bname, pattern = "[ACGT]{8}")
  tissue <- GetTissue(bname)
  chip <- GetChip(bname)
  biorep <- GetBioRep(bname)
  techrep <- GetTechRep(bname)
  experi <- GetExperiment(bname)
  cname.new <- paste(chip, tissue, biorep, techrep, experi, bc, sep="-")
  return(cname.new)
})

# make hash to match old targets to new names
name.hash <- hash(count.dat$targets, cnames)

count.dat$targets <- sapply(count.dat$targets, GetValueFromHash, name.hash)
colnames(count.dat$counts) <- sapply(colnames(count.dat$counts), GetValueFromHash, name.hash)
colnames(count.dat$stat) <- sapply(colnames(count.dat$stat), GetValueFromHash, name.hash)

save(count.dat, file = outpath)
# write count matrix as full textfile
# outpath.noext <- sub("^([^.]*).*", "\\1", outpath)
# write.table(as.matrix(count.dat$mat), file = matpath, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

print(Sys.time() - jstart)

