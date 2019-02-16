# Jake Yeung
# 1b-prepare_count_mat_but_dont_filter.R


jstart <- Sys.time() 

library(data.table)
library(Matrix)

# args <- commandArgs(trailingOnly=TRUE)
# inf <- args[[1]]
# metaobj <- args[[2]]
# outf <- args[[3]]

# jchip <- "H3K4me3"
# jchip <- "H3K27me3"
# jchip <- "H3K9me3"
jchips <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
# jchips <- c("H3K4me3")

jdist <- 100
for (jchip in jchips){

inf <- paste0("/hpc/hub_oudenaarden/avo/scChiC/metacell/K562-G1-", jchip, "-", jdist, "kb.txt")
outdir <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/count_mat_K562/count_mats_binned"

dir.create(outdir)

dat <- read.table(inf)

# keep all cells and peaks, filter them before running LDA
good.cells <- colnames(dat)
good.peaks <- rownames(dat)

rows.i <- which(rownames(dat) %in% good.peaks)
cols.i <- which(colnames(dat) %in% good.cells)

print(paste("Keeping ", length(rows.i), "rows"))
print(paste("Keeping  ", length(cols.i), "columns"))

dat.sub <- Matrix(as.matrix(dat[rows.i, cols.i]), sparse = TRUE)

# rename rownames to chr1:start-end
rnames <- rownames(dat.sub)

# add chr
rnames <- sapply(rnames, function(x) paste("chr", x, sep=""))
# rnames <- paste("chr", rnames, sep="")
print(head(rnames))
rnames.new <- sapply(rnames, function(x) paste(paste(strsplit(x, "_")[[1]][[1]], strsplit(x, "_")[[1]][[2]], sep=":"), strsplit(x, "_")[[1]][[3]], sep="-"))
print(head(rnames.new))

rownames(dat.sub) <- rnames.new

count.dat <- list(counts = dat.sub)  # make object for LDA

save(count.dat, file = file.path(outdir, paste0("K5562-G1-", jchip, ".no_filt.Robj")))

print(jstart - Sys.time())

}

