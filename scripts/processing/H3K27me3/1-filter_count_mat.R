# Jake Yeung
# 1-filter_count_mat.R
# Load metacell output from Dropbox, filter out count mat so LDA inputs same as metacell
# 2018-12-28

jstart <- Sys.time() 

library(data.table)
library(Matrix)

# args <- commandArgs(trailingOnly=TRUE)
# inf <- args[[1]]
# metaobj <- args[[2]]
# outf <- args[[3]]

inf <- "/hpc/hub_oudenaarden/avo/scChiC/metacell/BM-H3K27me3.txt"
metaobj <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/from_dropbox/FIG4_BM/H3K27me3.datadir_mc_f.Rda"
outf <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_analysis.h3k27me3.dropbox"

dir.create(outf)

# dat <- fread(inf, header=TRUE)
dat <- read.table(inf)

load(metaobj, v=T)  # loads object

# filter data table based on metaobj @mc slot

good.cells <- names(object@mc)
good.peaks <- rownames(object@mc_fp)

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

save(count.dat, file = file.path(outf, "BM-H3K27me3.AvO_filt.Robj"))

print(jstart - Sys.time())
