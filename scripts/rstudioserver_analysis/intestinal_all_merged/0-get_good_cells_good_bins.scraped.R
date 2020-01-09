# Jake Yeung
# Date of Creation: 2019-12-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/intestinal_all_merged/0-get_good_cells_good_bins.scraped.R
# 


rm(list=ls())

jstart <- Sys.time()

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)


library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(GGally)


# Load data ---------------------------------------------------------------

inf.meta <- "/home/jyeung/hpc/intestinal_scchic/metadata_all/OUD_to_filenames.txt"
meta <- fread(inf.meta, header = FALSE, col.names = c("Run", "Fastq")) %>%
  rowwise() %>%
  mutate(experi = strsplit(Fastq, "_")[[1]][[1]]) %>%
  dplyr::select(-Fastq)
meta <- meta[!duplicated(meta), ]

indir <- "/home/jyeung/hpc/intestinal_scchic/raw_data/HVG-intestines.tagged_bams.all_merged/countTables"

# jmark <- "k4me1"
jdate <- "2019-12-20"

# get LH counts
jname <- "Scraped"
jgrp.main <- "HVG_Scraped_Bl6_.*."
jgrp.rz <- paste0(jgrp.main, "LHcounts.csv")
(infs.lh <- list.files(path = indir, pattern = jgrp.rz, all.files = TRUE, full.names = TRUE))

# assertthat::assert_that(file.exists(infs.lh))
sapply(infs.lh, function(inf) assertthat::assert_that(file.exists(inf)))

jmarks <- sapply(infs.lh, function(x) strsplit(basename(x), split = "_")[[1]][[4]])
names(infs.lh) <- jmarks
names(jmarks) <- jmarks

blfile <- "/home/jyeung/data/from_rstudioserver/scchic/databases/blacklist/mm10.blacklist.bed.gz"
assertthat::assert_that(file.exists(blfile))



# Constants ---------------------------------------------------------------


counts.cutoff.low <- 500  # H3K4me3
counts.cutoff.high <- 1000   # other marks 
TA.cutoff <- 0.5
counts.cutoff.lst <- vector(mode = "list", length = length(jmarks))
names(counts.cutoff.lst) <- jmarks

pcutoff.bin <- 0.95


outmain <- "/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/quality_controls_intestines"
dir.create(outmain)
outpdf <- file.path(outmain, paste0("quality_controls_", jname, "_all_marks.", Sys.Date(), ".pdf"))

# Plot data ---------------------------------------------------------------

dats.rz <- ReadLH.SummarizeTA(infs.lh, remove.nones = FALSE, na.to.zero = TRUE, bind.rows = FALSE)

# label mark
dats.rz <- lapply(jmarks, function(jmark){
  jtmp <- dats.rz[[jmark]]
  jtmp$mark <- jmark
  return(jtmp)
}) %>%
  bind_rows() %>%
  rowwise() %>% 
  mutate(experi = ClipLast(samp, jsep = "_"))

# label runs 
dats.rz <- left_join(dats.rz, meta) %>%
  rowwise() %>%
  mutate(experiRun = ifelse(startsWith(experi, "HVG-"), paste(Run, experi, sep = "-"), experi))

# plot good cels 
m.all <- ggplot(dats.rz, aes(x = log10(total.count), y = TA.frac)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

m.dens <- ggplot(dats.rz, aes(x = log10(total.count), fill = mark)) + geom_density() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark, ncol = 1)

# split by Lgr5? 
for (jmark in jmarks){
  print(jmark)
  x <- table(subset(dats.rz, mark == jmark)$experiRun)
  print(x)
}

dat.sum <- dats.rz %>%
  group_by(experiRun) %>%
  summarise(ncells = length(samp))

dat.sum2 <- dats.rz %>%
  group_by(experi) %>%
  summarise(ncells = length(samp)) %>%
  arrange(desc(ncells))

hist(dat.sum$ncells)

# add empty wells
empty.wells <- GetEmptyWells()

dats.rz$cellindx <- sapply(dats.rz$samp, function(x) paste("cell", strsplit(x, "_")[[1]][[2]], sep = ""))

dats.rz$is.empty <- dats.rz$cellindx %in% empty.wells

marks.with.low.counts <- c("k4me3", "k36me3")
for (jmark in jmarks){
  counts.cutoff.lst[[jmark]] <- ifelse(jmark %in% marks.with.low.counts, counts.cutoff.low, counts.cutoff.high)
}

dats.rz <- dats.rz %>%
  rowwise() %>%
  mutate(is.good = total.count > counts.cutoff.lst[[mark]] & TA.frac > TA.cutoff & !is.empty)



table(subset(dats.rz, mark == "k27me3")$experiRun)


# Define cells.keep -------------------------------------------------------

cells.keep <- subset(dats.rz, is.good)$samp

# Get good bins -----------------------------------------------------------





# Load mats ---------------------------------------------------------------

# get LH counts
jgrp.counts <- paste0(jgrp.main, "countTable.csv")
# jgrp.counts <- paste0("HVG_Scraped_Bl6_", ".*", "_", jdate, ".countTable.csv")
infs.counts <- list.files(path = indir, pattern = jgrp.counts, all.files = TRUE, full.names = TRUE)

jmarks.counts <- sapply(infs.counts, function(x) strsplit(basename(x), "_")[[1]][[4]])
names(jmarks.counts) <- jmarks.counts

assertthat::assert_that(identical(jmarks, jmarks.counts))

names(infs.counts) <- jmarks.counts

sapply(infs.counts, function(inf) assertthat::assert_that(file.exists(inf)))

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
mats <- lapply(infs.counts, function(inf){
  jmat <- ReadMatSlideWinFormat(inf)
  # filter bad chromosomes
  rows.i <- sapply(rownames(jmat), function(x) strsplit(x, ":")[[1]][[1]] %in% jchromos)
  cols.i <- which(colnames(jmat) %in% cells.keep)
  if (length(cols.i) > 0){
    return(jmat[rows.i, cols.i])
  } else{
    return(NULL)
  }
})

rnames.all.common <- Reduce(intersect, lapply(mats, function(mat) rownames(mat)))


# Get bin means -----------------------------------------------------------

bin.lst <- lapply(jmarks, function(jmark){
  Matrix::rowMeans(mats[[jmark]])
})


# Filter blacklist --------------------------------------------------------

rnames.gr <- makeGRangesFromDataFrame(data.frame(seqnames = sapply(rnames.all.common, GetChromo, add.chr=FALSE),
                                                 start = sapply(rnames.all.common, GetStart),
                                                 end = sapply(rnames.all.common, GetEnd)))
bl.gr <- LoadBlacklist(inf = blfile)

overlaps <- findOverlaps(bl.gr, rnames.gr)

indx <- seq(length(rnames.gr))
bl.hits.i <- unique(subjectHits(overlaps))
bl.hits.l <- !indx %in% bl.hits.i
rnames.gr.filt <- rnames.gr[bl.hits.l]
rnames.gr.badbins <- rnames.gr[!bl.hits.l]


# Use filtered blacklist --------------------------------------------------

rnames.all.common.blfilt <- names(rnames.gr.filt)



# Plot bins ---------------------------------------------------------------

bin.long.nofilt <- lapply(jmarks, function(jmark){
  x <- bin.lst[[jmark]]
  dat <- data.frame(coord = names(x), bin.mean = x, mark = jmark, stringsAsFactors = FALSE)
}) %>%
  bind_rows()
bin.wide.nofilt <- spread(bin.long.nofilt, key = mark, value = bin.mean)

ggpairs(log10(bin.wide.nofilt %>% dplyr::filter(coord %in% names(rnames.gr.badbins)) %>% dplyr::select(-coord)),
        lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
  ggtitle("Blacklist bin signal")

# filter out blacklist before doing correlations
bin.long <- lapply(jmarks, function(jmark){
  x <- bin.lst[[jmark]]
  dat <- data.frame(coord = names(x), bin.mean = x, mark = jmark, stringsAsFactors = FALSE)
  dat <- subset(dat, coord %in% rnames.all.common.blfilt)
}) %>%
  bind_rows() %>%
  group_by(mark) %>%
  mutate(bin.rank = rank(bin.mean) / length(bin.mean))

bin.sum <- bin.long %>%
  group_by(coord) %>%
  summarise(bin.rank.mean = mean(bin.rank),
            bin.rank.min = min(bin.rank),
            bin.rank.max = max(bin.rank))


bin.wide <- spread(bin.long %>% dplyr::select(-bin.rank), key = mark, value = bin.mean)

ggpairs(log10(bin.wide %>% dplyr::select(-coord)), lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
  ggtitle("Before removing correlated bins")

# remove bad bins?
bins.high <- subset(bin.sum, bin.rank.mean >= pcutoff.bin)$coord
bins.low <- subset(bin.sum, bin.rank.mean <= (1-pcutoff.bin))$coord

bins.correlated <- c(bins.high, bins.low)
print(length(bins.correlated))

m.pairs.clean <- ggpairs(log10(bin.wide %>% dplyr::filter(!coord %in% bins.correlated) %>% dplyr::select(-coord)),
                         lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
  ggtitle("After removing correlated bins")


# check any systematic biases on the offset
bsize <- 100000
bin.sum$bin <- as.numeric(sapply(bin.sum$coord, GetStart))
bin.sum$offset <- sapply(bin.sum$bin, function(x) x %% bsize)

ggplot(bin.sum, aes(x = bin.rank.mean, group = offset)) + geom_histogram() + facet_wrap(~offset, ncol = 1) +
  theme_bw() + theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



plot(table(sapply(bins.high, GetChromo)))
plot(table(sapply(bins.low, GetChromo)))

rnames.blfilt.corrfilt <- (bin.wide %>% dplyr::filter(!coord %in% bins.correlated))$coord


unique(sapply(rnames.blfilt.corrfilt, function(x) strsplit(x, ":")[[1]][[1]]))

print(length(rnames.blfilt.corrfilt))



# Write correlated bins to output -----------------------------------------

# bins.correlated.dat <- data.frame(bin = bins.correlated)
bins.correlated.dat <- data.frame(chromo = paste0("chr", sapply(bins.correlated, GetChromo)),
                                  start = sapply(bins.correlated, GetStart),
                                  end = sapply(bins.correlated, GetEnd))

outbed <- file.path(outmain, paste0("correlated_bins.", jname, ".bincutoff_", pcutoff.bin, ".", Sys.Date(), ".bed"))


fwrite(bins.correlated.dat, file = outbed, col.names = FALSE, sep = "\t")

for (jmark in jmarks){
  cells.keep.dat <- data.frame(cell = cells.keep, stringsAsFactors = FALSE)
  fwrite(cells.keep.dat, file = file.path(outmain, paste0("cellskeep_", jname, ".txt")), col.names = FALSE, sep = "\t")
}



# Write matrix to output  -------------------------------------------------

# write pdf
pdf(file = file.path(outmain, paste0("qc_plots_", jname, ".TAcutoff_", TA.cutoff, ".countscutoff_", counts.cutoff.low, "_", counts.cutoff.high, ".binfilt_cellfilt.", Sys.Date(), ".pdf")))
print(m.all)
print(m.dens)
for (jmark in jmarks){
  m <- ggplot(dats.rz %>% filter(mark == jmark), aes(x = log10(total.count), y = TA.frac, color = is.good, size = is.empty, shape = is.empty)) + geom_point() + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~interaction(experiRun)) + ggtitle(jmark)
  print(m)
}
ggpairs(log10(bin.wide %>% dplyr::filter(!coord %in% bins.correlated) %>% dplyr::select(-coord)),
        lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
  ggtitle("After removing correlated bins")

ggpairs(log10(bin.wide %>% dplyr::select(-coord)), lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
  ggtitle("Before removing correlated bins")
dev.off()


# Do final bin filter on matrix -------------------------------------------

# filter out
mats.binfilt.lst <- lapply(mats, function(mat){
  return(mat[rnames.blfilt.corrfilt, ])  # cells already filtered
})

lapply(mats.binfilt.lst, dim)


# Measure sparsity --------------------------------------------------------

lapply(mats.binfilt.lst, function(mat){
  1 - Matrix::nnzero(mat) / length(mat)
})


# Write merged all --------------------------------------------------------


# WT+Lgr5
for (jmark in jmarks){
  print(jmark)
  mat <- mats.binfilt.lst[[jmark]]
  print(dim(mat))
  print(1 - Matrix::nnzero(mat) / length(mat))
  outf.base <- file.path(outmain, paste0("mat.", jname, ".AllMerged.", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", counts.cutoff.low, "_", counts.cutoff.high, ".binfilt_cellfilt.", Sys.Date()))
  saveRDS(mat, file = paste0(outf.base, ".rds"))
  writeMM(mat, sparse = TRUE, file = paste0(outf.base, ".mm"))
  write.table(rownames(mat), file = paste0(outf.base, ".rownames"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
  write.table(colnames(mat), file = paste0(outf.base, ".colnames"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
}

# WT
jgrp <- "Lgr5"
jgrp.name <- "Unenriched"
for (jmark in jmarks){
  outf.base <- file.path(outmain, paste0("mat.", jname, ".", jgrp.name, ".", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", counts.cutoff.low, "_", counts.cutoff.high, ".binfilt_cellfilt.", Sys.Date()))
  mat.tmp <- GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp, jgrp.name, paste0(outf.base, ".rds"), invert = TRUE)
  writeMM(mat.tmp, sparse = TRUE, file = paste0(outf.base, ".mm"))
  write.table(rownames(mat.tmp), file = paste0(outf.base, ".rownames"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
  write.table(colnames(mat.tmp), file = paste0(outf.base, ".colnames"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
}


# Lgr5
jgrp <- "Lgr5"
jgrp.name <- "StemCells"
jmarks.ignore <- c("k36me3", "k4me3")
for (jmark in jmarks){
  if (jmark %in% jmarks.ignore){
    print("Skipping k36me3, no stem cells done")
    next
  }
  print(jmark)
  outf.base <- file.path(outmain, paste0("mat.", jname, ".", jgrp.name, ".", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", counts.cutoff.low, "_", counts.cutoff.high, ".binfilt_cellfilt.", Sys.Date()))
  if (file.exists(paste0(outf.base, ".rds"))){
    print(paste0(outf.base, ".rds", " exists, skipping"))
    next
  }
  mat.tmp <- GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp, jgrp.name, paste0(outf.base, ".rds"))
  writeMM(mat.tmp, sparse = TRUE, file = paste0(outf.base, ".mm"))
  write.table(rownames(mat.tmp), file = paste0(outf.base, ".rownames"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
  write.table(colnames(mat.tmp), file = paste0(outf.base, ".colnames"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
}


print(Sys.time() - jstart)