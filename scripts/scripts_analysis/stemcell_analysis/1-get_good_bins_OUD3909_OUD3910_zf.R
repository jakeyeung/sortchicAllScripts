# Jake Yeung
# Date of Creation: 2019-10-22
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/1-get_good_cells_OUD3909_OUD3910.R

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(scchicFuncs)
library(GGally)


# Constants ---------------------------------------------------------------

pcutoff <- 0.95

prefix <- "_ZFbonemarrow"
outdir <- paste0("~/data/scchic/count_mat_binfilt_cellfilt_for_LDA", prefix)
dir.create(outdir)


# Filter ------------------------------------------------------------------

inf <- "/Users/yeung/data/scchic/quality_control_ZF/mats_filt_by.good_cells.sizecutoff_2.5.TAcutoff_0.5.rds"
mats.filt <- readRDS(inf)

# Find bad bins -----------------------------------------------------------

# keep only chromosomes that are worth keeping
bins <- rownames(mats.filt[[1]])

chromos <- unique(sapply(bins, function(bin) strsplit(bin, ":")[[1]][[1]], USE.NAMES = FALSE))
chromos.keep <- paste("chr", seq(25), sep = "")

bins.chromofilt.i <- which(!startsWith(bins, "chrK"))

mats.filt <- lapply(mats.filt, function(jmat){
  return(jmat[bins.chromofilt.i, ])
})


# get bin means across marks
jmarks <- unique(sapply(names(mats.filt), function(x) strsplit(x, "-")[[1]][[4]]))
names(jmarks) <- jmarks

# grep for each mats and then get bins for each mark
bin.means.long <- lapply(jmarks, function(jmark){
  mats.tmp <- mats.filt[grep(jmark, names(mats.filt))]
  x <- Matrix::rowMeans(do.call(cbind, mats.tmp))
  dat <- data.frame(coord = names(x), bin.mean = x, mark = jmark, stringsAsFactors = FALSE)
  return(dat)
}) %>%
  bind_rows()

bin.wide.nofilt <- spread(bin.means.long, key = mark, value = bin.mean)

ggpairs(log10(bin.wide.nofilt %>% dplyr::select(-coord)), lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic()

bin.sum <- bin.means.long %>%
  group_by(mark) %>%
  mutate(bin.rank = rank(bin.mean) / length(bin.mean)) %>%
  group_by(coord) %>%
  summarise(bin.rank.mean = mean(bin.rank),
            bin.rank.min = min(bin.rank),
            bin.rank.max = max(bin.rank))

# remove bad bins?
bins.high <- subset(bin.sum, bin.rank.mean >= pcutoff)$coord
bins.low <- subset(bin.sum, bin.rank.mean <= (1-pcutoff))$coord


bins.correlated <- c(bins.high, bins.low)
print(length(bins.correlated))

ggpairs(log10(bin.wide.nofilt %>% dplyr::filter(!coord %in% bins.correlated) %>% dplyr::select(-coord)),
        lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
  ggtitle("After removing correlated bins")

# check any systematic biases on the offset
bsize <- 100000
bin.sum$bin <- as.numeric(sapply(bin.sum$coord, GetStart))
bin.sum$offset <- sapply(bin.sum$bin, function(x) x %% bsize)

ggplot(bin.sum, aes(x = bin.rank.mean, group = offset)) + geom_histogram() + facet_wrap(~offset, ncol = 1) +
  theme_bw() + theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Write good cells and good bins ------------------------------------------

# write them in the right way

print(names(mats.filt))

gstrs.prefix <- c("-ZFWKM-", "-ZFWKMCD41plus-", "-")
gstrs.suffix <- "-"

for (jmark in jmarks){
  for (gstr.prefix in gstrs.prefix){
    count.dat <- list()
    # fname.prefix <- gsub(pattern = "-", replacement = "", gstr.prefix)
    if (gstr.prefix == ""){
      fname.prefix <- "Merged"
    } else {
      fname.prefix <- gstr.prefix
    }
    gstr <- paste0(gstr.prefix, jmark)
    plates.tmp <- grep(gstr, names(mats.filt), value = TRUE)
    print(paste(gstr, "has", length(plates.tmp), "plates"))
    mats.tmp <- mats.filt[plates.tmp]
    # merge together then write as .RData
    mats.tmp.merge <- do.call(cbind, mats.tmp)
    count.dat$counts <- mats.tmp.merge
    save(count.dat, file = file.path(outdir, paste0("ZF", fname.prefix, "", jmark, "_pcutoff_", pcutoff, "_binfilt_cellfilt.", Sys.Date(), ".RData")))
  }
}


