# Jake Yeung
# Date of Creation: 2019-09-05
# File: ~/projects/scchic/scripts/antibody_testing/new_adaptors.R
# Test new adaptors



library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(scchicFuncs)

library(ChIPseeker)
library(GenomicRanges)
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# library(org.Mm.eg.db)
# library(ChIPseeker)


library(here)

setwd(here::here())

# Functions ---------------------------------------------------------------

# 
# CollapseTaggedCountCnames <- function(dat, cnames.indx = 1:3){
#   # 2nd row often redundant, merge with first row
#   if (any(class(dat) == "data.table")){
#     colnames(dat)[cnames.indx] <- unlist(dat[1, ..cnames.indx], use.names = FALSE)
#   } else {
#     colnames(dat)[cnames.indx] <- unlist(dat[1, cnames.indx], use.names = FALSE)
#   }
#   dat <- dat[-1, ]
#   return(dat)
# }


# Load data ---------------------------------------------------------------

#' ## Check total counts

jtypes <- c("blockedadaptor", "forkedadaptor", "stdadaptor")
names(jtypes) <- jtypes

# jtype <- "blockedadaptor"
# inf <- paste0("/Users/yeung/data/scchic/from_cluster/new_adaptors/countTables/PZ-ChIC-K562-H3K27me3-", jtype, "24cells.countTable.csv")
# inf.rz <- paste0("/Users/yeung/data/scchic/from_cluster/new_adaptors/RZcounts/PZ-ChIC-K562-H3K27me3-", jtype, "24cells.RZ_counts.csv")

inf.counts <- lapply(jtypes, function(jtype){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/new_adaptors/countTables/PZ-ChIC-K562-H3K27me3-", jtype, "24cells.countTable.csv")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

inf.rzs <- lapply(jtypes, function(jtype){
  inf.rz <- paste0("/Users/yeung/data/scchic/from_cluster/new_adaptors/RZcounts/PZ-ChIC-K562-H3K27me3-", jtype, "24cells.RZ_counts.csv")
  assertthat::assert_that(file.exists(inf.rz))
  return(inf.rz)
})

dat.rz <- lapply(jtypes, function(jtype){
  inf.rz <- inf.rzs[[jtype]]
  dat.rz <- ReadDinuc(inf.rz)
  dat.rz$adaptor <- jtype
  return(dat.rz)
}) %>%
  bind_rows()

dat.sum <- dat.rz %>%
  group_by(samp, adaptor) %>%
  summarise(count.sum = sum(count, na.rm = TRUE))

#' Only one barcode should have reads, the counts from other barcodes should be bleed through or sequencing errors?
#' Outlier with >1 millino reads is probably the barcode that corresponds to the correct sequence.
ggplot(dat.sum, aes(x = count.sum)) + geom_histogram() + scale_x_log10() +  
  facet_wrap(~adaptor) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dat.TA <- dat.rz %>%
  filter(dinuc == "TA") %>%
  group_by(samp, dinuc, adaptor) %>%
  summarise(count = sum(count))

dat.merge <- left_join(dat.TA, dat.sum) %>%
  rowwise() %>%
  mutate(TA.frac = count / count.sum)



#' Only one barcode should have reads, the counts from other barcodes should be bleed through?
#' 
ggplot(dat.merge, aes(y = TA.frac, x = count.sum)) + 
  geom_point() +
  theme_bw() + 
  facet_wrap(~adaptor) + 
  scale_x_log10() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#' ## Get distribution of counts across bins

# need to remove blacklist 
dat.mat <- lapply(inf.counts, function(inf){
  dat <- as.data.frame(fread(inf))
  dat <- CollapseTaggedCountCnames(dat)
  dat$reference_name <- paste("chr", dat$reference_name, sep = "")
  rownames(dat) <- paste(dat$reference_name, paste(dat$start, dat$end, sep = "-"), sep = ":")
  dat[is.na(dat)] <- 0
  return(dat)
})

rnames.lst <- lapply(dat.mat, function(x) rownames(x))
rnames.all <- base::Reduce(f = intersect, x = rnames.lst)

rnames.gr <- makeGRangesFromDataFrame(data.frame(seqnames = sapply(rnames.all, GetChromo, add.chr=FALSE),
                                                 start = sapply(rnames.all, GetStart),
                                                 end = sapply(rnames.all, GetEnd))) 

bl.gr <- LoadBlacklist(inf = "/Users/yeung/data/databases/blacklist_regions/hg19-blacklist.v2.cut.bed")

overlaps <- findOverlaps(bl.gr, rnames.gr)

indx <- seq(length(rnames.gr))
bl.hits.i <- unique(subjectHits(overlaps))
bl.hits.l <- !indx %in% bl.hits.i

rnames.gr.filt <- rnames.gr[bl.hits.l]
rnames.gr.badbins <- rnames.gr[!bl.hits.l]

rnames.all.common.blfilt <- names(rnames.gr.filt)

dat.mat.filt <- lapply(dat.mat, function(x){
  # rows.keep <- which(rownames(x) %in% rnames.all.common.blfilt)
  dat.tmp <- x[rnames.all.common.blfilt, ]
  #  dat.tmp <- as.data.frame(x[rows.keep, ])
  # rownames(dat.tmp) <- rownames(x)[rows.keep]
  return(as.data.frame(dat.tmp[, -c(1,2,3)]))
}) %>%
  bind_cols()
rownames(dat.mat.filt) <- rnames.all.common.blfilt

# plot distribution for each 
cells.keep <- grepl("cells_0$", colnames(dat.mat.filt))
dat.mat.filt.cellfilt <- dat.mat.filt[, cells.keep]
dat.mat.filt.cellfilt$region <- rownames(dat.mat.filt.cellfilt)

dat.mat.filt.cellfilt.long <- dat.mat.filt.cellfilt %>%
  tidyr::gather(., key = "cell", value = "counts", -region)

dat.mat.filt.cellfilt.long$cname <- gsub("PZ-ChIC-K562-H3K27me3-", "", dat.mat.filt.cellfilt.long$cell)


#' Look at signal across bins genome wide. Before we would have a background signal and a foreground signal.
#' It looks like blocked adaptors do not have the background signal anymore (means higher background in blocked adaptors??)

jcutoff <- 30
m <- ggplot(dat.mat.filt.cellfilt.long, aes(x = counts)) + 
  geom_histogram(bins = 50) + 
  facet_wrap(~cname, ncol = 1) + 
  theme_bw(8) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + 
  geom_vline(xintercept = jcutoff, linetype = "dotted") + 
  ggtitle("H3K27me3 ChIC signal across bins (blacklist filtered)")
print(m)

#' ## Check if background bins are indeed lower signal in blocked adaptors

# bg.bins <- subset(dat.mat.filt.cellfilt.long, cname == "stdadaptor24cells_0" & counts < jcutoff)$region

# cname.ref <- "forkedadaptor24cells_0"
cname.ref <- "stdadaptor24cells_0"
bg.bins <- subset(dat.mat.filt.cellfilt.long, cname == cname.ref & counts < jcutoff)$region
assertthat::assert_that(length(bg.bins) > 0)

# plot pairs
dat.forpairs <- subset(dat.mat.filt.cellfilt, region %in% bg.bins, select = -region)

GGally::ggpairs(log10(dat.forpairs)) + ggtitle(paste("Reference:", cname.ref))

m.filt <- ggplot(dat.mat.filt.cellfilt.long %>% mutate(is.bg = region %in% bg.bins), aes(x = counts, fill = is.bg)) + 
  geom_histogram(bins = 50, position = "dodge") + 
  facet_wrap(~cname, ncol = 1) + 
  theme_bw(8) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  # scale_x_log10() + 
  geom_vline(xintercept = jcutoff, linetype = "dotted") + 
  ggtitle("Sackground bins (defined by stdadaptor)") + 
  xlim(c(0, 500))
print(m.filt)
