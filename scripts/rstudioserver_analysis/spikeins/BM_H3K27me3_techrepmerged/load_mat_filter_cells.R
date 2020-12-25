# Jake Yeung
# Date of Creation: 2020-12-08
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_techrepmerged/load_mat_filter_cells.R
# Filter cells 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


GetMarkFromSamp2 <- function(samp){
  if (grepl("K4me1", samp)){
    return("H3K4me1")
  } else if (grepl("K4me3", samp)){
    return("H3K4me3")
  } else if (grepl("K27me3", samp)){
    return("H3K27me3")
  } else if (grepl("K9me3", samp)){
    return("H3K9me3")
  }
}

l2r.min <- -4
ta.min <- 0.5
cuts.min <- 1000
nfrac.max <- 0.25


jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")
hubprefix <- "/home/jyeung/hub_oudenaarden"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq"

jbase <- "BM_rep2rep3reseq_H3K27me3.cleaned"
jsuffix <- paste0("l2rmin_", l2r.min, "tamin_", ta.min, "cutsmin_", cuts.min, ".fraczeromax_", nfrac.max)
fname <- paste0(jbase, jsuffix)

outrds <- file.path(outdir, paste0(fname, ".rds"))
outtxt1 <- file.path(outdir, paste0(fname, ".metadata_good_cells.txt"))
outtxt2 <- file.path(outdir, paste0(fname, ".metadat_all_cells.txt"))

outpdf <- file.path(outdir, paste0(fname, ".pdf"))

pdf(outpdf, useDingbats = FALSE)

# Load mat rep 3 ----------------------------------------------------------------


indir1 <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/tagged_bams/counts_tables")

infs1 <- list.files(indir1, pattern = "*.binsize_50000.csv", full.names = TRUE)
infs1.rz <- list.files(indir1, pattern = "*.RZ.csv", full.names = TRUE)
infs1.chromo <- list.files(indir1, pattern = "*.NoChromo.csv", full.names = TRUE)

mats1 <- lapply(infs1, function(inf){
  print(inf)
  ReadMatSlideWinFormat(inf)
})

rzs1 <- lapply(infs1.rz, function(infrz){
  ReadLH.SummarizeTA(infrz)
}) %>%
  bind_rows()



chromos1 <- lapply(infs1.chromo, function(infchromo){
  GetChromoCounts(infchromo) %>%
    filter(chromo == jspikeinchromo)
}) %>%
  bind_rows()

# Calculate TA fraction ---------------------------------------------------


ggplot(rzs1, aes(x = log10(total.count), y = TA.frac)) + 
  geom_point(alpha = 0.25)  + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load mat rep 2  ---------------------------------------------------------

indir2 <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links/countTablesAndRZr1only_ByChromo.NewFilters.blfix")
assertthat::assert_that(dir.exists(indir2))

infs2 <- list.files(indir2, pattern = paste0(".*rep2.*H3K27me3.*50000.csv"), full.names = TRUE)
infs2.rz <- list.files(indir2, pattern = paste0(".*rep2.*H3K27me3.*RZ.csv"), full.names = TRUE)
infs2.chromo <- list.files(indir2, pattern = paste0(".*rep2.*H3K27me3.*NoChromo.csv"), full.names = TRUE)

mats2 <- lapply(infs2, function(inf){
  print(inf)
  ReadMatSlideWinFormat(inf)
})


rzs2 <- lapply(infs2.rz, function(infrz){
  ReadLH.SummarizeTA(infrz)
}) %>%
  bind_rows()


chromos2 <- lapply(infs2.chromo, function(infchromo){
  GetChromoCounts(infchromo) %>%
    filter(chromo == jspikeinchromo)
}) %>%
  bind_rows()



# Merge together get good cells -------------------------------------------

rzs.all <- rbind(rzs1 %>% mutate(batch = "rep3"), rzs2 %>% mutate(batch = "rep2"))

ggplot(rzs.all, aes(x = log10(total.count), y = TA.frac)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  facet_wrap(~batch) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

chromos.all <- rbind(chromos1 %>% mutate(batch = "rep3"),
                     chromos2 %>% mutate(batch = "rep2"))


# Find empty cells --------------------------------------------------------


rzs.all <- rzs.all %>%
  rowwise() %>%
  mutate(rowcoord = AddPlateCoordinates(samp)$rowcoord,
         colcoord = AddPlateCoordinates(samp)$colcoord,
         is.empty = rowcoord <= 8 & colcoord == 1,
         mark = GetMarkFromSamp2(samp)) %>%
  left_join(., chromos.all %>% dplyr::select(samp, counts, experi, totalcounts, spikeincounts, chromocounts), by = "samp")


ggplot(rzs.all, aes(x = log10(total.count), y = TA.frac, color = is.empty)) + 
  geom_point(alpha = 0.25) + 
  facet_wrap(~experi) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(rzs.all, aes(x = log10(total.count), y = TA.frac, color = is.empty)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(rzs.all, aes(x = log2(chromocounts / spikeincounts), y = TA.frac, color = is.empty)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  facet_wrap(~experi) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Merge count mats --------------------------------------------------------

all.rows <- lapply(c(mats1, mats2), function(jmat){
  return(rownames(jmat))
}) %>%
  unlist() %>%
  unique() 

mats.merge <- cbind.fill.lst(c(mats1, mats2), all.rnames = all.rows)

nnzeros <- apply(mats.merge, MARGIN = 2, function(jcol) Matrix::nnzero(jcol) / length(jcol))

nnzeros.dat <- data.frame(samp = names(nnzeros), frac.nzeros = nnzeros, stringsAsFactors = FALSE)

rzs.all2 <- left_join(rzs.all, nnzeros.dat)

ggplot(rzs.all2, aes(x = log2(chromocounts / spikeincounts), y = frac.nzeros)) + 
  geom_point(alpha = 0.1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(rzs.all2, aes(x = frac.nzeros)) + 
  geom_density() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(rzs.all2, aes(x = log2(chromocounts / spikeincounts), y = TA.frac, color = is.empty, size = frac.nzeros)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  facet_wrap(~experi) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Find over-digested cells  -----------------------------------------------


# thresholds

ggplot(rzs.all2, aes(x = log2(chromocounts / spikeincounts), fill = is.empty)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + 
  facet_wrap(~experi) + 
  geom_vline(xintercept = l2r.min) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(rzs.all2, aes(x = log10(chromocounts), fill = is.empty)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + 
  facet_wrap(~experi) + 
  geom_vline(xintercept = log10(cuts.min)) +
  # geom_hline(yintercept = ta.min) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(rzs.all2, aes(x = log10(chromocounts), y = TA.frac, color = is.empty, size = frac.nzeros)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  facet_wrap(~experi) + 
  geom_vline(xintercept = log10(cuts.min)) +
  geom_hline(yintercept = ta.min) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(rzs.all2, aes(x = log2(chromocounts / spikeincounts), y = frac.nzeros, color = is.empty)) + 
  geom_point(alpha = 0.1) + 
  theme_bw() + 
  geom_hline(yintercept = nfrac.max) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(rzs.all2, aes(x = frac.nzeros)) + 
  geom_density() + 
  theme_bw() + 
  geom_vline(xintercept = nfrac.max) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Make filters and write output -------------------------------------------

rzs.all2 <- rzs.all2 %>%
  rowwise() %>%
  mutate(l2r = log2(chromocounts / spikeincounts)) 

rzs.filt <- rzs.all2 %>%
  rowwise() %>%
  mutate(is.good = TA.count > ta.min & total.count > cuts.min & l2r > l2r.min & frac.nzeros < nfrac.max & !is.empty)

dim(rzs.all2)
dim(subset(rzs.filt, is.good))


# Write to output ---------------------------------------------------------

cells.keep <- subset(rzs.filt, is.good)$samp

assertthat::assert_that(length(cells.keep) == length(unique(cells.keep)))

mat.clean <- mats.merge[, cells.keep]

saveRDS(mat.clean, file = outrds)
fwrite(subset(rzs.filt, is.good) %>% mutate(cell = samp), outtxt1, sep = "\t")
fwrite(rzs.filt %>% mutate(cell = samp), outtxt2, sep = "\t")


dev.off()

