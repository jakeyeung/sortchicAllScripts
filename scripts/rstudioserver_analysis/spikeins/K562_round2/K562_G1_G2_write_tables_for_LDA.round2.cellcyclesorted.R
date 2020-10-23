# Jake Yeung
# Date of Creation: 2020-08-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_G1_G2.R
# Find differences within G1 and G2, make data ready for LDA 
# and Poisson GLMPCA? 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

# wwrite to output
jsuffix <- "cellcyclefilt"
hubprefix <- "/home/jyeung/hub_oudenaarden"
outdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2")
dir.create(outdir)

outpdf <- file.path(outdir, paste0("K562_", jsuffix, "_filtering_and_qc.pdf"))

indir.chromo <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/tagged_bams/cellcycle_sorted/merged_bams/countTablesAndRZr1only_ByChromo.NewFilters")
indir.lh <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/tagged_bams/cellcycle_sorted/merged_bams/countTablesAndRZr1only_TAfrac.NewFilters")
indir.counts <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/tagged_bams/cellcycle_sorted/merged_bams/countTablesAndRZr1only_CountTables.NewFilters")

pdf(outpdf, useDingbats = FALSE)

# Load chromos ------------------------------------------------------------

# indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN6969/K562/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters")
infs.chromo <- list.files(indir.chromo, pattern = "K562-EtOH-.*.csv", full.names = TRUE)
# infs.chromo <- list.files(indir.chromo, pattern = "*.csv", full.names = TRUE)
assertthat::assert_that(length(infs.chromo) > 0)


jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")

dat.chromos <- lapply(infs.chromo, function(inf){
  dat.filt.long <- GetChromoCounts(inf, spikeinchromo = jspikeinchromo, chromos.keep = NA)
}) %>%
  bind_rows()


# Load LH counts ----------------------------------------------------------

# indir.lh <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN6969/K562/tagged_bams/RZcounts.NewFilters")
infs.lh <- list.files(indir.lh, pattern = "K562-EtOH-.*.csv", full.names = TRUE)

dat.lh <- lapply(infs.lh, function(inf){
  dat.filt.long.lh <- ReadLH.SummarizeTA(inf)
}) %>%
  bind_rows() %>%
  mutate(experi = ClipLast(samp, jsep = "_"))
 

chromocounts <- subset(dat.chromos, chromo == "1", select = c(samp, chromocounts, spikeincounts))

dat.lh <- left_join(dat.lh, chromocounts)

dat.lh$mark <- sapply(dat.lh$experi, function(x) strsplit(x, "-")[[1]][[3]])

ggplot(dat.lh, aes(x = log10(chromocounts), y = TA.frac)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.lh, aes(x = log10(chromocounts), y = TA.frac, color = mark)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~experi)

ggplot(dat.lh, aes(x = log10(chromocounts), y = TA.frac, color = mark)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)




# Filter bad cells  -------------------------------------------------------

fracmin <- 0.5
chromocountmin <-list(H3K4me1 = 500, H3K4me3 = 500, H3K27me3 = 1000, H3K9me3 = 1000)
log2fcmin <- 2.5
dat.cutoffs <- data.frame(mark = names(chromocountmin), chromocountsmin = unlist(chromocountmin), stringsAsFactors = FALSE)


ggplot(dat.lh, aes(x = log10(chromocounts), y = TA.frac, color = mark)) + 
  geom_point(alpha = 0.25) + 
  geom_vline(mapping = aes(xintercept = log10(chromocountsmin)), data = dat.cutoffs) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(dat.lh, aes(x = log2(chromocounts/spikeincounts), y = TA.frac, color = mark)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

dat.lh <- dat.lh %>%
  mutate(is.good = chromocounts > chromocountmin[[mark]] & TA.frac > fracmin & log2(chromocounts / spikeincounts) > log2fcmin)

dat.lh <- dat.lh %>%
  rowwise() %>%
  mutate(cellid = paste("cell", strsplit(samp, "_")[[1]][[2]], sep = ""),
         indx = strsplit(samp, "_")[[1]][[2]], 
         rowcoord = GetPlateCoord(cell = cellid, platecols = 24, is.zero.base = FALSE)[[1]],
         colcoord = GetPlateCoord(cell = cellid, platecols = 24, is.zero.base = FALSE)[[2]],
         is.empty = rowcoord <= 8 & colcoord == 1)

good.cells <- subset(dat.lh, is.good & !is.empty)$samp


ggplot(dat.lh, aes(x = log10(chromocounts), y = TA.frac, color = is.good)) + 
  geom_point(alpha = 0.25) + 
  geom_vline(mapping = aes(xintercept = log10(chromocountsmin)), data = dat.cutoffs) + 
  geom_hline(yintercept = fracmin) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(dat.lh, aes(x = log2(chromocounts/spikeincounts), y = TA.frac, color = is.good)) + 
  geom_point(alpha = 0.25) + 
  geom_vline(xintercept = log2fcmin) + 
  geom_hline(yintercept = fracmin) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

# Plot counts on the lpate  -----------------------------------------------


# ggplot(dat.lh.coord, aes(y = rowcoord, x = colcoord, size = log(spikeincounts), color = log(spikeincounts))) + 
ggplot(dat.lh, aes(y = rowcoord, x = colcoord, size = log(chromocounts), color = log(chromocounts))) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=16 / 24, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_continuous(trans = "reverse", breaks = seq(from = 0, to = 16, by = 4))  + 
  scale_x_continuous(breaks = seq(from = 0, to = 24, by = 4)) + 
  facet_wrap(~experi)


# look at G1, G2, S 

ggplot(dat.lh %>% filter(endsWith(experi, suffix = "G1-G2") & mark == "H3K9me3"), 
       aes(y = rowcoord, x = colcoord, size = log(chromocounts), color = log(chromocounts))) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=16 / 24, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_continuous(trans = "reverse", breaks = seq(from = 0, to = 16, by = 4))  + 
  scale_x_continuous(breaks = seq(from = 0, to = 24, by = 4)) + 
  facet_wrap(~experi)


cellcycle.lst = list("0" = "0_G1", "1" = "1_S", "2" = "2_G2/M")
cellcycle.hash <- hash::hash(cellcycle.lst)
jsub <- dat.lh %>% filter(endsWith(experi, suffix = "G1-G2")) %>%
  mutate(cellcycle = as.character(floor( ( colcoord - 1 ) / 8))) %>%
  rowwise() %>%
  mutate(cellcycle.str = cellcycle.hash[[cellcycle]]) %>%
  filter(samp %in% good.cells)

ggplot(jsub, aes(y = rowcoord, x = colcoord, color = cellcycle.str)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=16 / 24, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_continuous(trans = "reverse", breaks = seq(from = 0, to = 16, by = 4))  + 
  scale_x_continuous(breaks = seq(from = 0, to = 24, by = 4))

ggplot(jsub, aes(x = cellcycle.str, y = log10(chromocounts))) + 
  geom_boxplot() + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(jsub, aes(x = cellcycle.str, y = log2(chromocounts / spikeincounts))) + 
  geom_boxplot() + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(jsub, aes(x = log2(chromocounts / spikeincounts), fill = cellcycle.str)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)


# jsub.others <- dat.lh %>% filter(!endsWith(experi, suffix = "G1-G2")) %>%
#   mutate(cellcycle.str = "0_G1_only") %>%
#   filter(samp %in% good.cells) %>%
#   mutate(repl = strsplit(experi, split = "-")[[1]][[length(strsplit(experi, split = "-")[[1]])]])
# jsub.merge <- bind_rows(jsub, jsub.others)

ggplot(jsub, aes(x = cellcycle.str, y = log10(chromocounts))) + 
  geom_boxplot() + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(jsub, aes(x = cellcycle.str, y = log2(chromocounts / spikeincounts))) + 
  geom_boxplot() + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(jsub, aes(x = experi, y = log2(chromocounts / spikeincounts))) + 
  geom_boxplot() + 
  geom_point()  + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~mark)

ggplot(jsub, aes(x = log2(chromocounts / spikeincounts), fill = cellcycle.str)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(jsub, aes(x = cellcycle.str, y = log10(chromocounts))) + 
  geom_boxplot() + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)



# Load counts -------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

infs.counts <- list.files(indir.counts, pattern = paste0("K562-EtOH-.*.countTable.csv"), full.names = TRUE)
assertthat::assert_that(length(infs.counts) > 0)

dat.counts <- lapply(infs.counts, function(inf){
  print(inf)
  dat.filt.long.counts <- ReadMatSlideWinFormat(inf = inf, as.sparse = TRUE, sort.rnames = TRUE, add.chromo = TRUE)
})

# dat.counts.merge <- lapply(dat.counts, function(dat){
#   dat <- as.data.frame(as.matrix(dat))
#   dat$region <- rownames(dat)
#   return(dat)
# }) %>%
#   Reduce(f = full_join, x = .)

rows.all <- lapply(dat.counts, function(jdat){
  rownames(jdat)
})

rows.keep <- Reduce(intersect, rows.all)

chromos.keep <- paste("chr", seq(22), sep = "")
chromos.keep.grpstr <- paste(chromos.keep, collapse = "|")

rows.keep2 <- rows.keep[sapply(rows.keep, function(x) strsplit(x, split = ":")[[1]][[1]] %in% chromos.keep)]

dat.counts.merge <- lapply(dat.counts, function(jdat){
  good.cells.i <- colnames(jdat) %in% good.cells
  jdat[rows.keep2, good.cells.i]
})

dat.counts.merge <- do.call(cbind, dat.counts.merge)


# Write full table --------------------------------------------------------

dat.counts.merge.lst <- lapply(jmarks, function(jmark){
  cells.keep <- grepl(pattern = jmark, x = colnames(dat.counts.merge))
  dat.counts.merge[, cells.keep]
})


for (jmark in jmarks){
  print(jmark)
  outf <- file.path(outdir, paste0("K562_count_tables_50000.", jmark, ".", jsuffix, ".rds"))
  dat.tmp <- dat.counts.merge.lst[[jmark]]
  print(dim(dat.tmp))
  saveRDS(dat.tmp, file = outf)
}
# 
# # write output: G1,G2,S
# 
# jgrep <- "G1-G2$"
# for (jmark in jmarks){
#   print(jmark)
#   outf <- file.path(outdir, paste0("K562_count_tables_50000.", jmark, ".G1_G2_S.rds"))
#   dat.tmp <- dat.counts.merge.lst[[jmark]]
#   
#   cnames <- grepl(pattern = jgrep, sapply(colnames(dat.tmp), function(x) ClipLast(x, jsep = "_")))
#   dat.tmp <- dat.tmp[, cnames]
#   print(dim(dat.tmp))
#   
#   saveRDS(dat.tmp, file = outf)
# }
# 
# # write G1 only
# jgrep <- "G1-G2$"
# for (jmark in jmarks){
#   print(jmark)
#   outf <- file.path(outdir, paste0("K562_count_tables_50000.", jmark, ".G1only.rds"))
#   dat.tmp <- dat.counts.merge.lst[[jmark]]
#   
#   cnames <- !grepl(pattern = jgrep, sapply(colnames(dat.tmp), function(x) ClipLast(x, jsep = "_")))
#   dat.tmp <- dat.tmp[, cnames]
#   print(dim(dat.tmp))
#   
#   saveRDS(dat.tmp, file = outf)
# }


dev.off()


