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

jmark <- "H3K4me1"

pdfout <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_K562_cellcycle/genomewide_summaries_with_spikeins.", Sys.Date(), ".pdf")
pdf(pdfout, width = 1020/72, height = 815/72, useDingbats = FALSE)

# Load chromos ------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters")
infs.chromo <- list.files(indir, pattern = "K562-EtOH-.*.csv", full.names = TRUE)

jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")

dat.chromos <- lapply(infs.chromo, function(inf){
  dat.filt.long <- GetChromoCounts(inf, spikeinchromo = jspikeinchromo, chromos.keep = NA)
}) %>%
  bind_rows()


# Load LH counts ----------------------------------------------------------

indir.lh <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/tagged_bams/RZcounts.NewFilters")
infs.lh <- list.files(indir.lh, pattern = "K562-EtOH-.*.csv", full.names = TRUE)

dat.lh <- lapply(infs.lh, function(inf){
  print(inf)
  dat.filt.long.lh <- ReadLH.SummarizeTA(inf)
}) %>%
  bind_rows() %>%
  mutate(experi = ClipLast(samp, jsep = "_"))
 

chromocounts <- subset(dat.chromos, chromo == "1", select = c(samp, chromocounts, spikeincounts))

dat.lh <- left_join(dat.lh, chromocounts)

dat.lh$mark <- sapply(dat.lh$experi, function(x) strsplit(x, "-")[[1]][[3]])

ggplot(dat.lh, aes(x = log2(chromocounts), y = TA.frac)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.lh, aes(x = log2(chromocounts), y = TA.frac, color = mark)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~experi)

ggplot(dat.lh, aes(x = log2(chromocounts), y = TA.frac, color = mark)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)




# Filter bad cells  -------------------------------------------------------

fracmin <- 0.5
chromocountmin <-list(H3K4me1 = 500, H3K4me3 = 500, H3K27me3 = 1000, H3K9me3 = 1000)
log2fcmin <- 2.5
dat.cutoffs <- data.frame(mark = names(chromocountmin), chromocountsmin = unlist(chromocountmin), stringsAsFactors = FALSE)


ggplot(dat.lh, aes(x = log2(chromocounts), y = TA.frac, color = mark)) + 
  geom_point(alpha = 0.25) + 
  geom_vline(mapping = aes(xintercept = log2(chromocountsmin)), data = dat.cutoffs) + 
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


ggplot(dat.lh, aes(x = log2(chromocounts), y = TA.frac, color = is.good)) + 
  geom_point(alpha = 0.25) + 
  geom_vline(mapping = aes(xintercept = log2(chromocountsmin)), data = dat.cutoffs) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(dat.lh, aes(x = log2(chromocounts/spikeincounts), y = TA.frac, color = is.good)) + 
  geom_point(alpha = 0.25) + 
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

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

for (jmark in jmarks){
  m <- ggplot(dat.lh %>% filter(endsWith(experi, suffix = "G1-G2") & mark == jmark), 
             aes(y = rowcoord, x = colcoord, size = log2(chromocounts / spikeincounts), color = log2(chromocounts / spikeincounts))) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=16 / 24, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_y_continuous(trans = "reverse", breaks = seq(from = 0, to = 16, by = 4))  + 
    scale_x_continuous(breaks = seq(from = 0, to = 24, by = 4)) + 
    facet_wrap(~experi) + 
    ggtitle(jmark)
  print(m)
}



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

ggplot(jsub, aes(x = cellcycle.str, y = chromocounts)) + 
  geom_boxplot() + 
  geom_point()  + 
  theme_bw(18) + 
  xlab("") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  # theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark, nrow = 1) + 
  scale_y_log10()

ggplot(jsub, aes(x = cellcycle.str, y = log2(chromocounts / spikeincounts))) + 
  geom_boxplot() + 
  geom_point()  + 
  theme_bw(18) + 
  xlab("") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  # theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark, nrow = 1)

ggplot(jsub, aes(x = log2(chromocounts / spikeincounts), fill = cellcycle.str)) + 
  geom_density(alpha = 0.25) + 
  theme_bw(18) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~mark, nrow = 1)


jsub.others <- dat.lh %>% filter(!endsWith(experi, suffix = "G1-G2")) %>%
  mutate(cellcycle.str = "0_G1_only") %>%
  filter(samp %in% good.cells) %>%
  mutate(repl = strsplit(experi, split = "-")[[1]][[length(strsplit(experi, split = "-")[[1]])]])

jsub.merge <- bind_rows(jsub, jsub.others)

ggplot(jsub.merge, aes(x = cellcycle.str, y = log2(chromocounts))) + 
  geom_boxplot() + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(jsub.merge, aes(x = cellcycle.str, y = log2(chromocounts / spikeincounts))) + 
  geom_boxplot() + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(jsub.merge, aes(x = experi, y = log2(chromocounts / spikeincounts))) + 
  geom_boxplot() + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(jsub.merge, aes(x = log2(chromocounts / spikeincounts), fill = cellcycle.str)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(jsub.merge, aes(x = cellcycle.str, y = log2(chromocounts))) + 
  geom_boxplot() + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(jsub.others, aes(x = repl, y = log2(chromocounts))) + 
  geom_boxplot() + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(jsub.others, aes(x = repl, y = log2(chromocounts / spikeincounts))) + 
  geom_boxplot() + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(jsub.others, aes(x = log2(chromocounts / spikeincounts), fill = experi)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)



# Plate to plate variability  ---------------------------------------------


cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

jsub.filt <- subset(jsub.merge, cellcycle.str %in% c("0_G1", "0_G1_only"))
# jsub.filt <- subset(jsub.merge)

ggplot(jsub.filt, aes(x = experi, y = log2(chromocounts), fill = mark)) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  scale_fill_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ggtitle("G1 cells across all plates")

ggplot(jsub.filt, aes(x = experi, y = log2(chromocounts / spikeincounts), fill = mark)) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  scale_fill_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ggtitle("G1 cells across all plates")




dev.off()





