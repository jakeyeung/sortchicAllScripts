# Jake Yeung
# Date of Creation: 2020-11-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/14-check_hiddendomains_signal_by_distance.R
# Is there a distance-dependent signal from TSS for K4me1 and K27me3? 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

# load hidden domains count table -----------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmark <- "H3K27me3"
inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins/lda_outputs.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.binarize.FALSE/ldaOut.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.Robj"))
load(inf.lda, v=T)

rnames <- rownames(count.mat)
bins <- paste("chr", sapply(rnames, function(x) strsplit(x, ";")[[1]][[1]], USE.NAMES = FALSE), sep = "")

jinf.tss <- file.path(hubprefix, "jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed")
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

annots.out <- AnnotateCoordsFromList(coords.vec = bins, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)


# Load meta data ----------------------------------------------------------

# assign clusters
inf.meta <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins/cell_cluster_table_with_spikeins.", jmark, ".2020-11-18.dupfilt.txt"))
dat.meta <- fread(inf.meta)



# pseudobulk countmat -----------------------------------------------------

dat.meta.sub <- subset(dat.meta, cell %in% colnames(count.mat) & cluster != "Basophils" & !is.na(spikein_cuts))
cells.keep <- dat.meta.sub$cell
cluster.lst <- split(dat.meta.sub$cell, dat.meta.sub$cluster)

pbulk <- SumAcrossClusters(as.matrix(count.mat), cluster.lst) %>%
  bind_cols() %>%
  as.data.frame()
rownames(pbulk) <- bins


# Normalize by spikeins  --------------------------------------------------

# by cluster
dat.meta.sub.sum <- dat.meta.sub %>%
  group_by(cluster) %>%
  summarise(spikein_cuts = sum(spikein_cuts))

assertthat::assert_that(identical(colnames(pbulk), dat.meta.sub.sum$cluster))

pbulk.spikeinnorm <- sweep(pbulk, MARGIN = 2, STATS = dat.meta.sub.sum$spikein_cuts, FUN = "/")

plot(density(log2(pbulk.spikeinnorm$Bcells + 1)))
plot(density(log2(pbulk.spikeinnorm$Eryths + 1)))

# Normalize by length -----------------------------------------------------

GetLengthOfBin <- function(bin){
  # chr1:100067500-100097500 -> 30000
  jstart <- as.numeric(GetStart(bin))
  jend <- as.numeric(GetEnd(bin))
  return(jend - jstart)
}

jlengths <- sapply(bins, GetLengthOfBin)
jlengths.dat <- data.frame(bin = names(jlengths), jlengths, stringsAsFactors = FALSE)

pbulk.spikeinnorm.lengthnorm <- sweep(log2(pbulk.spikeinnorm * 10^6 + 1), MARGIN = 1, STATS = jlengths, FUN = "/")

plot(density(log2(pbulk.spikeinnorm.lengthnorm$Bcells)))
plot(density(log2(pbulk.spikeinnorm.lengthnorm$Eryths)))
plot(density(log2(pbulk.spikeinnorm.lengthnorm$HSPCs)))
plot(density(log2(pbulk.spikeinnorm.lengthnorm$NKs)))


# For Eryths is there a distance dependency -------------------------------

annots.dist <- subset(annots.out$regions.annotated, select = c(region_coord, distanceToTSS)) %>%
  left_join(., jlengths.dat, by = c("region_coord" = "bin"))

dat.pbulk.annot <- data.frame(bin = rownames(pbulk.spikeinnorm.lengthnorm), pbulk.spikeinnorm.lengthnorm, stringsAsFactors = FALSE) %>%
  left_join(., annots.dist, by = c("bin" = "region_coord")) 

ggplot(dat.pbulk.annot, aes(x = abs(distanceToTSS), y = Eryths)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pbulk.annot, aes(x = abs(distanceToTSS), y = Bcells)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pbulk.annot, aes(x = abs(distanceToTSS), y = HSPCs)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pbulk.annot, aes(x = abs(distanceToTSS), y = NKs)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pbulk.annot, aes(x = abs(distanceToTSS), y = pDCs)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pbulk.annot, aes(x = abs(distanceToTSS), y = DCs)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pbulk.annot, aes(x = Eryths, y = HSPCs)) + 
  geom_point(alpha = 0.1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Make pretety plot -------------------------------------------------------

mat <- fread(outtxt)

ggplot(mat, aes(x = H3K4me1_ChIP, y = H3K4me1_chic)) + 
  geom_bin2d(bins = 200) + 
  scale_fill_gradient(low = "grey90", high = "blue") + 
  theme_bw() + 
  geom_density_2d(color = "red", alpha = 0.25) + 
  theme(aspect.ratio=1) 

ggplot(mat, aes(x = H3K27me3_ChIP, y = H3K27me3_chic)) + 
  geom_bin2d(bins = 200) + 
  # scale_fill_continuous(type = "viridis") + 
  scale_fill_gradient(low = "grey90", high = "blue") + 
  theme_bw() + 
  geom_density_2d(color = "red", alpha = 0.25) + 
  theme(aspect.ratio=1) 

ggplot(mat, aes(x = H3K4me3_ChIP, y = H3K4me3_chic)) + 
  geom_bin2d(bins = 200) + 
  scale_fill_continuous(type = "viridis") + 
  theme_bw() + 
  geom_density_2d(color = "grey85", bins = 30) 

ggplot(mat, aes(x = H3K4me1_ChIP)) + 
  geom_density() + 
  theme_bw() 

ggplot(mat, aes(x = H3K4me3_chic)) + 
  geom_density() + 
  theme_bw() 

ggplot(mat, aes(x = H3K4me3_chic)) + 
  geom_density() + 
  theme_bw() 

ggplot(mat, aes(x = H3K4me1_chic)) + 
  geom_density() + 
  theme_bw() 
