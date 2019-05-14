# Jake Yeung
# Date of Creation: 2019-05-10
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/4-create_clusters_write_objects_H3K9me3.R
# H3K9me3


rm(list=ls())

library(JFuncs)
library(dplyr)
library(ggplot2)
library(data.table)
library(topicmodels)
library(tidytext)
library(umap)
library(ggrepel)
library(tidyr)

library(hash)
library(igraph)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/PlotFunctions.R")

# Constants ---------------------------------------------------------------

outdir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"

# Load bulk ---------------------------------------------------------------

inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"
dat <- fread(inf.bulkdat, sep = "\t")

colnames(dat) <- gsub(" ", "_", colnames(dat))

dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))

# Load files --------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# jbin <- TRUE; kstr <- "25_30_40_50"
jbin <- FALSE; kstr <- "30_40_50"
keep.top.genes <- 150

jmark <- "H3K9me3"
inf <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/ldaAnalysisBins_BinCellFilt2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.", jbin, ".no_filt/lda_out_meanfilt.B6_", jmark, "_pcutoff_0.CountThres0.K-", kstr, ".Robj")
assertthat::assert_that(file.exists(inf))

# Process LDA -------------------------------------------------------------

kchoose <- "auto"
# kchoose <- 50
out.objs <- LoadLDABins(jmark, jbin = NA, top.thres = 0.995, inf = inf, convert.chr20.21.to.X.Y = FALSE, add.chr.prefix = TRUE, choose.k = kchoose)
# load(inf, v=T)
print(paste("K:", out.objs$out.lda@k))

# Plot UMAP ---------------------------------------------------------------

jmetric.louv='euclidean'
jmindist.louv=0.25
jseed.louv=123

# nn.louv.new <- c(150, 100, 33, 31)

jmindist <- jmindist.louv
nn <- 47

custom.settings <- GetUmapSettings(nn, jmetric.louv, jmindist, jseed.louv)
dat.umap <- umap(out.objs$tm.result$topics, config = custom.settings)

dat.umap.long <- data.frame(cell = rownames(dat.umap$layout), umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], mark = jmark)

# PlotXYWithColor(dat.umap.long, xvar = "umap1", yvar = "umap2")

ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(nn)


# Do louvain clustering ---------------------------------------------------

# x <- 30
x <- 60
# x <- 48
jseed <- jseed.louv
# jsettings <- GetUmapSettings(x, jmetric.louv, jmindist.louv, jseed)
louv.settings <- custom.settings
louv.settings$n_neighbors <- x
louv.hash <- DoLouvain(out.objs$tm.result$topics, louv.settings)

# add to umap
dat.umap.long$louvain <- as.character(sapply(as.character(dat.umap.long$cell), function(x) louv.hash[[x]]))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "#D3D3D3")
m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + ggtitle(paste(x, jseed))
print(m)


# Save objects and write files  -------------------------------------------

pdf(file.path(outdir, paste0("umap_cluster_IDs.", jmark, ".pdf")))
print(m)
dev.off()

dat.sub <- subset(dat.umap.long, select = c(cell, louvain))
dat.sub$fname <- sapply(dat.sub$cell, function(x) paste0(x, ".sorted.bam"))
for (jclst in unique(base::sort(dat.sub$louvain))){
  dat.subsub <- subset(dat.sub, louvain == jclst)
  outf.sub <- file.path(outdir, paste0(paste("bamlist", jclst, jmark, sep = "-"), ".txt"))
  data.table::fwrite(dat.subsub %>% dplyr::select(fname), file = outf.sub, sep = "\t", col.names = FALSE)
  # if (file.exists(outf.sub)){
  #   next
  # } else {
  #   data.table::fwrite(dat.subsub %>% dplyr::select(fname), file = outf.sub, sep = "\t", col.names = FALSE)
  # }
}

save(out.objs, dat.umap.long, custom.settings, louv.settings,file = file.path(outdir, paste0("dat_umap_long_with_louvain.", jmark, ".RData")))



