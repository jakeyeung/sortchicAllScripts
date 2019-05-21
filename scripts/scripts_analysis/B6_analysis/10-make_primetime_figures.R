# Jake Yeung
# Date of Creation: 2019-05-15
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/10-make_primetime_figures.R
# Make primtime figures


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


# Load trajs --------------------------------------------------------------

load("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata", v=T)


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


jmark <- "H3K4me1"
keep.top.genes <- 150
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

colsvec <- list(H3K4me1 = "cyan1", H3K4me3 = "darkblue",  H3K27me3 = "darkorange1", H3K9me3 = "red1")

# jbin <- TRUE; kstr <- "25_30_40_50"
jbin.vec <- c(TRUE, TRUE, FALSE, FALSE)
kstr.vec <- c("25_30_40_50", "25_30_40_50", "30_40_50", "30_40_50")
names(jbin.vec) <- jmarks
names(kstr.vec) <- jmarks

# jmark <- "H3K4me1"

infs <- mapply(function(jbin, jmark, kstr) paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/ldaAnalysisBins_BinCellFilt2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.", jbin, ".no_filt/lda_out_meanfilt.B6_", jmark, "_pcutoff_0.CountThres0.K-", kstr, ".Robj"), 
                  jbin.vec, jmarks, kstr.vec)
lapply(infs, function(inf) assertthat::assert_that(file.exists(inf)))


kchoose <- "50"
out.objs <- lapply(jmarks, function(jmark) LoadLDABins(jmark, jbin = NA, top.thres = 0.995, inf = infs[[jmark]], convert.chr20.21.to.X.Y = FALSE, add.chr.prefix = TRUE, choose.k = kchoose))
# out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark, jbin = NA, top.thres = 0.995, inf = inf, convert.chr20.21.to.X.Y = FALSE, add.chr.prefix = TRUE, choose.k = kchoose), jmarks, infs, SIMPLIFY = FALSE)
# names(out.objs) <- jmarks
# print(paste("K:", out.objs$out.lda@k))

# save to output if not already
obj.outf <- "~/data/scchic/robjs/B6_objs/LDA_objects_all_marks.Rdata"
if (!file.exists(obj.outf)){
  save(out.objs, file = obj.outf)
}

# Plot the 4 umaps with proper colors  ------------------------------------

jsub <- dat.umap.long.trajs[[1]]
jcolor <- colsvec[[1]]

m <- PlotXYNoColor(jsub, xvar = "umap1", yvar = "umap2", jcol = jcolor, jsize = 1)
print(m)

mlst <- lapply(jmarks, function(jmark) PlotXYNoColor(dat.umap.long.trajs[[jmark]], xvar = "umap1", yvar = "umap2", jcol = colsvec[[jmark]], jsize = 1))




# Plot topics -------------------------------------------------------------

# eryth, bcells, nk cells, granu, progenitors

topics.vec <- c(7, 14, 40, 47, 48)  
hits <- c("Hbb-bs", "Il2ra", "Prf1", "S100a8", "Kit", "Gzmb", "Inpp4b", "Irf8", "S100a7a", "S100a6", "Tal1", "Sox6")
topname <- c("Erythroblasts", "B cells", "NK cells", "Granulocytes", "Progenitors")

 
# Plot hits for H3K4me1 ---------------------------------------------------

# impute dat
load("/Users/yeung/data/scchic/robjs/B6_objs/terms_filt_H3K4me1_bin_TRUE_k_50.RData", v=T)  # terms.filt
topics.mat <- as.data.frame(posterior(out.objs[[jmark]])$topics)
colnames(topics.mat) <- paste0("topic_", colnames(topics.mat))

imput.mat <- t(posterior(out.objs[[jmark]])$topics %*% posterior(out.objs[[jmark]])$terms)

terms.keep <- subset(terms.filt, gene %in% hits)$term

imput.sub <- tidyr::gather(data.frame(term = terms.keep, imput.mat[terms.keep, ]), key = "cell", value = "exprs", -term)
imput.sub$cell <- gsub("\\.", "-", imput.sub$cell)
imput.sub <- left_join(imput.sub, terms.filt)
imput.wide <- tidyr::spread(imput.sub %>% dplyr::select(gene, exprs, cell), key = "gene", value = "exprs")

topics.mat$cell <- rownames(topics.mat)

# for topics and gene
dat.merge <- left_join(dat.umap.long.trajs[[jmark]], topics.mat %>% dplyr::select(c("cell", paste0("topic_", topics.vec))))
dat.merge <- left_join(dat.merge, imput.wide)

pdf("/Users/yeung/data/scchic/pdfs/B6_figures/umaps_and_hits/umaps_and_hits.pdf", useDingbats = FALSE)
multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 4)

for (jtop in topics.vec){
  m1 <- PlotXYWithColor(dat.merge, xvar = "umap1", yvar = "umap2", cname = paste0("topic_", jtop), jsize = 2, jtitle = paste0("topic_", jtop), jcol = scales::muted("darkred"))
  print(m1)
}

for (hit in hits){
  print(hit)
  m1 <- PlotXYWithColor(dat.merge, xvar = "umap1", yvar = "umap2", cname = paste0("`", hit, "`"), jsize = 2, jtitle = hit, strip.ticks = TRUE)
  print(m1)
}


# Show H3K4me1 with louvain and trajectories ------------------------------

traj.jsize <- 2
ctypes <- c("eryth", "granu", "lymphoid", "mega", "nk", "Tcell")
cbPalette2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3")
m.trajs <- ggplot(dat.umap.long.trajs[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank()) +
  scale_color_manual(values=cbPalette2) + 
  xlab("") + ylab("")
for (ctype in ctypes){
  m.trajs <- m.trajs + geom_path(data = trajs[[jmark]][[ctype]], color = "gray25", size = traj.jsize, alpha = 0.5)
}
print(m.trajs)
dev.off()
