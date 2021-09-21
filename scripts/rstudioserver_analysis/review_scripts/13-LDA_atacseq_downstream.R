# Jake Yeung
# Date of Creation: 2021-08-28
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/13-LDA_atacseq_downstream.R
# 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(hash)
library(igraph)
library(umap)

hubprefix <- "/home/jyeung/hub_oudenaarden"


# Function ----------------------------------------------------------------

FixNames <- function(jpeak){
  ending <- strsplit(jpeak, "-")[[1]][[2]]
  endinghalf <- substr(ending, 1, nchar(ending) / 2)
  jstart <- strsplit(jpeak, "-")[[1]][[1]]
  jpeak.wrangled <- paste(jstart, endinghalf, sep = "-")
  jpeak.atac <- paste(jpeak.wrangled, jpeak.wrangled, sep = ";")
  return(jpeak.atac)
}

# Load dat ----------------------------------------------------------------


inf.atac <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisRevisions_countmat_peaks_filt.scATACseq.K-30/lda_outputs.countmat_peaks_filt.scATACseq.K-30.binarize.FALSE/ldaOut.countmat_peaks_filt.scATACseq.K-30.Robj"))
load(inf.atac, v=T)

tm.result <- posterior(out.lda)
dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)


# Add meta  ---------------------------------------------------------------

inf.meta <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Cusanovich_2018/metadata/cell_metadata.bonemarrow_only.withheader.txt")
dat.meta <- fread(inf.meta)
dat.meta.filt <- subset(dat.meta, select = c(cell, cell_label, tissue.replicate))
dat.meta.filt$barcode <- paste(dat.meta.filt$tissue.replicate, dat.meta.filt$cell, sep = ".")



# Annotate ----------------------------------------------------------------

dat.umap.annot <- left_join(dat.umap, dat.meta.filt)

ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = cell_label)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Check sortChIC ----------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3"); names(jmarks) <- jmarks

infs.chic.lst <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisRevisions_countmat_CusanovichPeaks.sortChIC.", jmark, ".K-30/lda_outputs.countmat_CusanovichPeaks.sortChIC.", jmark, ".K-30.binarize.FALSE/ldaOut.countmat_CusanovichPeaks.sortChIC.", jmark, ".K-30.Robj"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

out.lda.lst <- lapply(infs.chic.lst, function(jinf){
  load(jinf, v=T)
  return(out.lda)
})

count.mat.sortchic.lst <- lapply(infs.chic.lst, function(jinf){
  load(jinf, v=T)
  rownames(count.mat) <- sapply(rownames(count.mat), FixNames)
  return(count.mat)
})


tm.result.lst <- lapply(out.lda.lst, function(jout){
  posterior(jout)
})


dat.umap.lst <- lapply(tm.result.lst, function(jtm){
  DoUmapAndLouvain(jtm$topics, jsettings)
})


m.lst <- lapply(dat.umap.lst, function(jdat){
  ggplot(jdat, aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})



# Get meta ----------------------------------------------------------------

infs.meta <- lapply(jmarks, function(jmark){
  inf.meta <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt"))
  assertthat::assert_that(file.exists(inf.meta))
  return(inf.meta)
})

dat.meta.chic.lst <- lapply(infs.meta, function(jinf){
  dat.meta.tmp <- fread(jinf)
  
  dat.meta.tmp <- dat.meta.tmp %>%
    dplyr::select(c(cell, cluster, stypecol, clustercol, batch))
})

dat.umap.chic.annot.lst <- lapply(jmarks, function(jmark){
  left_join(dat.umap.lst[[jmark]], dat.meta.chic.lst[[jmark]])
})

m.celltypes.lst <- lapply(dat.umap.chic.annot.lst, function(jdat){
  ggplot(jdat, aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})


# Find granu peak ---------------------------------------------------------

jmarktmp <- "H3K4me1"
tm.mark <- tm.result.lst[[jmarktmp]]
tm.mark <- AddTopicToTmResult(tm.mark)
dat.topics <- OrderTopicsByEntropy(tm.result = tm.mark)

jtop <- "topic20"

jtops <- dat.topics$topic
for (jtop in jtops){
  topics.dat <- data.frame(cell = rownames(tm.mark$topics), topicweight = tm.mark$topics[, jtop], stringsAsFactors = FALSE)
  dat.annot <- left_join(dat.umap.chic.annot.lst[[jmarktmp]], topics.dat)
  m <- ggplot(dat.annot, aes(x = umap1, y = umap2, color = topicweight)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jtop) + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}

top.neutr <- "topic3"
# get top peaks
peak.weights <- sort(tm.mark$term[top.neutr, ], decreasing = TRUE)

jpeak <- names(peak.weights)[[2]]
ending <- strsplit(jpeak, "-")[[1]][[2]]
endinghalf <- substr(ending, 1, nchar(ending) / 2)
jstart <- strsplit(jpeak, "-")[[1]][[1]]
jpeak.wrangled <- paste(jstart, endinghalf, sep = "-")

# plot top peak 
dat.impute.log <- log2(t(tm.mark$topics %*% tm.mark$terms))

impute.peak <- data.frame(cell = colnames(dat.impute.log), peakimpute = dat.impute.log[jpeak, ], stringsAsFactors = FALSE)

dat.annot.impute <- left_join(dat.umap.chic.annot.lst[[jmarktmp]], impute.peak)

m <- ggplot(dat.annot.impute, aes(x = umap1, y = umap2, color = peakimpute)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jpeak) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)



# Check peak in atac ------------------------------------------------------

dat.impute.log.atac <- log2(t(tm.result$topics %*% tm.result$terms))

jpeak.atac <- paste(jpeak.wrangled, jpeak.wrangled, sep = ";")
impute.peak.atac <- data.frame(cell = colnames(dat.impute.log.atac), peakimpute = dat.impute.log.atac[jpeak.atac,  ])
dat.umap.atac.peakimpute <- left_join(dat.umap, impute.peak.atac)

m <- ggplot(dat.umap.atac.peakimpute, aes(x = umap1, y = umap2, color = peakimpute)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jpeak.atac) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)

# 
# tm.result.atac <- AddTopicToTmResult(tm.result)
# dat.topics.atac <- OrderTopicsByEntropy(tm.result.atac)
# 
# jtops <- dat.topics.atac$topic
# for (jtop in jtops){
#   topics.dat <- data.frame(cell = rownames(tm.result.atac$topics), topicweight = tm.result.atac$topics[, jtop], stringsAsFactors = FALSE)
#   dat.annot <- left_join(dat.umap.annot, topics.dat)
#   m <- ggplot(dat.annot, aes(x = umap1, y = umap2, color = topicweight)) + 
#     geom_point() + 
#     theme_bw() + 
#     ggtitle(jtop) + 
#     scale_color_viridis_c() + 
#     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#   print(m)
# }
# 
# # 11, 3, 9, 27 4
# jtop.atac <- "topic27"
# peak.weights.atac <- sort(tm.result.atac$terms[jtop.atac, ], decreasing = TRUE)
# jpeak.atac <- names(peak.weights.atac)[[1]]
# impute.peak.atac <- data.frame(cell = colnames(dat.impute.log.atac), peakimpute = dat.impute.log.atac[jpeak.atac,  ])



# Prove it  ---------------------------------------------------------------

m.umap.atac <- ggplot(dat.umap.annot %>% 
                        mutate(cell_label = ifelse(cell_label == "Hematopoietic progenitors", "CusanovichHSPCs", "notHSPCs")), 
                      aes(x = umap1, y = umap2, color = cell_label)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m.umap.atac.impute <- ggplot(dat.umap.atac.peakimpute, aes(x = umap1, y = umap2, color = peakimpute)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jpeak.wrangled) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m.umap.atac.impute)

m.umap.chic <- ggplot(dat.annot.impute%>% 
                        mutate(cluster = ifelse(cluster == "Granulocytes", "ZellerNeutrophils", "NotNeutrophils")), 
                      aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(paste(jmarktmp, "at scATACseq peaks")) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m.umap.chic)

m <- ggplot(dat.annot.impute, aes(x = umap1, y = umap2, color = peakimpute)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jpeak.wrangled) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m)




# Write for AvO -----------------------------------------------------------

outdir <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Cusanovich_2018/compare_atac_with_sortchic")

fwrite(dat.umap.annot, file = file.path(outdir, "Cusanovich_scATACseq_metadata.txt"))
saveRDS(tm.result.atac, file = file.path(outdir, "Cusanovich_scATACseq_LDA_at_ATAC_peaks.rds"))
saveRDS(count.mat, file = file.path(outdir, "Cusanovich_scATACseq_count_table_at_ATAC_peaks.rds"))

fwrite(dat.umap.chic.annot.lst$H3K4me1, file = file.path(outdir, "Zeller_sortChIC_H3K4me1_metadata.txt"))

# fix colnames of tm.mark
tm.mark$terms[1:5, 1:5]
fixed.names <- sapply(colnames(tm.mark$terms), function(peakname){
  FixNames(peakname)
})
colnames(tm.mark$terms) <- fixed.names

rownames(count.mat.sortchic.lst)

saveRDS(tm.mark, file = file.path(outdir, "Zeller_sortChIC_H3K4me1_LDA_at_ATAC_peaks.rds"))
saveRDS(count.mat.sortchic.lst$H3K4me1, file = file.path(outdir, "Zeller_sortChIC_H3K4me1_count_table_at_ATAC_peaks.rds"))
