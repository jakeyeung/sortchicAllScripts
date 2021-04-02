# Jake Yeung
# Date of Creation: 2021-02-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/4-annotate_high_bins_gc.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(hash)
library(topicmodels)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

make.plots <- FALSE

binsize <- 50000

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps"
outpdf <- file.path(outdir, paste0("High_bins_gc_and_distance_to_gene.", Sys.Date(), ".pdf"))

if (make.plots){
  pdf(file = outpdf, useDingbats = FALSE)
}

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"

# merge some celltypes
ctypes <- list("Eryths" = "Erythroid", 
               "Bcells" = "Lymphoid", 
               "NKs" = "Lymphoid", 
               "Granulocytes" = "Myeloid",
               "Basophils" = "Myeloid", 
               "pDCs" = "Lymphoid",
               "DCs" = "Myeloid", 
               "HSPCs" = "HSPCs",
               "Erythroid" = "Erythroid",
               "Lymphoid" = "Lymphoid",
               "Myeloid" = "Myeloid")


merge.ctypes.by.lineage <- FALSE

# Load metadata -----------------------------------------------------------

indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage")

dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  fread(file.path(indir.meta, fname)) %>%
    rowwise() %>%
    mutate(lineage = ctypes[[cluster]])
}) 

# cluster to col
cluster2col <- hash(dat.metas$H3K4me1$cluster, dat.metas$H3K4me1$clustercol)
cluster2col[["Erythroid"]] <- "#0072B2"
cluster2col[["Lymphoid"]] <- "#56B4E9"
cluster2col[["Myeloid"]] <- "#D55E00"



# Load high bins ----------------------------------------------------------

indir.bins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned"
dat.high.bins <- lapply(jmarks, function(jmark){
  fname <- paste0("High_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.jmidbug_fixed.", jmark, ".2021-02-19.txt")
  inf <- file.path(indir.bins, fname)
  fread(inf)
})

bins.lst <- lapply(dat.high.bins, function(jdat){
  jdat$CoordOriginal
})

# Load GCs ----------------------------------------------------------------


# load GCs 
load("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_gc_analysis/gcs_genomewide.RData", v=T)

gr.gc.dat.filt.lst <- lapply(jmarks, function(jmark){
  jbins <- bins.lst[[jmark]]
  jout <- subset(gr.gc.dat, bname %in% jbins) %>%
    mutate(mark = jmark)
  jout <- jout[!duplicated(jout$bname), ]
}) 
  
gr.gc.dat.filt.long <- gr.gc.dat.filt.lst %>%
  bind_rows() 

gr.gc.dat.filt.long$mark <- factor(gr.gc.dat.filt.long$mark, levels = jmarks)

ggplot(gr.gc.dat.filt.long, aes(x = mark, y = gc)) + 
  geom_boxplot(outlier.alpha = 0.1) + 
  # geom_point() + 
  theme_bw()  + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# distance to nearest gene
dat.high.bins.long <- dat.high.bins %>%
  bind_rows()
dat.high.bins.long$mark <- factor(dat.high.bins.long$mark, levels = jmarks)

ggplot(dat.high.bins.long, aes(x = mark, y = abs(distanceToTSS) + 1)) + 
  geom_boxplot(outlier.alpha = 0.1) +
  scale_y_log10() + 
  geom_hline(yintercept = 25000, linetype = "dotted") + 
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# do bootstrapped estimates
set.seed(0)
nbootstraps.vec <- seq(1000); names(nbootstraps.vec) <- nbootstraps.vec
nsamps <- 5000

# dat.gc.bins.clean <- lapply(gr.gc.dat.filt.lst, function(jdat) jdat %>% dplyr::select(bname, gc))

dat.high.bins.subsamp.lst <- lapply(jmarks, function(jmark){
  jdat <- gr.gc.dat.filt.lst[[jmark]]
  jdat.subsamp.lst <- lapply(nbootstraps.vec, function(i){
    if (i %% 100 == 0){
      print(i)
    }
    rows.subsamp <- sample(x = seq(nrow(jdat)), size = nsamps, replace = FALSE)
    jdat.subsamp <- jdat[rows.subsamp, ]
    jdat.subsamp.sum <- jdat.subsamp %>%
      group_by(mark) %>%
      summarise(gc.mean = mean(gc)) %>%
      mutate(i = i)
  })
})


prob.low <- 0.001
prob.high <- 0.999
dat.vars.long.subsamp.long.sum.lst <- lapply(jmarks, function(jmark){
  dat.sum.subsamp <- dat.high.bins.subsamp.lst[[jmark]] %>%
    bind_rows() %>%
    group_by(mark) %>%
    summarise(jmedian = median(gc.mean),
              jlower = quantile(gc.mean, probs = prob.low),
              jupper = quantile(gc.mean, probs = prob.high))
})

dat.vars.long.subsamp.long.sum.lst

if (make.plots){
  dev.off()
}

# Plot --------------------------------------------------------------------




# Check GCs  --------------------------------------------------------------

