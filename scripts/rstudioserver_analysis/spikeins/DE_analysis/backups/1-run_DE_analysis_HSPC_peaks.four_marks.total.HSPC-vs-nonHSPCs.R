# Jake Yeung
# Date of Creation: 2021-01-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/DE_analysis/1-run_DE_analysis_HSPC_peaks.four_marks.spikeins.HSPC-vs-nonHSPCs.R
# description
rm(list=ls())

library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(hash)
library(igraph)
library(umap)

library(JFuncs)

library(scchicFuncs)

jstart <- Sys.time()

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123



# Load LDA (contains countmat)  ---------------------------------------------------------------

ncores <- 8
hubprefix <- "/home/jyeung/hub_oudenaarden"
# jtype <- "bins"
# jdist <- "TSS"
jtype <- "HSPCpeaks"

# outdir <- "/home/jyeung/data/from_rstudioserver/spikein_fits_BM_poisson"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again.peaks"
dir.create(outdir)

# jmark <- "H3K4me1"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmarks <- c("H3K9me3"); names(jmarks) <- jmarks

indir.metas <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned"

inf.lst <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(indir.metas, "cell_cluster_table_with_spikeins.H3K4me1.2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
})

# inf.k4me1 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.H3K4me1.txt")
# inf.k4me3 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.H3K4me3.txt")
# inf.k27me3 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.H3K27me3_techrepmerged/BM_rep2_rep3reseq_H3K27me3.2020-12-10.txt")
# inf.k9me3 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/BM_celltypes.bincutoff_0.binskeep_1000.byplate.szname_none.niter_500.reorder_rownames.dupfilt.2020-11-23.H3K9me3.txt")
# assertthat::assert_that(file.exists(inf.k4me1))
# assertthat::assert_that(file.exists(inf.k4me3))
# assertthat::assert_that(file.exists(inf.k27me3))
# assertthat::assert_that(file.exists(inf.k9me3))
# 
# inf.lst <- list(H3K4me1 = inf.k4me1,
#                 H3K4me3 = inf.k4me3,
#                 H3K27me3 = inf.k27me3,
#                 H3K9me3 = inf.k9me3)

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

indir.cuts <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/bams_split_by_cluster/bams_merged_by_cluster/countmat_in_HSPCs_peaks")
for (jmark in jmarks){
  outf <- file.path(outdir, paste0("poisson_fit_", jtype, ".",  jmark, ".", Sys.Date(), ".newannot2.witherrors.total.newmetas.HSPCs_vs_nonHSPCs.RData"))
  if (file.exists(outf)){
    print(paste("File exists, skipping", outf))
    next
  }
  
  # load raw cuts
  infs.cuts <- list.files(path = indir.cuts, pattern = paste0(jmark, ".*.csv"), all.files = TRUE, full.names = TRUE)
  
  dat.csv.lst <- lapply(infs.cuts, function(inf.tmp){
    print(inf.tmp)
    mat.tmp <- ReadMatTSSFormat(inf.tmp, as.sparse = TRUE, add.coord = TRUE, sort.rnames = FALSE)
    rownames(mat.tmp) <- paste("chr", rownames(mat.tmp), sep = "")
    return(mat.tmp)
  })
  
  rnames.all <- sort(unique(unlist(lapply(dat.csv.lst, function(jmat) rownames(jmat)))))
  
  count.mat <- cbind.fill.lst(dat.csv.lst, all.rnames = rnames.all)
  
  inf.annot <- inf.lst[[jmark]]
  dat.annot <- fread(inf.annot) %>%
    rowwise() %>%
    dplyr::rename(plate.nbr = plate)
  
  # if (jmark == "H3K27me3"){
  #   # use cluster.fewer
  #   dat.annot <- dat.annot %>%
  #     dplyr::rename(cluster.more = cluster,
  #                   cluster = cluster.fewer)
  # }
  
  dat.annot <- dat.annot %>%
    rowwise() %>%
    mutate(plate = ClipLast(x = cell,jsep = "_"))
  
  cells.keep <- colnames(count.mat)
  dat.annot.filt <- subset(dat.annot, cell %in% cells.keep)
  
  # dat.umap.merge <- left_join(dat.umap, subset(dat.annot, select = c(cell, cluster)))
  # ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cluster)) + 
  #   geom_point() +  
  #   scale_color_manual(values = cbPalette) + 
  #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  # Run fits gene by gene ---------------------------------------------------
  
  dat.annot.filt2 <- subset(dat.annot.filt, !is.na(cluster))
  
  cells.keep <- dat.annot.filt2$cell
  
  print(jmark)
  jmat.mark <- count.mat[, cells.keep]
  
  dat.annots.filt.mark <- dat.annot.filt2 %>%
    filter(cell %in% cells.keep) %>%
    # mutate(Cluster = ifelse(cluster == "HSPCs", "aHSPCs", cluster)) %>%
    mutate(Cluster = ifelse(cluster == "HSPCs", "aHSPCs", "notHSPCs")) %>%
    rowwise() %>%
    mutate(batch = IsRound1(cell, mark = jmark)) %>%
    mutate(Plate = ifelse(jrep == "rep1old", "Round1", plate))
  
  # dat.annots.filt.mark <- dat.annot.filt2 %>%
  #   filter(cell %in% cells.keep) %>%
  #   mutate(Cluster = ifelse(cluster == "HSPCs", "aHSPCs", cluster)) %>%
  #   rowwise() %>%
  #   mutate(batch = IsRound1(cell, mark = jmark)) %>%
  #   mutate(Plate = ifelse(batch == "Round2", plate, "Round1"))
  
  print(unique(dat.annots.filt.mark$Plate))
  
  # ncuts.for.fit.mark <- data.frame(cell = colnames(jmat.mark), ncuts.total = colSums(jmat.mark), stringsAsFactors = FALSE)
  ncuts.for.fit.mark <- data.frame(cell = colnames(jmat.mark), ncuts.total = colSums(jmat.mark), stringsAsFactors = FALSE)
  # ncuts.for.fit.mark <- data.frame(cell = dat.annot.filt2$cell, ncuts.total = dat.annot.filt2$spikein_cuts, stringsAsFactors = FALSE)
  cnames <- colnames(jmat.mark)
  
  jrow.names <- rownames(jmat.mark)
  names(jrow.names) <- jrow.names
  
  
  jfits.lst <- parallel::mclapply(jrow.names, function(jrow.name){
    jrow <- jmat.mark[jrow.name, ]
    jout <- FitGlmRowClustersPlate(jrow, cnames, dat.annots.filt.mark, ncuts.for.fit.mark, jbin = NULL, returnobj = FALSE, with.se = TRUE)
    return(jout)
  }, mc.cores = ncores)
  
  
  # Ssave outputs -----------------------------------------------------------
  # saveRDS(jfits.lst, outf)
  save(jfits.lst, dat.annots.filt.mark, ncuts.for.fit.mark, jmat.mark, file = outf)
  
  print(Sys.time() - jstart)
}

print("Done")
print(Sys.time() - jstart)


