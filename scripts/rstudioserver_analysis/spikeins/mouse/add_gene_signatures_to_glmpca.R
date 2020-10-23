# Jake Yeung
# Date of Creation: 2020-09-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse/add_gene_signatures_to_glmpca.R
# 

rm(list=ls())

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Load gene sets ----------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
inf.annot <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData")
load(inf.annot, v=T)

# jsuffix <- "TSS_to_TES"
jsuffix <- "TSS.bsize_10000"
outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_mouse.round2/geneset_glmpca/celltyping_by_genesets.", jsuffix, ".pdf")

# pdf(outpdf, useDingbats = FALSE)

# Losd spikeins -----------------------------------------------------------



# Load glmpca -------------------------------------------------------------


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks


# Load data  --------------------------------------------------------------

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_mouse_spikein_VAN5046_varfilt"
assertthat::assert_that(dir.exists(indir))

jvarfilt <- 1
jmaxiter <- 1000
jtol <- as.character("1e-6")


dat.umap.lst <- lapply(jmarks, function(jmark){
  
  print(jmark)
  inf <- file.path(indir, paste0("count_mat_", jmark, "_l2r_filt.2020-09-13.minl2r_-1.varfilt_", jvarfilt, ".glmpcaout.penalty_1.maxiter_", jmaxiter, ".stochastic.avagrad.tol_", jtol, ".devfilt_5000.by_plate.RData"))
  print(basename(inf))
  load(inf, v=T)
  
  dat.umap.long <- DoUmapAndLouvain(glmpcaout$factors, jsettings = jsettings)
  
  dat.umap.long$dim1 <- glmpcaout$factors[, 1]
  dat.umap.long$dim2 <- glmpcaout$factors[, 2]
  dat.umap.long$dim3 <- glmpcaout$factors[, 3]
  dat.umap.long$dim4 <- glmpcaout$factors[, 4]
  
  dat.umap.long$mark <- jmark
  
  return(dat.umap.long)
})

cells.keep.lst <- lapply(dat.umap.lst, function(jdat) jdat$cell)


# Load spikeins -----------------------------------------------------------


inf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.chromo2spikeinfilt.VAN5046/spikein_info.txt"
dat.spikeins.mat <- fread(inf.spikeins)

dat.merge <- left_join(dat.umap.lst %>% bind_rows(), subset(dat.spikeins.mat, select = c(samp, spikeincounts, chromocounts)), by = c("cell" = "samp")) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"), 
         plate = as.numeric(strsplit(experi, "-")[[1]][[6]]),
         rowcoord = AddPlateCoordinates(cell)$rowcoord,
         colcoord = AddPlateCoordinates(cell)$colcoord,
         stype = AnnotateSortFromLayout(plate, rowcoord, colcoord))

dat.merge$stype <- gsub("^LSK", "aLSK", dat.merge$stype)

# Load count tables -------------------------------------------------------

indir.counts <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM/tagged_bams/merged_bams/TSS_count_tables.NewFilters")
assertthat::assert_that(dir.exists(indir.counts))
infs <- lapply(jmarks, function(jmark){
  fname <- paste0("PZ-ChIC-mouse-BM-", jmark, "-merged.sorted.tagged.countTable.", jsuffix, ".csv")
  # fname <- paste0("PZ-ChIC-mouse-BM-", jmark, "-merged.sorted.tagged.countTable.TSS.bsize_10000.csv")
  inf <- file.path(indir.counts, fname)
  assertthat::assert_that(file.exists(inf))
  return(inf)
})
  
  
mats.lst <- lapply(infs, function(inf) ReadMatTSSFormat(inf, as.sparse = TRUE, add.coord = FALSE, sort.rnames = FALSE))
mats.filt.lst <- lapply(jmarks, function(jmark) mats.lst[[jmark]][ , cells.keep.lst[[jmark]]])

cuts.gbodies <- lapply(mats.filt.lst, function(jmat) data.frame(cell = colnames(jmat), genecounts = colSums(jmat)))

genes.full.all <- Reduce(intersect, lapply(mats.lst, function(jmat) rownames(jmat)))
genes.all <- sapply(genes.full.all, function(g){
  gsplit <- strsplit(strsplit(g, "\\..")[[1]][[3]], "_")[[1]][[1]]
  return(gsplit)
})

ens.all <- sapply(genes.all, function(g) AssignHash(g, g2e.hash, null.fill = NA), USE.NAMES = FALSE)

# Get granulocytes --------------------------------------------------------

# jmark <- "H3K27me3"

# jmark <- "H3K4me1"

jgsets <- names(de.ens.sorted.stringent)
names(jgsets) <- jgsets


# Get histograms ----------------------------------------------------------

jmark <- "H3K27me3"

jcnames.keep.lst <- lapply(split(dat.merge, f = dat.merge$stype), function(x) x$cell)
jmat.pbulk <- SumAcrossClusters(mats.filt.lst[[jmark]], cnames.keep.lst = jcnames.keep.lst)

jmat.pbulk.long <- bind_rows(jmat.pbulk) %>%
  as.data.frame()
jmat.pbulk.long$gene <- names(jmat.pbulk$aLSK) 
jmat.pbulk.long <- jmat.pbulk.long %>%
  melt()
colnames(jmat.pbulk.long) <- c("gene", "stype", "cuts")

dat.pbulk <- dat.merge %>%
  filter(mark == jmark) %>%
  group_by(mark, stype) %>%
  summarise(spikeincounts = sum(spikeincounts),
            chromocounts = sum(chromocounts))

jmat.pbulk.long.annot <- left_join(jmat.pbulk.long, dat.pbulk)  

# Plot marker genkes ------------------------------------------------------




pdf(outpdf, useDingbats = FALSE)

ggplot(jmat.pbulk.long.annot, aes(x = log2(cuts / spikeincounts), fill = stype)) + 
  geom_histogram()  + 
  facet_wrap(~stype, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jmat.pbulk.long.annot, aes(x = log2(cuts / spikeincounts), fill = stype)) + 
  geom_density()  + 
  facet_wrap(~stype, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


for (jgset in jgsets){
  for (jmark in jmarks){
    gset <- as.character(de.ens.sorted.stringent[[jgset]])
    ens.keep <- ens.all[which(ens.all %in% gset)]
    rnames.keep <- names(ens.keep)
    
    exprs.dat <- data.frame(cell = colnames(mats.filt.lst[[jmark]]), cuts = colSums(mats.filt.lst[[jmark]][rnames.keep, ]), stringsAsFactors = FALSE) %>%
      left_join(., dat.merge, by = "cell") %>%
      left_join(., cuts.gbodies[[jmark]], by = "cell")
    
    m1 <- ggplot(exprs.dat, aes(x = umap1, y = umap2, color = log2(cuts / genecounts))) + 
      geom_point() + 
      scale_color_viridis_c() + 
      facet_wrap(~stype) + 
      ggtitle(paste(jmark, jgset))  + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    
    m2 <- ggplot(exprs.dat, aes(x = umap1, y = umap2, color = log2(cuts / spikeincounts))) + 
      geom_point() + 
      scale_color_viridis_c() + 
      facet_wrap(~stype) + 
      ggtitle(paste(jmark, jgset))  + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m1)
    print(m2)
  }
}
dev.off()

