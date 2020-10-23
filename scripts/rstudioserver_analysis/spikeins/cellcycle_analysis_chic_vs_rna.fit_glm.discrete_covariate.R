# Jake Yeung
# Date of Creation: 2020-08-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/cellcycle_analysis_chic_vs_rna.R
# Find cell cycle regulated genes? 


rm(list=ls())

jstart <- Sys.time()


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


library(ChIPseeker)
library(GenomicRanges)



# Load ChIC data  ---------------------------------------------------------




# Load data  --------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
# jmark <- "H3K4me1"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# outdir <- "/home/jyeung/data/from_rstudioserver/spikein_fits_cellcycle"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/spikein_fits_cellcycle.discrete_covariate"
dir.create(outdir)

# jsuffix <- "topn_5000.glmpcaout.penalty_5.RData"
jsuffix <- "topn_5000.glmpcaout.penalty_5.by_plate.RData"
jprefix <- ".G1_G2_S."

# Assign genes to bins  ---------------------------------------------------

# hubprefix 
jmark.test <- "H3K27me3"
jinf.tss <- file.path(hubprefix, "jyeung/data/databases/refseq/GRCh38/GCF_000001405.39_GRCh38.p13_genomic.parsed.TSS.txt.merged.chromorenamed.4columns.bed")
jinf.test <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/K562_count_tables_50000.", jmark.test, jprefix, "rds"))
mat.test <- readRDS(jinf.test)


jchromos.keep <- paste("chr", c(seq(22), "X", "Y"), sep = "")
annot.out <- AnnotateCoordsFromList.GeneWise(coords.vec = rownames(mat.test), inf.tss = jinf.tss, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene, annodb = "org.Hs.eg.db", chromos.keep = jchromos.keep)
annot.out$out2.df$geneclean <- sapply(annot.out$out2.df$gene, function(x) strsplit(x, "\\.\\.")[[1]][[2]])

# annot.out <- AnnotateBins2(terms.mat = , top.thres = 1, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene, annodb = "org.Hs.eg.db", inf.tss = jinf.tss)
  

mclapply(jmarks, function(jmark){
  print(jmark)
  
  print(Sys.time() - jstart)
  
  outpdf <- file.path(outdir, paste0("fit_cellcycle_pseudotime.", jmark, ".", Sys.Date(), ".pdf"))
  outrdata <- file.path(outdir, paste0("fit_cellcycle_pseudotime.", jmark, ".", Sys.Date(), ".RData"))
  
  # if (file.exists(outpdf)){
  #   next
  # }
  # if (file.exists(outrdata)){
  #   next
  # }
  
  
  
  # inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_K562_spikein/lda_outputs.K562_count_tables_50000.", jmark, jprefix, "K-30.binarize.FALSE/ldaOut.K562_count_tables_50000.", jmark, jprefix, "K-30.Robj"))
  # assertthat::assert_that(file.exists(inf.lda))
  
  # inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein/K562_count_tables_50000.", jmark, ".G1_G2_S.glmpcaout", jsuffix))
  inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein/K562_count_tables_50000.", jmark, jprefix, jsuffix))
  assertthat::assert_that(file.exists(inf.glmpca))
  
  
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/K562_count_tables_50000.", jmark, jprefix, "rds"))
  inf.spike <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/spikein_counts/spikein_counts_all.RData")
  
  assertthat::assert_that(file.exists(inf))
  assertthat::assert_that(file.exists(inf.spike))
  
  mat <- readRDS(inf)
  
  load(inf.spike, v=T)
  
  load(inf.glmpca, v=T)
  
  dat.spikeins.mat <- dat.spikeins.mat %>%
    rowwise() %>%
    mutate(mark = strsplit(samp, "-")[[1]][[3]]) %>%
    filter(mark == jmark)
  
  
  dat.spikeins.mat <- scchicFuncs::AddCellCycleLabel.bydat(dat.spikeins.mat)
  
  totalcounts <- data.frame(cell = colnames(mat), totalcounts = colSums(mat), stringsAsFactors = FALSE)
  
  dat.spikeins.mat <- left_join(dat.spikeins.mat, totalcounts, by = c("samp" = "cell")) %>%
    filter(!is.na(totalcounts))
  
  
  
  
  # Plot PCA ?  -------------------------------------------------------------
  
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- 0.1
  jsettings$random_state <- 123
  
  dat.umap <- DoUmapAndLouvain(glmpcaout$factors, jsettings)
  
  dat.umap$dim1 <- glmpcaout$factors$dim1
  dat.umap$dim2 <- glmpcaout$factors$dim2 
  
  
  dat.umap <- left_join(dat.umap, dat.spikeins.mat, by = c("cell" = "samp"))
  
  # add pseudotme
  dat.umap$pseudotime <- dat.umap$dim1
  dat.umap$pseudotime <- dat.umap$pseudotime - min(dat.umap$pseudotime)
  dat.umap$pseudotime <- dat.umap$pseudotime / max(dat.umap$pseudotime)
  
  pdf(outpdf, useDingbats = FALSE)
  
  # decide whether to flip pseudotime
  m <- ggplot(dat.umap, aes(x = pseudotime, y = log2(totalcounts / spikeincounts))) + 
    geom_point() + 
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  jcheck <- lm(log2(totalcounts / spikeincounts) ~ pseudotime, dat.umap)
  if (coef(jcheck)[["pseudotime"]] < 0){
    dat.umap$pseudotime <- dat.umap$pseudotime * -1
    dat.umap$pseudotime <- dat.umap$pseudotime + 1
  }
  
  m <- ggplot(dat.umap, aes(x = pseudotime, y = log2(totalcounts / spikeincounts), color = cellcycle.str)) + 
    geom_point() + 
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  m <- ggplot(dat.umap, aes(x = dim1, y = dim2, color = log2(totalcounts / spikeincounts))) + 
    geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c()
  print(m)
  
  
  # some interesting regions  ------------------------------------------------
  
  loadings <- glmpcaout$loadings$dim1
  names(loadings) <- rownames(glmpcaout$loadings)
  
  loadings <- sort(loadings, decreasing = FALSE)
  
  jbin <- names(loadings)[[100]]
  
  dat.counts <- data.frame(count = mat[jbin, ], cell = colnames(mat), stringsAsFactors = FALSE) %>%
    left_join(., dat.umap)
  
  m <- ggplot(dat.counts, aes(x = dim1, y = log2(count / spikeincounts))) + 
    geom_point() + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.counts, aes(x = dim1, y = log2(count / totalcounts))) + 
    geom_point() + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  
  # Load RNA data  ----------------------------------------------------------
  
  inf.rna <- file.path(hubprefix, "jyeung/data/public_data/TableS1.csv")
  
  dat.rna <- fread(inf.rna, skip = 1)
  
  # load GLMPCA ? 
  
  # pick straegy
  g2bin <- hash::hash(annot.out$out2.df$geneclean, annot.out$out2.df$region_coord)
  # genes.g2 <- subset(dat.rna, strategy_group == "B")$gene_symbol
  genes.g2 <- subset(dat.rna, strategy_group == "G")$gene_symbol
  
  
  # genes.g2 <- c("PLK1", "CENPE", "KIF11", "UBE2C", "TOP2A", "KIF2C", "CDK1", "TTK")
  bins.g2 <- sapply(genes.g2, function(x) AssignHash(x, g2bin))
  bins.g2 <- bins.g2[which(!is.na(bins.g2))]
  
  # bins.g2 <- sample(rownames(mat), size = 100)
  # bins.g2 <- names(loadings[50:100])
  
  dat.counts.merged <- data.frame(count = colSums(mat[bins.g2, ]), cell = colnames(mat), stringsAsFactors = FALSE) %>%
    left_join(., dat.umap)
  
  
  m <- ggplot(dat.counts.merged, aes(x = dim1, y = log2(count / spikeincounts))) + 
    geom_point() + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.counts.merged, aes(x = dim1, y = log2(count / totalcounts))) + 
    geom_point() + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.counts.merged, aes(x = pseudotime, y = log2(count / spikeincounts))) + 
    geom_point() + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  dev.off()
  
  # Fit each TSS counts  ----------------------------------------------------
  
  # set up metadata 
  rnames <- rownames(mat)
  names(rnames) <- rnames
  
  jfits <- lapply(rnames, function(rname){
    dat.counts.merged <- data.frame(count = mat[rname, ], cell = colnames(mat), stringsAsFactors = FALSE) %>%
      left_join(., dat.umap)
    jfit <- glm(formula = count ~ cellcycle.str + offset(log(spikeincounts)), family = "poisson", data = dat.counts.merged)
  })
  save(jfits, dat.umap, file = outrdata)
  
}, mc.cores = length(jmarks))

# for (jmark in jmarks){
# }



