# Jake Yeung
# Date of Creation: 2021-02-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/14b-check_DE_pval_vs_topn.R
# description

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

library(forcats)

library(topicmodels)
library(JFuncs)
library(scchicFuncs)

options(scipen=0)

make.plots <- FALSE

pvalcutoff <- 1e-10
# padjcutoff <- 1e-9
padjcutoff <- 1e-50

# Constants ---------------------------------------------------------------

jkeeptop <- 150
jlow.in.k9 <- TRUE

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load MATs  --------------------------------------------------------------

# ctypes <- c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs")
ctypes <- c("Eryths", "Bcells", "NKs", "pDCs", "Granulocytes", "DCs", "Basophils", "HSPCs")
ctypes.k9me3 <- c("Eryths", "Bcells", "Granulocytes", "HSPCs")

dat.metas <- lapply(jmarks, function(jmark){
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final/umaps_final_get_marks.", jmark, ".metadata.2020-12-17.txt")
  fread(inf)
})

dat.metas <- lapply(jmarks, function(jmark){
  dat.metas.tmp <- dat.metas[[jmark]]
  if (jmark != "H3K9me3"){
    dat.metas.tmp$jrep2 <- sapply(dat.metas.tmp$jrep, function(x) ifelse(x == "rep1old", "zold", "anew"))
  } else {
    dat.metas.tmp$jrep2 <- sapply(dat.metas.tmp$jrep, function(x) ifelse(x != "rep1old", "zold", "anew"))
  }
  return(dat.metas.tmp)
})

jmetas.pretty.lst <- lapply(jmarks, function(jmark){
  jmeta <- dat.metas[[jmark]]
  if (jmark != "H3K9me3"){
    jmeta$cluster <- factor(jmeta$cluster, levels = ctypes)
  } else { 
    jmeta$cluster <- factor(jmeta$cluster, levels = ctypes.k9me3)
  }
  jmeta <- jmeta %>% arrange(cluster, jrep)
})

cells.keep.lst <- lapply(jmetas.pretty.lst, function(jdat){
  jdat$cell
})


# Get DE outputs ----------------------------------------------------------

jfits.lst.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jinf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.RData"))
  load(jinf, v=T)
  return(jfits.lst)
})


params.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfits.lst <- jfits.lst.lst[[jmark]]
  jnames <- names(jfits.lst); names(jnames) <- jnames
  params.dat.all <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xkeep <- grepl("^Cluster.*.Estimate$", x = names(x))
    jparams <- x[xkeep]
    data.frame(bin = jname, param = names(jparams), estimate = unlist(jparams), mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  if (jmark == "H3K9me3"){
    params.dat.all <- params.dat.all %>% 
      mutate(param = gsub("Eryth", "Eryths", param),
             param = gsub("Lymphoid", "Bcells", param))
  }
  # make params more readable
  params.dat.all$ctype <- params.dat.all$param
  params.dat.all$ctype <- gsub("Cluster", "", params.dat.all$ctype)
  params.dat.all$ctype <- gsub(".Estimate", "", params.dat.all$ctype)
  return(params.dat.all)
})


pvals.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfits.lst <- jfits.lst.lst[[jmark]]
  jnames <- names(jfits.lst); names(jnames) <- jnames
  lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, padj = p.adjust(xvec), stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
})

pvals.lst2 <- lapply(jfits.lst.lst$H3K9me3, function(x) x$pval)
k9.bins <- which(pvals.lst2 < pvalcutoff)

pvals.adj.lst <- lapply(pvals.lst, function(jdat){
  jdat$padj <- p.adjust(jdat$pval, method = "BH")
  return(jdat)
})

pfilt <- lapply(pvals.adj.lst, function(jdat){
  subset(jdat, bin %in% names(k9.bins))
})

# plot(-log10(pfilt$H3K9me3$pval), -log10(pfilt$H3K9me3$padj), log = "xy")

signif.bins.lst <- lapply(jmarks, function(jmark){
  jdat <- pvals.lst[[jmark]]
  if (jmark != "H3K9me3"){
    jdat.filt <- subset(jdat, padj < padjcutoff)
  } else {
    jdat.filt <- subset(jdat, pval < pvalcutoff)
  }
  coord <- jdat.filt$bin
  jstart <- sapply(coord, JFuncs::GetStart, returnAsInt = TRUE)
  jend <- sapply(coord, JFuncs::GetEnd, returnAsInt = TRUE)
  jchromo <- sapply(coord, JFuncs::GetChromo, add.chr = FALSE)
  jmid <- jstart + round((jstart - jend) / 2)
  # new coord
  coordnew <- paste(jchromo, paste(jmid - 1, jmid + 1, sep = "-"), sep = ":")
  
  jdat <- data.frame(Coord = coordnew, 
                     Coordold = coord,
                     Chromo = jchromo, 
                     Start = jstart, End = jend, 
                     Mid = jmid, 
                     stringsAsFactors = FALSE)  %>%
    ungroup() %>%
    mutate(Dist = End - Start,
           mark = jmark)
  
  
})


inf.lda.check <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM_dbl_reseq.varfilt.k4_k9_dynamic_bins/lda_outputs.count_name.H3K9me3.k4_k9_dynamic_bins.2021-01-30.K-30.binarize.FALSE/ldaOut.count_name.H3K9me3.k4_k9_dynamic_bins.2021-01-30.K-30.Robj")
load(inf.lda.check, v=T)

keepn <- 6085
bins.filt.lst <- lapply(pvals.lst, function(pvals.dat){
  pvals.dat[1:keepn, ]$bin
})

rnames.true <- sapply(rownames(count.mat)[!grepl(";NM", rownames(count.mat))], function(x) strsplit(x, ";")[[1]][[2]])
rnames.dynamic <- signif.bins.lst$H3K9me3$Coordold

length(intersect(rnames.true, rnames.dynamic))

jmark <- "H3K4me1"
inf.bedcheck1 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.", jmark, ".2021-01-30.bad.txt")
inf.bedcheck2 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_bins_all_marks_top_6085_dists_to_TSS.annot_table.", jmark, ".2021-02-15.bed")
inf.bedcheck3 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_bins_all_marks_top_6085_sorted_by_pval.", jmark, ".2021-02-19.bed")
bedcheck1 <- fread(inf.bedcheck1) %>%
  left_join(., subset(signif.bins.lst[[jmark]], select = c(Coord, Coordold)), by = c("region_coord" = "Coord"))
bedcheck2 <- fread(inf.bedcheck2)
bedcheck3 <- fread(inf.bedcheck3)

# length(intersect(bedcheck1$Coordold, rnames.true))
# length(intersect(bedcheck2$V4, rnames.true))
length(intersect(bedcheck2$V4, bedcheck1$Coordold))
length(intersect(bedcheck3$V4, bedcheck2$V4))

dat.raw.check <- read.table(paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps/table.", jmark, ".DE_bins_all_marks_top_6085_dists_to_TSS.annot_table.ward.D2.1000.2021-02-19.txt"), header = TRUE, row.names = 1)
dat.raw.check2 <- read.table(paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps/table.", jmark, ".dynamic_bins.50kb.ward.D2.1000000.2021-02-19.txt"), header = TRUE, row.names = 1)
rnames.raw <- sapply(rownames(dat.raw.check), function(x) strsplit(x, ";")[[1]][[2]], USE.NAMES = FALSE)
rnames.raw2 <- sapply(rownames(dat.raw.check2), function(x) strsplit(x, ";")[[1]][[2]], USE.NAMES = FALSE)

length(intersect(rnames.raw, bedcheck3$V4))
length(intersect(rnames.raw2, bedcheck3$V4))
length(intersect(rnames.raw2, bedcheck1$Coordold))

pvals.dat <- pvals.adj.lst[[jmark]]

jcheck.dat2 <- subset(pvals.dat, bin %in% rnames.raw2)
jcheck.dat1 <- subset(pvals.dat, bin %in% rnames.raw)

range(jcheck.dat2$pval)
range(jcheck.dat1$pval)

# signif.bins.lst$H3K9me3$