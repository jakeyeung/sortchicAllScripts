# Jake Yeung
# Date of Creation: 2020-12-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/12-summarize_bins_by_DE_parameters.with_mRNA.R
# Make panelsto summarize Giladi and DE of bins

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(JFuncs)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

library(topicmodels)
library(DescTools)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)



# Constants ---------------------------------------------------------------



jkeeptop <- 150
# jlow.in.k9 <- FALSE
jlow.in.k9 <- TRUE


outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/H3K9me3_analysis"
dir.create(outdir)
outpdf <- file.path(outdir, paste0("heatmap_allmarks_signif_bins_k9.lowink9_", jlow.in.k9, ".", Sys.Date(), ".AllInOne.pdf"))
outprefix <- file.path(outdir, paste0("heatmap_allmarks_signif_bins_k9.lowink9_", jlow.in.k9, ".", Sys.Date(), ".AllInOne.params"))

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load MATs  --------------------------------------------------------------

ctypes <- c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs")
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



inf.mat.adj <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_H3K4me1_H3K9me3_analysis/batch_corrected_imputed_values.bins.all_marks.mat.namesfix.2020-12-20.H3K27me3rep2rep3reseq.RData"
load(inf.mat.adj, v=T)

mat.adj.lst <- lapply(mat.adj.lst, function(jmat){
  rownames(jmat) <- jmat$rname
  jmat$rname <- NULL
  return(jmat)
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
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
})

pvals.lst2 <- lapply(jfits.lst.lst$H3K9me3, function(x) x$pval)
k9.bins <- which(pvals.lst2 < 1e-10)


# Get k9 bins and plot  ---------------------------------------------------


pval.k9.sub <- subset(pvals.lst$H3K9me3, pval < 1e-10) %>%
  arrange(desc(pval))

k9.bins.names <- unique(pval.k9.sub$bin)
ctypes.keep <- c("Eryths", "Bcells", "Granulocytes")
params.keep <- paste("Cluster", ctypes.keep, ".Estimate", sep = "")


params.dat.wide.lst <- lapply(jmarks, function(jmark){
  jsub <- subset(params.lst[[jmark]], bin %in% k9.bins.names & param %in% params.keep) %>%
    group_by(bin) %>% filter(max(abs(estimate)) < 5)
  jdat <- GetParamsWideFormat(jsub, jvalue.var = "estimate")
  # # keep only effect cnames ( do this later maybe ? )
  # cnames.keep.i <- grep("effect$", colnames(jdat))
  # cnames.new <- paste(colnames(jdat)[cnames.keep.i], jmark, sep = "_")
  # colnames(jdat)[cnames.keep.i] <- cnames.new
  # cnames.keep.bin.i <- grep("bin", colnames(jdat))
  # cnames.keep.merged.i <- c(cnames.keep.bin.i, cnames.keep.i)
  # jdat.filt <- jdat[, cnames.keep.merged.i]
  return(jdat)
})

bins.keep.lst <- GetK9CelltypeBins(params.dat.wide.lst$H3K9me3, low.in.k9 = jlow.in.k9, keeptop = jkeeptop)
bnames <- names(bins.keep.lst); names(bnames) <- bnames

jbins.eryth <- bins.keep.lst[["Eryths"]]
jbins.bcell <- bins.keep.lst[["Bcells"]]
jbins.granu <- bins.keep.lst[["Granulocytes"]]
jbins.hspcs <- bins.keep.lst[["HSPCs"]]

bins.keep <- c(jbins.eryth, jbins.bcell, jbins.granu, jbins.hspcs)


bins.common <- Reduce(f = intersect, x = lapply(mat.adj.lst, rownames))
bins.keep.common <- bins.keep[bins.keep %in% bins.common]

# Get colors for bins  ----------------------------------------------------



ctype2col <- hash::hash(jmetas.pretty.lst$H3K4me1$cluster, jmetas.pretty.lst$H3K4me1$colorcode)
names(bins.keep) <- c(rep("Eryths", jkeeptop), rep("Bcells", jkeeptop), rep("Granulocytes", jkeeptop), rep("HSPCs", jkeeptop))
colsvec <- sapply(names(bins.keep), function(x) AssignHash(x, jhash = ctype2col, null.fill = NA))
bin2col <- hash::hash(bins.keep, colsvec)

# Giladi ------------------------------------------------------------------

inf.public <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/giladi_pseudobulk_exprs_data.rds"
dat.public <- readRDS(inf.public)

m2c <- MarkerToCelltype()
dat.public$celltype2 <- sapply(as.character(dat.public$celltype), function(x) m2c[[x]])

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
bins.uniq <- unique(bins.keep)
dat.annot <- AnnotateCoordsFromList(coords.vec = bins.uniq, inf.tss = "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.100000.bed", txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)


# Make subset matrices winsorize -------------------------------------------------


jmat.lst <- lapply(jmarks, function(jmark){
  jmat <- mat.adj.lst[[jmark]][bins.keep.common, cells.keep.lst[[jmark]]]
  jmat <- apply(jmat, 2, function(jcol) Winsorize(jcol, probs = c(0.01, 0.99)))
  jmat <- t(apply(jmat, 1, function(jrow) Winsorize(jrow, probs = c(0.01, 0.99))))
  return(jmat)
})



# Make plots  -------------------------------------------------------------




print("Making heatmaps")
pdf(outpdf, useDingbats = FALSE)
for (jmark in jmarks){
  heatmap3::heatmap3(jmat.lst[[jmark]], Rowv = NA, Colv = NA, scale = "row", ColSideColors = jmetas.pretty.lst[[jmark]]$clustercol, RowSideColors = sapply(bins.keep.common, AssignHash, jhash = bin2col),  revC = TRUE, main = paste0(jmark, " 50kb bins"), margins = c(5, 8))
}

print("Plotting DE parameters")

# make boxplots
for (jmark in jmarks){
  for (jset in bnames){
    m.boxplot <- ggplot(params.lst[[jmark]] %>% 
                          filter(bin %in% bins.keep.lst[[jset]]) %>%
                          group_by(bin) %>% 
                          filter(max(abs(estimate)) < 5), 
                        aes(x = forcats::fct_reorder(.f = ctype, .x = estimate, .fun = median, .desc = TRUE), y = estimate)) + 
      ggtitle(jset, jmark) + 
      ylab("log2FC relative to HSPCs") + 
      geom_boxplot() + 
      geom_point() + 
      theme_bw(24) + 
      geom_hline(yintercept = 0, linetype = "dotted") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    print(m.boxplot)
  }
}

# make boxplots
for (jmark in jmarks){
  for (jset in bnames){
    m.boxplot <- ggplot(params.lst[[jmark]] %>% 
                          filter(bin %in% bins.keep.lst[[jset]]) %>%
                          group_by(bin) %>% 
                          filter(max(abs(estimate)) < 5), 
                        aes(x = ctype, y = estimate)) + 
      ggtitle(jset, jmark) + 
      ylab("log2FC relative to HSPCs") + 
      geom_boxplot() + 
      geom_point() + 
      theme_bw(24) + 
      geom_hline(yintercept = 0, linetype = "dotted") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    print(m.boxplot)
  }
}




print("Correlating with Giladi")

dat.params <- lapply(bnames, function(bname){
  jbins <- bins.keep.lst[[bname]]
  jgenes <- subset(dat.annot$out2.df.closest, region_coord %in% jbins)$gene
  dat.public.sub <- subset(dat.public, gene %in% jgenes) %>%
    ungroup() %>%
    mutate(celltype2 = forcats::fct_reorder(.f = celltype2, .x = zscore, .fun = median, .desc = TRUE),
           celltype = forcats::fct_reorder(.f = celltype, .x = zscore, .fun = median, .desc = TRUE))
  m <- ggplot(dat.public.sub, aes(x = interaction(celltype, celltype2), y = zscore)) + 
    ggtitle(paste(bname, "Ngenes:", length(jgenes))) +  
    geom_boxplot(outlier.shape = NA) + 
    geom_point(alpha = 0.25) + 
    theme_bw(24) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
  m <- ggplot(dat.public.sub, aes(x = interaction(celltype, celltype2), y = exprs)) + 
    ggtitle(paste(bname, "Ngenes:", length(jgenes))) +  
    geom_boxplot(outlier.shape = NA) + 
    geom_point(alpha = 0.25) + 
    theme_bw(24) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
  # get annotated params
  params.sub <- subset(params.dat.wide.lst$H3K9me3, bin %in% jbins) %>%
    mutate(label = bname, "LowInK9me3" = jlow.in.k9)
  return(params.sub)
}) %>%
  bind_rows()

dat.params <- lapply(bnames, function(bname){
  jbins <- bins.keep.lst[[bname]]
  jgenes <- subset(dat.annot$out2.df.closest, region_coord %in% jbins)$gene
  dat.public.sub <- subset(dat.public, gene %in% jgenes) %>%
    ungroup() %>%
    mutate(celltype2 = forcats::fct_reorder(.f = celltype2, .x = zscore, .fun = median, .desc = TRUE),
           celltype = forcats::fct_reorder(.f = celltype, .x = zscore, .fun = median, .desc = TRUE))
  
  m <- ggplot(dat.public.sub, aes(x = celltype, y = zscore)) + 
    ggtitle(paste(bname, "Ngenes:", length(jgenes))) +  
    geom_boxplot(outlier.shape = NA) + 
    geom_point(alpha = 0.25) + 
    theme_bw(24) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
  m <- ggplot(dat.public.sub, aes(x = celltype, y = exprs)) + 
    ggtitle(paste(bname, "Ngenes:", length(jgenes))) +  
    geom_boxplot(outlier.shape = NA) + 
    geom_point(alpha = 0.25) + 
    theme_bw(24) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
  # get annotated params
  params.sub <- subset(params.dat.wide.lst$H3K9me3, bin %in% jbins) %>%
    mutate(label = bname, "LowInK9me3" = jlow.in.k9)
  return(params.sub)
}) %>%
  bind_rows()



dev.off()



# Write list  -------------------------------------------------------------

# matrix, and gene annotations

for (jmark in jmarks){
  outftxt <- paste0(outprefix, ".", jmark, ".txt")
  jtmp <- left_join(params.dat.wide.lst[[jmark]], pvals.lst[[jmark]], by = "bin")
  fwrite(x = jtmp, file = outftxt, sep = "\t")
}


