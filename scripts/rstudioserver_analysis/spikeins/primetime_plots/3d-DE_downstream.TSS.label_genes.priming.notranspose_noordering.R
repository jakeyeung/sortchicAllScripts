# Jake Yeung
# Date of Creation: 2021-03-08
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/3d-DE_downstream.TSS.label_genes.priming.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(ggrepel)
library(JFuncs)
library(hash)

hubprefix <- "/home/jyeung/hub_oudenaarden"

# ctypes <- c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "DCs", "pDCs", "HSPCs")
ctypes <- c("Eryths", "Bcells", "NKs", "pDCs", "DCs", "Basophils", "Granulocytes")
ctypes.rename <- c("Eryths", "Bcells", "NKs", "pDCs", "DCs", "Baso/Eosino", "Neutrophils")
ctypes.hash <- hash::hash(ctypes, ctypes.rename)


# PDFs --------------------------------------------------------------------

# mark1 <- "H3K4me3"; mark2vec <- c("H3K4me1", "H3K27me3")
mark1 <- "H3K4me1"; mark2vec <- c("H3K4me3", "H3K27me3")
# mark1 <- "H3K27me3"; mark2vec <- c("H3K4me1", "H3K4me3")

pdfout <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_TSS_genesets.refmark_", mark1, ".labels.", Sys.Date(), ".priming.NoTransposeNoOrdering.pdf")

make.plots <- TRUE



# Load metadata -----------------------------------------------------------

jmarksall <- c(mark1, mark2vec)
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/shuffled_cells"
dat.metas <- lapply(jmarksall, function(jmarktmp){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.shuffled.", jmarktmp, ".2021-02-19.txt")
  fread(file.path(indir, fname))
})

cluster2col <- hash::hash(dat.metas[[1]]$cluster, dat.metas[[1]]$clustercol)

# Load gene annots --------------------------------------------------------

# inf.genes <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-09.from_LDA_topics.condensed.heatmap.famousgenes.keepn_400.refmark_H3K4me3.2020-12-09.txt")
# dat.genes <- fread(inf.genes)
dat.genes <- GetGeneSets()


# try Eryth-specific genes
jsets.vec <- as.character(unique(dat.genes$jset))
names(jsets.vec) <- jsets.vec

ctypespec.lst <- lapply(jsets.vec, function(jjset){
  subset(dat.genes, jset == jjset)$gene
})

ctypespec.hash <- hash::hash(ctypespec.lst)

stypecols <- c("grey", "blue", "red")


# Load data ---------------------------------------------------------------



jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

dat.params.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  # inf.rdata <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again/poisson_fit_bins.", jmark, ".2020-12-12.newannot2.witherrors.MoreBins.RData"))
  inf.rdata <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again/poisson_fit_TSS.", jmark, ".2020-12-12.newannot2.witherrors.TSS.RData"))
  print(inf.rdata)
  load(inf.rdata, v=T)
  jfits.lst1 <- jfits.lst
  
  params.lst1 <- lapply(jfits.lst1, function(x){
    xkeep <- grepl("^Cluster.*.Estimate$", x = names(x))
    x[xkeep]
  })

  # pvals.vec <- unlist(lapply(jfits.lst1, function(x) x$pval))
  pvals.dat <- data.frame(bin = names(jfits.lst1), pval = unlist(lapply(jfits.lst1, function(x) x$pval)), stringsAsFactors = FALSE)
  
  jnames <- names(jfits.lst1)
  names(jnames) <- jnames
  
  params.dat1 <- lapply(jnames, function(jname){
    jparams <- params.lst1[[jname]]
    if (is.null(jparams)){
      return(data.frame(NULL))
    }
    # assertthat::assert_that(!is.null(jparams))
    data.frame(bin = jname, param = names(jparams), estimate = unlist(jparams), mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows() %>%
    left_join(., pvals.dat)
  if (jmark == "H3K9me3"){
   params.dat1 <- params.dat1 %>% 
    mutate(param = gsub("Eryth", "Eryths", param),
           param = gsub("Lymphoid", "Bcells", param))
  }
  return(params.dat1)
  
})

# handle k27me3 rownames --------------------------------------------------

bins.all <- unique(dat.params.lst[["H3K4me1"]]$bin)
coords <- sapply(bins.all, function(b) paste("chr", strsplit(b, ";")[[1]][[1]], sep = ""))
coord2rname <- hash::hash(coords, bins.all)

dat.params.lst[["H3K27me3"]]$bin <- sapply(dat.params.lst[["H3K27me3"]]$bin, function(b) AssignHash(x = b, jhash = coord2rname, null.fill = b))


# Set mark  ---------------------------------------------------------------

if (make.plots){
  pdf(pdfout, useDingbats = FALSE)
}
                 
jmark <- mark1


dat.params.wide.lst <- lapply(dat.params.lst, function(jdat){
  reshape2::dcast(data = jdat, formula = bin ~ param + mark, value.var = "estimate")
})

dat.params.wide.joined <- Reduce(f = full_join, x = dat.params.wide.lst)

jjsets <- unique(as.character(dat.genes$jset))




for (jmarktmp in jmarksall){
  print(jmarktmp)
  for (jjset in jjsets){
    print(jjset)
    jbins.tmp <- subset(dat.genes, jset == jjset)$gene
    jtmp <- dat.params.lst[[jmarktmp]] %>% filter(bin %in% jbins.tmp)
    print(dim(jtmp))
    jtmp$param <- gsub("Cluster", "", jtmp$param)
    jtmp$param <- gsub(".Estimate", "", jtmp$param)
    jtmp$jcol <- sapply(jtmp$param, function(x) cluster2col[[x]])
    jtmp$param <- sapply(jtmp$param, function(x) ctypes.hash[[x]])
    jtmp$param <- factor(jtmp$param, levels = ctypes.rename)
    print(levels(jtmp$param))
    m <- ggplot(jtmp, 
                aes(x = param, 
                    y = estimate / log(2), 
                    fill = jcol)) + 
      geom_point() + 
      scale_fill_identity() + 
      ggtitle(jjset) + 
      geom_boxplot() + 
      theme_bw() + 
      xlab("") + 
      ylab("log2FC relative to HSPCs") + 
      ggtitle(paste(jmarktmp, jjset)) + 
      theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
      # coord_flip(ylim = c(-5, 5)) + 
      # scale_y_continuous(breaks = c(-5, 0, 5)) + 
      coord_cartesian(ylim = c(-5, 5)) + 
      scale_y_continuous(breaks = c(-5, 0, 5)) + 
      geom_hline(yintercept = 0, linetype = "dotted") 
    print(m)
  }
}

if (make.plots){
  dev.off()
}

