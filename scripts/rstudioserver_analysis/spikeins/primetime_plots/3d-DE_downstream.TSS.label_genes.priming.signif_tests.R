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


# PDFs --------------------------------------------------------------------

# mark1 <- "H3K4me3"; mark2vec <- c("H3K4me1", "H3K27me3")
mark1 <- "H3K4me1"; mark2vec <- c("H3K4me3", "H3K27me3")
# mark1 <- "H3K27me3"; mark2vec <- c("H3K4me1", "H3K4me3")

pdfout <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_TSS_genesets.refmark_", mark1, ".labels.", Sys.Date(), ".priming.with_ks_tests.K27me3_only.pdf")

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

# ggplot(dat.params.lst[[jmark]] %>% filter(abs(estimate) < 5), aes(x = estimate, fill = param)) +
#        ggtitle(jmark) + 
#        geom_density() + 
#        facet_wrap(~param) + 
#        theme_bw() + 
#        geom_vline(xintercept = 0, linetype = "dotted") + 
#        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.params.wide.lst <- lapply(dat.params.lst, function(jdat){
  reshape2::dcast(data = jdat, formula = bin ~ param + mark, value.var = "estimate")
})

dat.params.wide.joined <- Reduce(f = full_join, x = dat.params.wide.lst)

# fill NAs with zeros 
# dat.params.wide.joined[is.na(dat.params.wide.joined)] <- 0
# 
# ggplot(dat.params.wide.joined, aes(x = ClusterEryths.Estimate_H3K27me3, y = ClusterEryths.Estimate_H3K9me3)) + 
#   geom_point(alpha = 0.25) + 
#   theme_bw() + 
#   geom_density_2d() + 
#   geom_vline(xintercept = 0, linetype = "dotted") + 
#   geom_hline(yintercept = 0, linetype = "dotted") + 
#   coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(dat.params.wide.joined, aes(x = ClusterBcells.Estimate_H3K4me1, y = ClusterBcells.Estimate_H3K27me3)) + 
#   geom_point(alpha = 0.25) + 
#   theme_bw() + 
#   geom_density_2d() + 
#   geom_vline(xintercept = 0, linetype = "dotted") + 
#   geom_hline(yintercept = 0, linetype = "dotted") + 
#   coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# show priming? 

# jjset <- "Bcells"
jjsets <- unique(as.character(dat.genes$jset))

q25 <- function(x, p = 0.25){
  quantile(x, probs  = p)
}

jjset <- "NKs"
jmarktmp <- "H3K27me3"

# for (jmarktmp in jmarksall){
#   print(jmarktmp)
  for (jjset in jjsets){
    print(jjset)
    
    if (jjset == "HSPCs"){
      print("Skipping")
      next
    }
    
    
    jbins.tmp <- subset(dat.genes, jset == jjset)$gene
    jtmp <- dat.params.lst[[jmarktmp]] %>% filter(bin %in% jbins.tmp)
    print(dim(jtmp))
    jtmp$param <- gsub("Cluster", "", jtmp$param)
    jtmp$param <- gsub(".Estimate", "", jtmp$param)
    jtmp$jcol <- sapply(jtmp$param, function(x) cluster2col[[x]])
    jtmp <- jtmp %>%
      mutate(param = forcats::fct_reorder(.f = param, .x = estimate, .fun = q25, .desc = TRUE))
    
    jjset1 <- levels(jtmp$param)[[1]]
    jjset2 <- levels(jtmp$param)[length(levels(jtmp$param))]
    # check signif
    # ctypes.bg <- subset(jtmp, param != jjset & param != "Eryths")$estimate / log(2)
    ctypes.bg <- subset(jtmp, param != jjset2)$estimate / log(2)
    ctypes.fg <- -1 * subset(jtmp, param == jjset2)$estimate / log(2)
    ks.out <- ks.test(ctypes.bg, ctypes.fg, alternative = "less")
    D.out <- signif(ks.out$statistic, digits = 2)
    pval.out <- signif(ks.out$p.value, digits = 2)
      
    m <- ggplot(jtmp, aes(x = param, y = estimate / log(2), fill = jcol)) + 
      geom_point() + 
      scale_fill_identity() + 
      coord_cartesian(ylim = c(-5, 5)) + 
      ggtitle(jjset) + 
      geom_boxplot() + 
      geom_hline(yintercept = 0, linetype = "dotted") + 
      theme_bw() + 
      xlab("") + 
      ylab("log2FC relative to HSPCs") + 
      ggtitle(paste(jmarktmp, jjset, "\n", jjset1, "vs", jjset2, D.out, pval.out)) + 
      theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    print(m)
    
    m.abs <- ggplot(jtmp, aes(x = param, y = abs(estimate / log(2)), fill = jcol)) + 
      geom_point() + 
      scale_fill_identity() + 
      coord_cartesian(ylim = c(0, 5)) + 
      ggtitle(jjset) + 
      geom_boxplot() + 
      geom_hline(yintercept = 0, linetype = "dotted") + 
      theme_bw() + 
      xlab("") + 
      ylab("log2FC relative to HSPCs") + 
      ggtitle(paste(jmarktmp, jjset, "\n", jjset1, "vs", jjset2, D.out, pval.out)) + 
      theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    print(m.abs)
    
  }
# }

if (make.plots){
  dev.off()
}

