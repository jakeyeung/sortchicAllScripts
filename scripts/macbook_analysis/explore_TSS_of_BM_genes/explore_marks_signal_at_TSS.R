# Jake Yeung
# Date of Creation: 2020-03-12
# File: ~/projects/scchic/scripts/macbook_analysis/explore_TSS_of_BM_genes/explore_H3K4me3_signal_at_TSS.R
# 

rm(list=ls())

# load libs ---------------------------------------------------------------


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(reticulate)

library(DescTools)

library(purrr)
library(ggrastr)
library(pracma)

library(hash)
library(JFuncs)

library(forcats)
library(ggrepel)

reticulate::source_python("/Users/yeung/projects/scchic/scripts/python_functions/parse_dictionary_text.py")

MarkerToCelltype <- function(){
  # get correspondance between marker to celltype
  x <- list("Car1" = "Erythroblast",
            "core" = "HSCs",
            "Vpreb1" = "Bcell",
            "Siglech" = "pDendritic",
            "Prg2" = "Eosinophil",
            "Gstm1" = "Neutrophil",
            "Ly86" = "Monocyte",
            "Ccl5" = "NKcell",
            "Prss34" = "Basophil",
            "Cd74" = "cDendritic",
            "Pf4" = "Megakaryocyte",
            "Fcrla" = "Bcell",
            "Fcnb" = "Neutrophil",
            "Hba.a2" = "Erythroblast",
            "Ltf" = "Neutrophil")
  return(x)
}
annot.lst <- MarkerToCelltype()

# Load annots -------------------------------------------------------------

inf.annot <- "/Users/yeung/data/scchic/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat.celltypes_filt.rds"
dat.annots <- readRDS(inf.annot)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3")
jdirs <- c("lessthan", "greaterthan")
# jmark <- "H3K4me1"
# jdir <- "lessthan"
for (jmark in jmarks){
  for (jdir in jdirs){
    print(jmark)
    print(jdir)
    
    if (jdir == "greaterthan"){
      logfcmin <- 0.5
    } else if (jdir == "lessthan"){
      logfcmin <- -0.5
    } else {
      warning("jdir must be greaterthan or lessthan:", jdir)
    }
    pvalmax <- 0.01
    outpdf <- paste0("/Users/yeung/data/scchic/pdfs/marker_genes_Giladi_TSS_signal/", jmark, "_GiladiMarkerGenes.", jdir, ".pvallmax_", pvalmax, ".logfcmin_", logfcmin, ".pdf")
    
    
    # Load bed ----------------------------------------------------------------
    
    inf.bed <- "/Users/yeung/data/databases/gene_tss/gene_tss/giladi_filtered.seurat/gene_tss_winsize.2.diff_exprs_Giladi_seurat.celltypes_filt.bed"
    dat.bed <- fread(inf.bed)
    
    coord.hash <- hash::hash(sapply(dat.bed$V4, function(x) strsplit(x, ";")[[1]][[1]]), sapply(dat.bed$V4, function(x) strsplit(x, ";")[[1]][[2]]))
    
    # Load data  --------------------------------------------------------------
    
    inf <- paste0("/Users/yeung/data/scchic/from_cluster/bigwig_outputs/merged_bams.deeptools_outputs.tss.MAPQ_40.dist_10000.allctypes_from_seurat.", jmark, ".bsize_100/computeMatrix.MAPQ_40.", jmark, ".gene_tss_winsize.2.diff_exprs_Giladi_seurat.celltypes_filt.tab.gz")
    assertthat::assert_that(file.exists(inf))
    
    meta <- inf2dic(inf)
    mat.all <- fread(inf, header = FALSE, skip = 1, sep = "\t", quote = "")
    
    # name columns 
    nbins <- (unique(meta$downstream) + unique(meta$upstream)) / unique(meta[["bin size"]])
    nsamps <- length(meta$sample_labels)
    
    # Stack matrix vertically  ------------------------------------------------
    
    mat <- mat.all[, -seq(6)]
    coords <- mat.all[, c(seq(6))]
    colnames(coords) <- c("chromo", "start", "end", "coord", "meta1", "meta2")
    assertthat::assert_that(ncol(mat) == nsamps * nbins)
    
    # split into a list of matrices, then rbind
    sampids <- ceiling(seq(ncol(mat)) / nbins)
    sampids.uniq <- as.character(unique(sort(sampids)))
    names(sampids.uniq) <- meta$sample_labels
    
    colids <- ( (seq(ncol(mat)) - 1) %% nbins ) + 1  # - 1 and +1 so first element is 1, last element is 10
    
    colnames(mat) <- as.character(sampids)
    
    mats.lst <- lapply(sampids.uniq, function(sampid){
      cols.i <- which(colnames(mat) == sampid)
      mat.sub <- mat[, ..cols.i]
      colnames(mat.sub) <- paste("bin", seq(nbins), sep = "")
      return(mat.sub) 
    })
    
    # make long dataframe of bed locations
    coords.lst <- lapply(sampids.uniq, function(sampid){
      jtmp <- coords[, c(1,2,3,4)]
      jtmp$sampid <- sampid
      return(jtmp)
    })  
    
    mats.long <- do.call(rbind, mats.lst)
    coords.long <- do.call(rbind, coords.lst)
    
    mat.long.merge <- cbind(mats.long, coords.long)
    
    
    # Make gene exprs ---------------------------------------------------------
    
    mats.lst.clean <- lapply(mats.lst, function(jmat){
      jmat <- as.matrix(jmat)
      # jmat.win <- Winsorize(jmat, minval = 0, maxval = quantile(jmat, 0.98))
      jmat.win <- jmat
      return(jmat.win)
    }) 
    
    samp.remove <- names(mats.lst.clean)[grepl("BM_AllMerged..sorted.100", names(mats.lst.clean))]
    
    for (jsamp in samp.remove){
      mats.lst.clean[[jsamp]] <- NULL
    }
    
    gene.exprs <- lapply(mats.lst.clean, function(jmat){
      rowMeans(jmat)
    })
    
    cpm.mat <- do.call(cbind, gene.exprs)
    
    rownames(cpm.mat) <- coords$coord
    
    cpm.mat.long <- data.frame(coord = rownames(cpm.mat), gene = sapply(rownames(cpm.mat), function(x) AssignHash(x, coord.hash, NA)), cpm.mat, stringsAsFactors = FALSE)
    cpm.mat.long <- tidyr::gather(cpm.mat.long, key = "pseudobulk", value = "cpm", -c(coord, gene)) %>%
      group_by(coord) %>%
      # mutate(cpm = log2(cpm * 10^6 + 1)) %>%
      mutate(zscore = scale(cpm, center = TRUE, scale = TRUE))
    
    
    
    # get diff exprs genes ----------------------------------------------------
    
    jclust <- "Ltf"
    jclust <- "Fcrla"
    
    
    jclusts <- unique(dat.annots$cluster)
    
    pdf(outpdf, useDingbats = FALSE)
    for (jclust in jclusts){
      print(jclust)
      
      if (jdir == "greaterthan"){
        jsub <- subset(dat.annots, cluster == jclust & p_val_adj <= pvalmax & avg_logFC >= logfcmin)
      }  else if (jdir == "lessthan"){
        jsub <- subset(dat.annots, cluster == jclust & p_val_adj <= pvalmax & avg_logFC <= logfcmin)
      }
      
      plot(density(jsub$avg_logFC), main = paste("AvgLogFC for chosen genes:", jclust, "N=", length(jsub$gene)))
      
      # select gene with most variability across clusters
      
      cpm.varfilt.long <- cpm.mat.long %>%
        group_by(coord, gene) %>% 
        summarise(jsd = sd(cpm)) %>%
        group_by(gene) %>%
        filter(jsd == max(jsd))
      
      coords.keep <- cpm.varfilt.long$coord
      
      cpm.long.filt <- subset(cpm.mat.long, gene %in% unique(jsub$gene) & coord %in% coords.keep) %>%
        ungroup()
      
      
      # shorten names
      cpm.long.filt$pseudobulk <- gsub(".sorted.100$", "", cpm.long.filt$pseudobulk)
      cpm.long.filt$pseudobulk <- gsub(".BM_AllMerged", "", cpm.long.filt$pseudobulk)
      
      cpm.mat.filt <- as.data.frame(tidyr::spread(cpm.long.filt %>% dplyr::select(-coord, -cpm), key = pseudobulk, value = zscore))
      rownames(cpm.mat.filt) <- cpm.mat.filt$gene
      cpm.mat.filt$gene <- NULL
      
      pca.out <- prcomp(t(cpm.mat.filt), center = FALSE, scale. = FALSE)
      
      # plot(pca.out)
      
      pca.dat <- data.frame(pseudobulk = rownames(pca.out$x), PC1 = pca.out$x[, 1], PC2 = pca.out$x[, 2], stringsAsFactors = FALSE) 
      
      jtitle <- paste("Marker from scRNAseq:", jclust, "(", annot.lst[[jclust]], ")\nMark:", jmark, "Direction", jdir)
      
      m.pca <- ggplot(pca.dat, aes(x = PC1, y = PC2, label = pseudobulk)) + geom_point() + geom_text_repel() + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
        ggtitle(jtitle) 
      print(m.pca)
      
      # Show boxplot ------------------------------------------------------------
      
      m <- ggplot(cpm.long.filt, aes(x = forcats::fct_reorder(.f = pseudobulk, .x = zscore, .fun = median, .desc = TRUE), y = zscore)) + geom_boxplot() + 
        theme_bw() +
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
        xlab("") + 
        ggtitle(jtitle)
      print(m)
      m.exprs <- ggplot(cpm.long.filt, aes(x = forcats::fct_reorder(.f = pseudobulk, .x = cpm, .fun = median, .desc = TRUE), y = cpm)) + geom_boxplot() + 
        theme_bw() +
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
        xlab("") + 
        ggtitle(jtitle)
      print(m.exprs)
    }
    dev.off()
  }
}

dev.off()