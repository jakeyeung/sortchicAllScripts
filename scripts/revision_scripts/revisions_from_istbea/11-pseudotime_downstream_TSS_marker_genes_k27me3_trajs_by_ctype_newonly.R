# Jake Yeung
# Date of Creation: 2022-02-20
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/11-pseudotime_downstream_TSS_marker_genes_k27me3_trajs_by_ctype_newonly.R
# description

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)

library(topicmodels)
library(DescTools)


# Load trajs --------------------------------------------------------------

inf <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/trajs_dims_ctypes_list_output.2022-02-08.rds"
annot.trajs <- readRDS(inf)

# Define marks ------------------------------------------------------------

# jmarks <- c("k27me3"); names(jmarks) <- jmarks
jmarks <- c("k27me3", "k4me1", "k4me3"); names(jmarks) <- jmarks


# Run ---------------------------------------------------------------------

markergenes <- c("Ly6c2", "Ly6g", "S100a8", "S100a2", "Chil3", "Sox6", "Tal1", "Gata1", "Ebf1", "Cd79a", "Cd79b", "Hoxa9", "Meis1", "Runx2", "Kit", "Hlf", "Erg", "Cd34", "Stat4", "Tcf7", "Cebpe")
markergenes <- c("Ly6c2", "Ly6g", "S100a8", "S100a2", "Chil3", "Sox6", "Tal1", "Gata1", "Ebf1", "Cd79a", "Cd79b", "Hoxa9", "Meis1", "Runx2", "Kit", "Hlf", "Erg", "Cd34", "Stat4", "Tcf7", "Cebpe")
# topn <- 20
indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/ptimefits_TSS"
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/ptimefits_plots_top_hits_TSS_markergenes_colcode_byctype_newonly"
dir.create(outdir)

for (jmark in jmarks){
  
  
  ctypes.lst <- annot.trajs[[jmark]]$ctypes.lst
  ctypes.names <- annot.trajs[[jmark]]$ctypes.names
  dims.lst <- annot.trajs[[jmark]]$dims.lst
  
  
  # Load imputes ------------------------------------------------------------
  
  inf.lda <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_TSS/ldaAnalysis_fripfilt_TSS/ldaOut.count_mat_TSS_combined.", jmark, ".2022-02-07.Robj")
  load(inf.lda, v=T)
  
  count.mat.from.lda <- count.mat
  tm <- posterior(out.lda)
  mat.impute.log <- t(log2(tm$topics %*% tm$terms))
  
  # Load metadata -----------------------------------------------------------
  
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping/rename/metadata_", jmark, ".txt")
  inf.meta <- file.path(paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up/metadata.", jmark, ".txt"))
  dat.meta.LL <- fread(inf.meta)
  
  # inf.traj <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/rename/trajs_outputs.", jmark, ".rds")
  inf.traj <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/by_ctype/trajs_outputs.", jmark, ".2022-02-20.rds")
  dat.traj <- readRDS(inf.traj)
  
  # inf.meta.color <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_with_colors.", jmark, ".2022-02-13.txt")
  # dat.meta.color <- fread(inf.meta.color) %>%
  #   dplyr::select(cell, colcode)
  
  # dat.meta.LL <- left_join(dat.meta.LL, dat.meta.color)
  
  
  
  # Load raw counts ---------------------------------------------------------
  
  # inf.glmpca <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/varfilt/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_batch.szname_none.niter_1000.RData")
  
  # Load fits ---------------------------------------------------------------
  
  for (jctype in ctypes.names){
    
    inf.glmpca <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/glmpca_outputs_split_by_trajs/glmpca.", jctype, ".", jmark, ".bincutoff_0.binskeep_0.platename_batch.szname_none.niter_1000.RData")
    load(inf.glmpca, v=T)
    dat.totalcounts <- data.frame(cell = colnames(count.mat.from.lda), totalcounts = colSums(count.mat.from.lda), stringsAsFactors = FALSE)
    dat.glmpca <- data.frame(cell = rownames(glm.out$factors), glm.out$factors, stringsAsFactors = FALSE) %>%
      left_join(., dat.meta.LL)
  
    print(jctype)
    
    ctypes.tmp <- ctypes.lst[[jctype]]
    dims.tmp <- dims.lst[[jctype]]
    jdim1 <- dims.tmp[[1]]
    jdim2 <- dims.tmp[[2]]
    
    outpdf.tmp <- file.path(outdir, paste0("plot_top_ptime_hits.", jctype, ".", jmark, ".", Sys.Date(), ".pdf")) 
    outf.tmp <- file.path(outdir, paste0("fits_with_annotations.", jctype, ".", jmark, ".", Sys.Date(), ".txt"))
    
    inf.fits <- file.path(indir, paste0("ptime_glm_fits_parallel.", jctype, ".", jmark, ".2022-02-08.rds"))
    dat.fits <- readRDS(inf.fits)
    
    genes <- names(dat.fits); names(genes) <- genes
    
    dat.fits.long <- lapply(genes, function(gene){
      dat.coefs <- dat.fits[[gene]]$coefs1 %>%
        mutate(rname = gene, 
               gene = strsplit(rname, split = ";")[[1]][[2]])
    }) %>%
      bind_rows() %>%
      ungroup() %>%
      arrange(pval10, desc(estimate)) %>%
      filter(abs(estimate) <= 5)
    
    dat.traj.ctype <- dat.traj[[jctype]] %>%
      filter(!is.na(ptime)) %>%
      left_join(., subset(dat.meta.LL, select = c(cell, colcode)))
    
    
    # Plot examples  ----------------------------------------------------------
    
    pdf(outpdf.tmp, useDingbats = FALSE)
    
    m.volcano <- ggplot(dat.fits.long, aes(x = estimate, y = -log10(pval10))) + 
      geom_point(alpha = 0.25) + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m.volcano)
    
    # jgene <- dat.fits.long[order(dat.fits.long$pval10), ]$gene[[1]]
    
    # bins.annot <- AnnotateCoordsFromList.GeneWise(coords.vec = dat.fits.long$gene, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb ="org.Mm.eg.db", chromos.keep = jchromos)
    # bins.annot <- AnnotateCoordsFromList(coords.vec = dat.fits.long$gene, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb ="org.Mm.eg.db", chromos.keep = jchromos)
    
    # get number of reads for each hit
    # jrname <- "7:23212331-23222331;NM_001166850.1..Vmn1r165"
    # frac.zeros.vec <- sapply(dat.fits.long$rname[1:100], function(jrname){
    #   frac.zeros <- nnzero(count.mat.from.lda[jrname, ]) / length(count.mat.from.lda[jrname, ])
    # })
    # dat.frac.zeros <- data.frame(rname = names(frac.zeros.vec), frac.nonzero = frac.zeros.vec, stringsAsFactors = FALSE)
    # 
    # dat.fits.long.check <- left_join(dat.frac.zeros, dat.fits.long)
    # 
    # ggplot(dat.fits.long.check, aes(x = estimate, y = frac.nonzero)) + 
    #   geom_point() + 
    #   theme_bw() + 
    #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    # dat.fits.long$genesymbol <- sapply(datfits.long$gene, function(x) strsplit(x, split = "\\.")[[1]][[3]])
    dat.fits.long$genesymbol <- sapply(dat.fits.long$gene, function(x) strsplit(x, split = "\\.")[[1]][[4]])
    
    # topbins.genes <- subset(dat.fits.long, genesymbol %in% markergenes)$gene
    topbins <- sort(unique(subset(dat.fits.long, genesymbol %in% markergenes)$rname))
    topbins.genes <- sapply(topbins, function(x) strsplit(x, ";")[[1]][[2]])
    
    # topbins <- (dat.fits.long %>% arrange(desc(estimate), pval10))$rname[1:topn]
    names(topbins) <- topbins
      
    # topbins.genes <- (dat.fits.long %>% arrange(desc(estimate), pval10))$gene[1:topn]
    names(topbins.genes) <- topbins
    
      
    # integrate fits with bin annots
    
    dat.fits.long.annot <- dat.fits.long %>% dplyr::rename(region_coord = rname) %>%
      dplyr::arrange(desc(estimate), pval10)
    
    # write output
    fwrite(dat.fits.long.annot, file = outf.tmp, sep = "\t")
    
    m <- ggplot(dat.traj[[jctype]], aes_string(x = jdim1, y = jdim2, color = "ptime")) + 
      geom_point() + 
      scale_color_viridis_c(na.value = "grey85") + 
      theme_bw() + 
      ggtitle(paste(jmark, jctype)) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
    
    cells.keep <- dat.meta.LL$cell
    cells.keep.i <- colnames(count.mat.from.lda) %in% cells.keep
    assertthat::assert_that(identical(colnames(count.mat.from.lda), colnames(mat.impute.log)))
    
    for (jbin in topbins){
      jgene <- topbins.genes[[jbin]]
      
      jrow <- count.mat.from.lda[jbin, cells.keep.i]
      jrow.impute <- mat.impute.log[jbin, cells.keep.i]
      
      cells.keep.after <- names(jrow)
      
      dat.counts.row <- data.frame(cell = cells.keep.after, bincounts = jrow, impute = jrow.impute, stringsAsFactors = FALSE) %>%
        left_join(dat.totalcounts) %>%
        left_join(dat.traj.ctype) %>%
        ungroup() %>%
        mutate(logbinfracs = log(bincounts / totalcounts + 1), 
               logbinfracs.wins = DescTools::Winsorize(x = logbinfracs, probs = c(0.01, 0.99))) %>%
        filter(batch == "New")
      
      jtitle <- paste(jbin, jgene)
      
      m <- ggplot(dat.counts.row, aes_string(x = jdim1, y = jdim2, color = "colcode")) + 
        geom_point() + 
        ggtitle(jbin) + 
        # scale_color_viridis_c(na.value = "grey85") + 
        scale_color_identity() + 
        ggtitle(jtitle) + 
        theme_bw() + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      print(m)
      
      m <- ggplot(dat.counts.row, aes_string(x = jdim1, y = jdim2, color = "impute")) + 
        geom_point() + 
        ggtitle(jtitle) + 
        scale_color_viridis_c(na.value = "grey85") + 
        theme_bw() + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      print(m)
      
      dat.counts.row$ctype.from.LL <- factor(dat.counts.row$ctype.from.LL, levels = ctypes.tmp)
      
      m <- ggplot(dat.counts.row %>% arrange(ctype.from.LL), aes(x = ptime, y = impute, color = colcode)) + 
        geom_point(alpha = 0.75) + 
        ggtitle(jtitle) + 
        # scale_color_viridis_d(na.value = "grey85") + 
        scale_color_identity() + 
        theme_bw() + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      print(m)
      
    }
    dev.off()
    
    
  }
  
  
  
  
}


