# Jake Yeung
# Date of Creation: 2022-02-05
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/11-pseudotime_downstream.R
# 

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

jinf.tss <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databases/gene_tss/gene_tss_winsize.50000.bed"
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

# Define marks ------------------------------------------------------------


jmark <- "k4me1"


# Define trajs ------------------------------------------------------------


hpcs <- c("HSCs", "LT", "ST", "MPPs")
# granu tracj
ctypes.granus <- c(hpcs, "CMP", "GMP", "Granulocytes")
dims.granus <- c("dim1", "dim2")
ctypes.bcells <- c(hpcs, "Bcells")
dims.bcells <- c("dim7", "dim8")

ctypes.eryths <- c(hpcs, "MEP", "Eryths")
dims.eryths <- c("dim2", "dim3")

ctypes.pdcs <- c(hpcs, "pDCs")
dims.pdcs <- c("dim8", "dim9")

ctypes.dcs <- c(hpcs, "DCs")
dims.dcs <- c("dim1", "dim2")

ctypes.basos <- c(hpcs, "CMP", "GMP", "Basophils")
dims.basos <- c("dim1", "dim2")

ctypes.nks <- c(hpcs, "NKs")
dims.nks <- c("dim4", "dim5")

ctypes.monos <- c(hpcs, "CMP", "GMP", "Monocytes")
dims.monos <- c("dim1", "dim2")

ctypes.lst <- list(ctypes.granus, ctypes.bcells, ctypes.eryths, ctypes.pdcs, ctypes.dcs, ctypes.basos, ctypes.nks, ctypes.monos)
names(ctypes.lst) <- sapply(ctypes.lst, function(x) x[length(x)])

ctypes.names <- names(ctypes.lst); names(ctypes.names) <- ctypes.names

dims.lst  <- list(dims.granus, dims.bcells, dims.eryths, dims.pdcs, dims.dcs, dims.basos, dims.nks, dims.monos)
names(dims.lst) <- ctypes.names


# Load imputes ------------------------------------------------------------

inf.lda <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28/ldaOut.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28.Robj")
load(inf.lda, v=T)

tm <- posterior(out.lda)

mat.impute.log <- t(log2(tm$topics %*% tm$terms))

# Load metadata -----------------------------------------------------------

# inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping/metadata_reannotate_from_LLmat_dynamicbins.", jmark, ".2022-02-02.txt")
inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping/metadata_reannotate_from_LLmat_dynamicbins.", jmark, ".2022-02-03.txt")
dat.meta.LL <- fread(inf.meta)

inf.traj <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/trajs_outputs.", jmark, ".2022-02-03.rds")
dat.traj <- readRDS(inf.traj)

# Load raw counts ---------------------------------------------------------

inf.glmpca <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/varfilt/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_batch.szname_none.niter_1000.RData")
load(inf.glmpca, v=T)
count.mat.from.glmpca <- glm.inits$Y.filt

dat.totalcounts <- data.frame(cell = colnames(count.mat.from.glmpca), totalcounts = colSums(count.mat.from.glmpca), stringsAsFactors = FALSE)

dat.glmpca <- data.frame(cell = rownames(glm.out$factors), glm.out$factors, stringsAsFactors = FALSE) %>%
  left_join(., dat.meta.LL)



# Load fits ---------------------------------------------------------------


topn <- 20
indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/ptimefits"
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/ptimefits_plots_top_hits"
dir.create(outdir)
for (jctype in ctypes.names){
  print(jctype)
  
  ctypes.tmp <- ctypes.lst[[jctype]]
  dims.tmp <- dims.lst[[jctype]]
  jdim1 <- dims.tmp[[1]]
  jdim2 <- dims.tmp[[2]]
  
  outpdf.tmp <- file.path(outdir, paste0("plot_top_ptime_hits.", jctype, ".", jmark, ".", Sys.Date(), ".pdf")) 
  outf.tmp <- file.path(outdir, paste0("fits_with_annotations.", jctype, ".", jmark, ".", Sys.Date(), ".txt"))
  
  inf.fits <- file.path(indir, paste0("ptime_glm_fits_parallel.", jctype, ".", jmark, ".2022-02-04.rds"))
  dat.fits <- readRDS(inf.fits)
  
  genes <- names(dat.fits); names(genes) <- genes
  dat.fits.long <- lapply(genes, function(gene){
    dat.coefs <- dat.fits[[gene]]$coefs1 %>%
      mutate(gene = gene)
  }) %>%
    bind_rows()
  
  
  dat.traj.ctype <- dat.traj[[jctype]] %>%
    filter(!is.na(ptime)) 
  
  
  # Plot examples  ----------------------------------------------------------
  
  pdf(outpdf.tmp, useDingbats = FALSE)
  
  m.volcano <- ggplot(dat.fits.long, aes(x = estimate, y = -log10(pval10))) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.volcano)
  
  # jgene <- dat.fits.long[order(dat.fits.long$pval10), ]$gene[[1]]
  
  # bins.annot <- AnnotateCoordsFromList.GeneWise(coords.vec = dat.fits.long$gene, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb ="org.Mm.eg.db", chromos.keep = jchromos)
  bins.annot <- AnnotateCoordsFromList(coords.vec = dat.fits.long$gene, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb ="org.Mm.eg.db", chromos.keep = jchromos)
  
  topbins <- (dat.fits.long %>% arrange(desc(estimate), pval10))$gene[1:topn]
  names(topbins) <- topbins
  
  # integrate fits with bin annots
  
  dat.fits.long.annot <- left_join(dat.fits.long %>% dplyr::rename(region_coord = gene), 
                                   bins.annot$out2.df.closest %>% dplyr::select(region_coord, gene, tssname, dist.to.tss, seqnames, start, end, midpt.1)) %>%
    dplyr::arrange(desc(estimate), pval10)
 
  # write output
  fwrite(dat.fits.long.annot, file = outf.tmp, sep = "\t")
  
  topbins.genes <- sapply(topbins, function(b){
    jsub <-  subset(bins.annot$out2.df.closest %>% arrange(abs(dist.to.tss)), region_coord == b)
    if (nrow(jsub) > 0){
      (jgene <- jsub$gene)[[1]]
    } else {
      jgene <- "NotFound"
    }
    return(jgene)
  })
  
  m <- ggplot(dat.traj[[jctype]], aes_string(x = jdim1, y = jdim2, color = "ptime")) + 
    geom_point() + 
    scale_color_viridis_c(na.value = "grey85") + 
    theme_bw() + 
    ggtitle(paste(jmark, jctype)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  for (jbin in topbins){
    jgene <- topbins.genes[[jbin]]
    
    
    jrow <- count.mat.from.glmpca[jbin, dat.traj.ctype$cell]
    jrow.impute <- mat.impute.log[jbin, dat.traj.ctype$cell]
    
    dat.counts.row <- data.frame(cell = dat.traj.ctype$cell, bincounts = jrow, impute = jrow.impute, stringsAsFactors = FALSE) %>%
      left_join(dat.totalcounts) %>%
      left_join(dat.traj.ctype) %>%
      ungroup() %>%
      mutate(logbinfracs = log(bincounts / totalcounts + 1), 
             logbinfracs.wins = DescTools::Winsorize(x = logbinfracs, probs = c(0.01, 0.99)))
    
    jtitle <- paste(jbin, jgene)
    
    m <- ggplot(dat.counts.row, aes_string(x = jdim1, y = jdim2, color = "logbinfracs.wins")) + 
      geom_point() + 
      ggtitle(jbin) + 
      scale_color_viridis_c(na.value = "grey85") + 
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
    
    m <- ggplot(dat.counts.row %>% arrange(ctype.from.LL), aes(x = ptime, y = impute, color = ctype.from.LL)) + 
      geom_point(alpha = 0.75) + 
      ggtitle(jtitle) + 
      scale_color_viridis_d(na.value = "grey85") + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
    
  }
  dev.off()
  
  
}


