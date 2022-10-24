# Jake Yeung
# Date of Creation: 2022-07-29
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/11-plot_raw_counts_gene_lists.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(DescTools)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarks

jinf.tss <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databases/gene_tss/gene_tss_winsize.50000.bed")
assertthat::assert_that(file.exists(jinf.tss))

# Load new colors ---------------------------------------------------------

inf.colors.fixed <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/ctypes_on_umap_batch_corrected_colors_fixed/dat_colors_DC_monocyte_fixed.2022-05-17.txt"
dat.colors.fixed <- fread(inf.colors.fixed)


# Load meta ----------------------------------------------------------------

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt")
  dat.meta <- fread(inf.meta) %>%
    left_join(., dat.colors.fixed) %>%
    rowwise() %>%
    mutate(colcode = colcodenew)
  # replace colcode with colcodenew
})

dat.meta.colors <- subset(dat.meta.lst$k4me1, select = c(ctype.from.LL, colcode))
ctype2col <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)




# Load TSS mats -----------------------------------------------------------

count.mat.norm.lst <- lapply(jmarksold, function(jmark){
  print(jmark)
  # inf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_progs_only/lda_outputs.countmat_TSS_new_progenitors_only.", jmark, ".2022-07-27/ldaOut.countmat_TSS_new_progenitors_only.", jmark, ".2022-07-27.Robj")
  inf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_TSS_only/lda_outputs.count_mat_allmerged_for_LDA_TSS_only.", jmark, ".2022-07-20/ldaOut.count_mat_allmerged_for_LDA_TSS_only.", jmark, ".2022-07-20.Robj")
  load(inf, v=T)
  count.mat.norm <- sweep(count.mat, MARGIN = 2, STATS = colSums(count.mat), FUN = "/")
  return(count.mat.norm)
})



# Load gene lists  --------------------------------------------------------

marker.genes.c1 <- c("Tgm2", "Mmrn1", "Trim47", "Pdzk1ip1", "Mllt3", "Mecom", "Esam", "Cish", "Hes1", "Mycn")
marker.genes.c2 <- c("Dntt", "Gm5111", "Flt3", "9030619P08Rik", "Il1r1", "Vldlr")
marker.genes.c3 <- c("Uhrf1", "Rmr2", "Lig1", "Tipin", "Rad54l", "Fabp5", "Nrm", "Gatm", "Hn1l", "Fam111a")
marker.genes.mk <- c("Pf4", "Gata1", "Slc14a1", "F2r", "Itga2b", "Zfp385a", "Zfpm1", "Plek", "Cd9", "Zeb2")

marker.genes <- list(marker.genes.c1, marker.genes.c2, marker.genes.c3, marker.genes.mk)
names(marker.genes) <- c("c1", "c2", "c3", "mk")

indir.markers.mpp <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/gene_list_from_Carmago"

mpp.names <- c("MPP2", "MPP3", "MPP4"); names(mpp.names) <- mpp.names

marker.genes.mpp <- lapply(mpp.names, function(jname){
  readr::read_lines(file.path(indir.markers.mpp, paste0(jname, "_genes.txt")))
})

marker.genes.merged <- c(marker.genes, marker.genes.mpp)

# Write outputs -----------------------------------------------------------

jnames <- names(marker.genes.merged)
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/primetime_plots"
# jname <- "c2"
for (markref in jmarks){
  print(markref)
  outpdf <- file.path(outdir, paste0("marker_genes_merged.", markref, ".", Sys.Date(), ".pdf"))
  pdf(outpdf, useDingbats = FALSE)
  
  m.ctype <- ggplot(dat.meta.lst[[markref]], aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point() + 
    ggtitle(markref) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode, guide = "legend") +
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.ctype)
  for (jname in jnames){
    genes.vec <- marker.genes.merged[[jname]]
    print(length(genes.vec))
    genes.grep <- paste0(genes.vec, collapse = "|")
    count.mat <- count.mat.norm.lst[[markref]]
    count.mat.filt <- count.mat[grepl(genes.grep, rownames(count.mat)), ]
    print(dim(count.mat.filt))
    jsignal <- colMeans(count.mat.filt)
    # winsorize
    jsignal.win <- DescTools::Winsorize(jsignal, probs = c(0.05, 0.95))
    dat.signal <- data.frame(cell = names(jsignal), frac.counts = jsignal, frac.counts.win = jsignal.win, stringsAsFactors = FALSE) %>%
      left_join(., dat.meta.lst[[markref]])
    
    m <- ggplot(dat.signal, aes(x = umap1, y = umap2, color = log2(frac.counts.win * 1000 + 1))) + 
      geom_point() + 
      ggtitle(paste(markref, jname, "ngenes:", length(genes.vec))) + 
      scale_color_viridis_c() + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m)
    
    m.facet <- ggplot(dat.signal, aes(x = forcats::fct_reorder(.f = ctype.from.LL, .x = frac.counts.win * 1000 + 1, .fun = median, .desc = TRUE), 
                                      y = log2(frac.counts.win * 1000 + 1), 
                                      fill = colcode)) + 
      geom_boxplot() + 
      ggtitle(paste(markref, jname, "ngenes:", length(genes.vec))) + 
      scale_fill_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode, guide = "legend") +
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    print(m.facet)
  }
  dev.off()
}



