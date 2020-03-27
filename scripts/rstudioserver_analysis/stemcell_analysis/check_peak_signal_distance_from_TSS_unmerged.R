# Jake Yeung
# Date of Creation: 2020-03-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/1-check_peak_signal.R
# Pseudobulk calling peaks. What is the genome distribution? 

rm(list=ls())

jstart <- Sys.Date()

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)
library(ggrepel)


outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/peaks_analysis.unmerged"
dir.create(outdir)
assertthat::assert_that(dir.exists(outdir))

# fname <- paste0("peaks_counts_vs_total.", jmark, ".pdf")
# outf <- file.path(outdir, fname)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
# jmark <- "H3K4me3"






dat.annot.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  # inf.bed <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output/hd_merged.", jmark, ".minlength_1000/merged.", jmark, ".minlength_1000.cutoff_analysis.merged.withchr.annotated.bed")
  inf.bed <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output/hd_merged.", jmark, ".minlength_1000/merged.", jmark, ".minlength_1000.cutoff_analysis.annotated.bed")
  assertthat::assert_that(file.exists(inf.bed))
  dat.annot <- fread(inf.bed)
  colnames(dat.annot) <- c("chromo", "start", "end", "gene", "dist")
  plot(density(dat.annot$dist), main = jmark)
  plot(density(log10(abs(dat.annot$dist))), main = jmark)
  dat.annot$mark <- jmark
  dat.annot$coord <- paste(dat.annot$chromo, paste(dat.annot$start, dat.annot$end, sep = "-"), sep = ":")
  return(dat.annot)
})

dat.annot.long <- dat.annot.lst %>%
  bind_rows()

ggplot(dat.annot.long, aes(x = abs(log10(dist)))) + geom_density() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark, ncol = 1, scales = "free_y")

ggplot(dat.annot.long, aes(x = abs(dist))) + geom_density() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark, ncol = 1, scales = "free_y")

dat.annot.sum <- dat.annot.long %>%
  group_by(mark) %>%
  summarise(dist = median(abs(dist)))

ggplot(dat.annot.long, aes(y = log10(abs(dist)), x = mark, group = mark)) + geom_boxplot() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


 
for (jmark in jmarks){
  

  # Check raw count ---------------------------------------------------------
  
  # inf.mat <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output/count_tables_from_hiddenDomains.rds_format/merged.", jmark, ".minlength_1000.cutoff_analysis.merged.withchr.annotated.rds")
  inf.mat <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output/count_tables_from_hiddenDomains_unmerged/", jmark, "-BM_AllMerged.merged_by_clusters_with_NAs.txt")
  # peak.count.mat <- readRDS(inf.mat)
  peak.count.mat <- ReadMatTSSFormat(inf.mat, as.sparse=TRUE, add.coord = TRUE)
  
  # Load UMAP annnots --------------------------------------------------------
  
  inf.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
  assertthat::assert_that(file.exists(inf.annot))
  load(inf.annot, v=T)
  
  
  # Where are most reads located?  ------------------------------------------
  
  peak.count.mat[1:5, 1:5]
  
  # order pseudobulks by distance? 
  
  jsplit <- split(dat.umap.glm.fillNAs, dat.umap.glm.fillNAs$cluster)
  
  cnames.keep.lst <- lapply(jsplit, function(x){
    x$cell
  })
  
  peak.pbulk.lst <- SumAcrossClusters(peak.count.mat, cnames.keep.lst)
  peak.pbulk <- do.call(cbind, peak.pbulk.lst)
  rownames(peak.pbulk) <- sapply(rownames(peak.pbulk), function(x) strsplit(x, ";")[[1]][[1]])
  
  # take a pseudobulk, sort by Neutorphils-specific peaks?
  
  jctypes <- colnames(peak.pbulk)
  
  print(jctypes)
  
  # jctype <- "Neutrophils_topic2"
  for (jctype in jctypes){
    
    
    inf.bed.ctype <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output/hd_clusters.", jmark, ".minlength_1000/", jmark, "-BM_AllMerged.", jctype, ".sorted.minlength_1000/", jmark, "-BM_AllMerged.", jctype, ".sorted.minlength_1000_analysis.bed")
    assertthat::assert_that(file.exists(inf.bed.ctype))
    
    dat.peaks.ctype <- fread(inf.bed.ctype)
    
    coord <- paste(dat.peaks.ctype$V1, paste(dat.peaks.ctype$V2, dat.peaks.ctype$V3, sep = "-"), sep = ":")
    coord.dat <- data.frame(coord = coord, stringsAsFactors = FALSE) %>%
      left_join(., dat.annot.lst[[jmark]])
    
    jcommon <- intersect(rownames(peak.pbulk), rownames(peak.pbulk))
    # coord.keep <- intersect(rownames(peak.pbulk), coord)  # some peaks got merged.. so these would be "specific" neutros?
    peak.pbulk.ctypefilt <- as.data.frame(peak.pbulk[jcommon, ])
    
    # Take neutros, sort by reads  --------------------------------------------
    
    dat.pbulk.ctype <-  data.frame(coord = rownames(peak.pbulk.ctypefilt), ctype.count = peak.pbulk.ctypefilt[[jctype]], stringsAsFactors = TRUE) %>%
      arrange(desc(ctype.count)) %>%
      left_join(., dat.annot.lst[[jmark]]) %>%
      ungroup() %>%
      mutate(jrank = rank(-ctype.count, ties.method = "random"),
             gene.label = ifelse(jrank <= 50, gene, NA))
    
    jtitle <- paste0(jctype, " pseudobulk,\ncounts under peak vs distance: ", jmark)
    
    m1 <- ggplot(dat.pbulk.ctype, aes(x = ctype.count, y = dist, label = gene.label)) + geom_point(alpha = 0.25) + theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
      geom_text_repel() + ggtitle(jtitle)
    
    m2 <- ggplot(dat.pbulk.ctype, aes(x = ctype.count, y = sqrt(abs(dist)), label = gene.label)) + geom_point(alpha = 0.25) + theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
      geom_text_repel() + ggtitle(jtitle)
    
    
    pdf(file.path(outdir, paste0("peaks_counts_by_distance.", jctype, ".", jmark, ".pdf")), useDingbats = FALSE)
      print(m1)
      print(m2)
    dev.off()
    print(Sys.time() - jstart)
  }
}

