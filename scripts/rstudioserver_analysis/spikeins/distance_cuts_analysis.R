# Jake Yeung
# Date of Creation: 2020-08-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/distance_cuts_analysis.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(zoo)

hubprefix <- "/home/jyeung/hub_oudenaarden"

jwidth <- 35

# jmark <- "H3K27me3"
# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
jmarks <- c("H3K9me3"); names(jmarks) <- jmarks
jstrands <- c("pos", "neg"); names(jstrands) <- jstrands
cc.vec <- c("0_G1", "1_S", "2_G2M")

for (jmark in jmarks){
  print(jmark)
  
  pdfout <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/cut_distances.K562_cellcycle/K562_cellcycle.", jmark, ".pdf")
  if (file.exists(pdfout)){
    print(paste("Skipping:", pdfout))
    next
  }
  
  dat.smooth.lst <- lapply(cc.vec, function(cc){
    print(cc)
    cutsout <- lapply(jstrands, function(jstrand){
      print(jstrand)
      inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/bams_split_by_clusters/K562-EtOH-", jmark, "-G1-G2/cut_positions_TSS_rad_2000_", jstrand, ".K562_spikeins/K562-EtOH-", jmark, "-G1-G2.sorted.tagged.", cc, ".sorted.TSS_cuts.radius_2000.", jstrand, ".mat.gz"))
      print(inf)
      # inf.sum <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/bams_split_by_clusters/K562-EtOH-", jmark, "-G1-G2/cut_positions_TSS_rad_2000_neg.K562_spikeins/K562-EtOH-", jmark, "-G1-G2.sorted.tagged.", cc, ".sorted.TSS_cuts.radius_2000.", jstrand, ".summary.gz"))
      mat <- fread(inf)
      return(list(mat = mat))
    })
    
    mats <- rbind(as.matrix(cutsout$pos$mat), as.matrix(cutsout$neg$mat[, ncol(cutsout$neg$mat):1]))
    
    smoothcounts <- rollapply(colMeans(mats), jwidth, sum)
    
    jdist <- (ncol(mats) - 1) / 2
    dat.smooth <- data.frame(x = seq(from = -jdist, by = 1, length.out = length(smoothcounts)), 
                             SmoothedCounts = smoothcounts, 
                             cellcycle = cc,
                             stringsAsFactors = FALSE)
    return(dat.smooth)
  })
  
  dat.smooth <- dat.smooth.lst %>%
    bind_rows() 
  
  # ggplot(dat.smooth, aes(x = x, y = SmoothedCounts, color = cellcycle)) + 
  #   geom_point(alpha = 0.25) +  
  #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #   ggtitle(jmark) + facet_wrap(~cellcycle, nrow = 1) + 
  #   xlab("Distance from TSS") + 
  #   coord_cartesian(xlim = c(-500, 500))
  
  # Normalize by spikeins  --------------------------------------------------
  
  inf.spike <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/spikein_counts/spikein_counts_all.RData")
  assertthat::assert_that(file.exists(inf.spike))
  
  load(inf.spike, v=T)
  
  dat.annot <- AddCellCycleLabel.bydat(dat.spikeins.mat) %>%
    rowwise() %>%
    mutate(mark = strsplit(samp, split = "-")[[1]][[3]])
  
  dat.annot.sum <- dat.annot %>%
    filter(mark == jmark) %>%
    mutate(cellcycle.str = gsub("/", "", cellcycle.str)) %>%
    group_by(mark, cellcycle.str) %>%
    summarise(spikeincounts = sum(spikeincounts))
  
  dat.smooth.norm <- left_join(dat.smooth, dat.annot.sum, by = c("cellcycle" = "cellcycle.str"))
  
  
  pdf(file = pdfout)
  
  m1 <- ggplot(dat.smooth.norm, aes(x = x, y = log2(SmoothedCounts / spikeincounts), color = cellcycle)) + 
    geom_point(alpha = 0.25)  + 
    theme_bw() + theme(aspect.ratio=0.33, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
    ggtitle(jmark) + 
    xlab("Distance from TSS")
  
  m2 <- ggplot(dat.smooth.norm, aes(x = x, y = log2(SmoothedCounts / spikeincounts), color = cellcycle)) + 
    geom_point(alpha = 0.25)  + 
    theme_bw() + theme(aspect.ratio=0.33, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
    ggtitle(jmark) + 
    xlab("Distance from TSS") + 
    coord_cartesian(xlim = c(-500, 500))
  
  print(m1)
  print(m2)
  
  dev.off()
}



