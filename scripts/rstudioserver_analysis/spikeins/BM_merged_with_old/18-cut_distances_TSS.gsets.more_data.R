# Jake Yeung
# Date of Creation: 2020-12-06
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/18-cut_distances_TSS.gsets.more_data.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(zoo)


# Load DS cuts ------------------------------------------------------------


hubprefix <- "/home/jyeung/hub_oudenaarden"

# cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
jcol1 <- "#0072B2"
jcol1 <- "#D55E00"
jcoltxt <- gsub("#", "", jcol1)
cbPalette <- c(jcol1, "gray65", "#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

jsmooth <- 35
# jrad <- 500
jrad <- 2000
# inname <- "H3K9me3-BM_AllMerged.HSCs_topic22.sorted.cleaned.TSS_cuts.mat.gz"

# "/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/DedupR1onlyNoAltHits/cut_positions_TSS_rad_500_pos"

gsets <- c("Neutrophil", "HSCs", "Bcell", "Erythroblast")

# jpairs <- list(c("Erythroblasts", "Erythroblast"), c("Granulocytes", "Neutrophil"), c("HSPCs", "HSCs"), c("Bcells", "Bcell"))
jpairs <- list(c("Eryths", "Erythroblast"), c("Granulocytes", "Neutrophil"), c("HSPCs", "HSCs"), c("Bcells", "Bcell"))

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
jstrands <- c("pos", "neg"); names(jstrands) <- jstrands

# ntopics <- "1000"
indir.geneset <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations/MouseBMFromRNAseq")
# outdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/pdfs_from_WKM_BM_merged/cut_distances_gsets.radius_", jrad, ".BM.2020-08-10")
outdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/cut_distances_gsets_BM_all_merged")
dir.create(outdir)

for (jmark in jmarks){
  
  
  for (jpair in jpairs){
    jctype <- jpair[[1]]
    gset.ref <- jpair[[2]]
    
    
    outname <- paste0(paste(jmark, jctype, jcoltxt, sep = "_"), ".pdf")
    outpdf <- file.path(outdir, outname)
    
    jtitle <- paste(jmark, jctype)
    
    # jstrand <- "pos"
    cutsout <- lapply(jstrands, function(jstrand){
      inname <- paste0("BM_round1_round2_merged_", jmark, "_", jctype, ".TSS_cuts.radius_", jrad, ".", jstrand, ".mat.gz")
      insum <- gsub("mat.gz", "summary.gz", inname)
      # indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/DedupR1onlyNoAltHits/cut_positions_TSS_rad_", jrad, "_" , jstrand))
      indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds/cut_positions/cut_positions_TSS_rad_", jrad, "_", jstrand))
      assertthat::assert_that(dir.exists(indir))
      inf <- file.path(indir, inname)
      infsum <- file.path(indir, insum)
      assertthat::assert_that(file.exists(inf))
      assertthat::assert_that(file.exists(infsum))
      
      mat <- fread(inf)
      smry <- fread(infsum, col.names = c("chromo", "jstart", "jend", "jstrand", "jname"))
      smry$gene <- sapply(smry$jname, function(x) strsplit(x, "\\.\\.")[[1]][[2]])
      # add chromo? 
      smry$chromo <- paste("chr", smry$chromo, sep = "")
      return(list(mat = mat, smry = smry))
    })
    
    mat <- rbind(as.matrix(cutsout$pos$mat), as.matrix(cutsout$neg$mat[, ncol(cutsout$neg$mat):1]))
    smry <- rbind(cutsout$pos$smry, cutsout$neg$smry)
    
    
    dat.smooth.long <- lapply(gsets, function(jgset){
      print(jgset)
      # inf.geneset <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations/MouseBMFromTopics.1000/MouseBM_TSS_FromTopics.", jgset, ".bsize_2.bed")
      # inf.geneset <- file.path(indir.geneset, paste0("MouseBM_TSS_FromTopics.", jgset, ".bsize_2.bed"))
      inf.geneset <- file.path(indir.geneset, paste0("MouseBM_TSS.", jgset, ".bsize_2.bed"))
      annot <- fread(inf.geneset, col.names = c("chromo", "jstart", "jend", "gene"))
      genes <- annot$gene
      # remove Ig* genes
      genes <- genes[sapply(genes, function(g) !startsWith(g, prefix = "Ig"))]
      indx <- which(smry$gene %in% genes)
      mat.sub <- mat[indx, ]
      smoothcounts <- rollapply(colMeans(mat.sub), jsmooth, sum)
      dat.smooth <- data.frame(x = seq(length(smoothcounts)), 
                               NormCounts = smoothcounts, 
                               geneset = jgset, 
                               stringsAsFactors = FALSE)
      return(dat.smooth)
    }) %>%
      bind_rows() %>% 
      rowwise() %>%
      mutate(geneset.ref = ifelse(geneset == gset.ref, gset.ref, paste0("zNot", gset.ref))) %>%
      mutate(x = x - jrad)
    
    # plot all
    dat.all <- data.frame(x = seq(ncol(mat)) - jrad, NormCounts = colMeans(mat), stringsAsFactors = FALSE)
    
    m.all <- ggplot(dat.all, aes(x = x, y = NormCounts)) + 
      geom_point() + 
      theme_bw() + theme(aspect.ratio=0.3, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle(jtitle, "no smoothing") +  
      xlab("Position from TSS") + 
      geom_vline(xintercept = 0, linetype = "dotted")
    
    
    m.gsets <- ggplot(dat.smooth.long, aes(x = x, y = NormCounts, color = geneset)) + 
      geom_point(alpha = 0.25) + 
      theme_bw() + theme(aspect.ratio=0.3, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle(jtitle, paste0("Smoothing Window:", jsmooth)) + 
      scale_color_manual(values = cbPalette) + 
      xlab("Position from TSS") + 
      geom_vline(xintercept = 0, linetype = "dotted")
    
    
    m.gsets.collapsed <- ggplot(dat.smooth.long %>% group_by(x, geneset.ref) %>% mutate(NormCounts = mean(NormCounts)), aes(x = x, y = NormCounts, color = geneset.ref)) + 
      geom_point(alpha = 0.5) + 
      theme_bw() + theme(aspect.ratio=0.3, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle(jtitle, paste0("Smoothing Window:", jsmooth)) + 
      scale_color_manual(values = cbPalette) + 
      xlab("Position from TSS") + 
      geom_vline(xintercept = 0, linetype = "dotted")
    
    pdf(file = outpdf, useDingbats = FALSE)
    print(m.all)
    print(m.gsets)
    print(m.gsets.collapsed)
    dev.off()
    
    
    
  }
  
  
}


