# Jake Yeung
# Date of Creation: 2020-03-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/create_bed_annotated_regions.R
# Create a list of TSS, enhancer, etc and make count table out of it.

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

winsizes <- c(500L, 1000L, 5000L, 10000L)

# winsize <- 500L  # 5kb left and right
for (winsize in winsizes){
  outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/ENCODE"
  inf.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/ENCODE/genomic_regions_10kbpromoters.txt"
  dat.annot <- fread(inf.annot) %>%
    mutate(start = as.integer(start),
           stop = as.integer(stop))
  
  
  proms <- subset(dat.annot, type == "promotor") 
  
  enhs <- subset(dat.annot, type == "enhancer")
  
  genes <- subset(dat.annot, type == "gene")
  
  
  # Writer promoters in a strand specific way for deeptools ------------------------------
  
  proms.pos <- subset(proms, strand == "+") %>%
    rowwise() %>%
    mutate(chr2 = paste("chr", chr, sep = ""),
           name2 = paste(name, strand, sep = "")) %>%
    dplyr::select(chr2, start, stop, name2)
  
  proms.neg <- subset(proms, strand == "-") %>%
    rowwise() %>%
    mutate(chr2 = paste("chr", chr, sep = ""),
           name2 = paste(name, strand, sep = "")) %>%
    dplyr::select(chr2, start, stop, name2)
  
  enhs.ref <- enhs %>%
    rowwise() %>%
    mutate(chr2 = paste("chr", chr, sep = ""),
           name2 = paste(name, strand, sep = "")) %>%
    dplyr::select(chr2, start, stop, name2)
  
  # bed file should be chromo, start, stop, geneid
  fwrite(proms.pos, file = file.path(outdir, paste0("promotersPos_winsize_0.", Sys.Date(), ".bed")), sep = "\t", col.names = FALSE, row.names = FALSE)
  fwrite(proms.neg, file = file.path(outdir, paste0("promotersNeg_winsize_0.", Sys.Date(), ".bed")), sep = "\t", col.names = FALSE, row.names = FALSE)
  fwrite(enhs.ref, file = file.path(outdir, paste0("enhancers_winsize_0.", Sys.Date(), ".bed")), sep = "\t", col.names = FALSE, row.names = FALSE)
  
  
  # Writer promoters with window  -------------------------------------------
  
  proms.window <- proms %>%
    rowwise() %>%
    mutate(start = as.integer(start - winsize / 2),
           stop = as.integer(stop + winsize / 2),
           chr2 = paste("chr", chr, sep = ""),
           name2 = paste(name, strand, sep = "")) %>%
    
    dplyr::select(chr2, start, stop, name2)
  
  enhs.window <- enhs %>%
    rowwise() %>%
    mutate(start = as.integer(start - winsize / 2),
           stop = as.integer(stop + winsize / 2),
           chr2 = paste("chr", chr, sep = ""),
           name2 = paste(name, strand, sep = "")) %>%
    dplyr::select(chr2, start, stop, name2)
  
  promsenhs.window <- rbind(proms.window, enhs.window) 
  
  genes.window <- genes %>%
    rowwise() %>%
    mutate(start = as.integer(start),
           stop = as.integer(stop),
           chr2 = paste("chr", chr, sep = ""),
           name2 = paste(name, strand, sep = "")) %>%
    dplyr::select(chr2, start, stop, name2)
  
  
  fwrite(proms.window, file = file.path(outdir, paste0("promoters_winsize_", winsize, ".", Sys.Date(), ".bed")), sep = "\t", col.names = FALSE, row.names = FALSE)
  fwrite(enhs.window, file = file.path(outdir, paste0("enhancers_winsize_", winsize, ".", Sys.Date(), ".bed")), sep = "\t", col.names = FALSE, row.names = FALSE)
  fwrite(promsenhs.window, file = file.path(outdir, paste0("promsenhs_winsize_", winsize, ".", Sys.Date(), ".bed")), sep = "\t", col.names = FALSE, row.names = FALSE)
  fwrite(genes.window, file = file.path(outdir, paste0("genes.", Sys.Date(), ".bed")), sep = "\t", col.names = FALSE, row.names = FALSE)
  
  
  # Write all? --------------------------------------------------------------
  
  jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
  all.window <- subset(dat.annot, type != "gene") %>%
    ungroup() %>%
    mutate(start = as.integer(start - winsize / 2),
           stop = as.integer(stop + winsize / 2),
           chr2 = paste("chr", chr, sep = ""),
           name2 = paste(name, strand, sep = "")) %>%
    dplyr::select(chr2, start, stop, name2) %>%
    filter(chr2 %in% jchromos)
  
  
  fwrite(all.window, file = file.path(outdir, paste0("all_winsize_", winsize, ".", Sys.Date(), ".bed")), sep = "\t", col.names = FALSE, row.names = FALSE)
  
  
  
}



