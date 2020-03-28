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

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/ENCODE"
winsize <- 10000L  # 5kb left and right

inf.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/ENCODE/genomic_regions_10kbpromoters.txt"
dat.annot <- fread(inf.annot) %>%
  mutate(start = as.integer(start),
         stop = as.integer(stop))


proms <- subset(dat.annot, type == "promotor") 

enhs <- subset(dat.annot, type == "enhancer")


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


# bed file should be chromo, start, stop, geneid
fwrite(proms.pos, file = file.path(outdir, paste0("promotersPos_winsize_0.bed")), sep = "\t", col.names = FALSE, row.names = FALSE)
fwrite(proms.neg, file = file.path(outdir, paste0("promotersNeg_winsize_0.bed")), sep = "\t", col.names = FALSE, row.names = FALSE)
fwrite(enhs, file = file.path(outdir, paste0("enhancers_winsize_0.bed")), sep = "\t", col.names = FALSE, row.names = FALSE)


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

fwrite(proms.window, file = file.path(outdir, paste0("promoters_winsize_", winsize, ".bed")), sep = "\t", col.names = FALSE, row.names = FALSE)
fwrite(enhs.window, file = file.path(outdir, paste0("enhancers_winsize_", winsize, ".bed")), sep = "\t", col.names = FALSE, row.names = FALSE)
fwrite(promsenhs.window, file = file.path(outdir, paste0("promsenhs_winsize_", winsize, ".bed")), sep = "\t", col.names = FALSE, row.names = FALSE)


  # Write all? --------------------------------------------------------------

all.window <- subset(dat.annot, type != "gene") %>%
  ungroup() %>%
    mutate(start = as.integer(start - winsize / 2),
           stop = as.integer(stop + winsize / 2),
           chr2 = paste("chr", chr, sep = ""),
           name2 = paste(name, strand, sep = "")) %>%
    dplyr::select(chr2, start, stop, name2)


fwrite(all.window, file = file.path(outdir, paste0("all_winsize_", winsize, ".bed")), sep = "\t", col.names = FALSE, row.names = FALSE)



