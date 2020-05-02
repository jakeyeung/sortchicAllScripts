# Jake Yeung
# Date of Creation: 2020-04-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/get_gene_TSS_zebrafish_again_from_server.R
# /home/jyeung/projects/scchic/scripts/scripts_analysis/biomart_scripts/get_gene_TSS_zebrafish.R 


rm(list=ls())

library(biomaRt)
library(dplyr)
library(GenomicRanges)


# Constants ---------------------------------------------------------------

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/zebrafish"
winsizes <- c(1000L, 5000L, 10000L)

for (winsize in winsizes){
  print(winsize)
  # Biomart init ------------------------------------------------------------
  
  species <- "drerio"
  mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = paste0(species, '_gene_ensembl'))
  
  gos <- getBM(
    attributes=c("external_gene_name", "chromosome_name", "transcription_start_site", "start_position", "end_position", "strand", "gene_biotype"),
    mart=mart.obj)
  
  # 
  print(unique(gos$chromosome_name))
  
  chromos <- c(seq(25))
  chromos.withprefix <- paste("chr", chromos, sep = "")
  
  # Processs ----------------------------------------------------------------
  
  # get biotypes
  igvtypes <- "IG_C_gene
IG_D_gene
IG_J_gene
IG_LV_gene
IG_V_gene
TR_C_gene
TR_J_gene
TR_V_gene
TR_D_gene"
  igvtypes <- strsplit(igvtypes, "\n")[[1]]
  
  biotypes <- c("protein_coding", igvtypes)
  
  gos.filt <- gos %>% filter(chromosome_name %in% chromos) %>%
    filter(gene_biotype %in% biotypes) %>%
    filter(!grepl("^Gm", external_gene_name)) %>%
    group_by(chromosome_name) %>%
    arrange(transcription_start_site) %>%
    group_by(external_gene_name) %>%
    mutate(isoform = seq(length(external_gene_name)))
  
  # check how far TSSs are apart from each other
  
  gos.sum <- gos.filt %>%
    group_by(external_gene_name) %>%
    summarise(jmean = mean(transcription_start_site),
              jmad = mad(transcription_start_site),
              ntx = length(transcription_start_site)) %>%
    arrange(desc(jmad))
  
  # add isoform number based on TSS
  
  # add window
  
  dat.win <- gos.filt %>%
    dplyr::rename(gene = external_gene_name, seqnames = chromosome_name, tss = transcription_start_site, gstart = start_position, gend = end_position) %>%
    dplyr::select(gene, seqnames, tss, isoform) %>%
    mutate(start = as.integer(tss - winsize / 2), 
           end = as.integer(tss + winsize / 2),
           seqnames = paste("chr", seqnames, sep = "")) %>%
    dplyr::select(-tss) %>%
    arrange(seqnames, start) %>%
    mutate(peak.name = paste(paste(seqnames, paste(start, end, sep = "-"), sep = ":"), gene, isoform, sep = ";")) %>%
    ungroup() %>%
    dplyr::select(seqnames, start, end, peak.name)
  # mutate(dist = 0)
  
  #
  # set negatives to zero
  dat.win <- dat.win %>%
    rowwise() %>%
    mutate(start = ifelse(start < 0, 0, start))
  
  dat.win.nochr <- dat.win %>%
    ungroup() %>%
    mutate(seqnames = gsub("^chr", "", seqnames))
  
  # make sure it's all integers
  
  # Write -------------------------------------------------------------------
  
  # write bedfile to table
  outf <- file.path(outdir, paste0("gene_tss.winsize_", winsize, ".species_", species, ".bed"))
  outf.nochr <- file.path(outdir, paste0("gene_tss.winsize_", winsize, ".species_", species, ".nochr.bed"))
  data.table::fwrite(dat.win, file = outf, sep = "\t", col.names = FALSE, row.names = FALSE)
  data.table::fwrite(dat.win.nochr, file = outf.nochr, sep = "\t", col.names = FALSE, row.names = FALSE)
}

