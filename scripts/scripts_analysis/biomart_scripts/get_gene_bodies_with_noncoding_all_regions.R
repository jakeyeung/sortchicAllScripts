# Jake Yeung
# Date of Creation: 2019-07-28
# File: ~/projects/scchic_gastru/scripts/scripts_analysis/scchicseq/get_gene_body/get_gene_bodies.R

rm(list=ls())

library(biomaRt)
library(dplyr)
library(GenomicRanges)

# Constants ---------------------------------------------------------------


# winsize <- 100000L
# use gene start and gene body 


upstream <- 2000L
downstream <- upstream

withchr <- TRUE

pfilt <- 0.8

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/genebodies/mm10_all_genes"
assertthat::assert_that(dir.exists(outdir))

if (withchr){
  outf <- file.path(outdir, paste0("gene_start_end.all_genes.", Sys.Date(), ".pfilt_", pfilt, ".up_", upstream, ".down_", downstream, ".bed"))
} else {
  outf <- file.path(outdir, paste0("gene_start_end.all_genes.", Sys.Date(), ".pfilt_", pfilt, ".up_", upstream, ".down_", downstream, ".nochr.bed"))
}
assertthat::assert_that(!file.exists(outf))


# Biomart init ------------------------------------------------------------


mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl')

gos <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", 
                          "chromosome_name", "transcription_start_site", 
                          "transcript_start", "transcript_end", "strand", "gene_biotype"),
             mart=mart.obj)

# 
print(unique(gos$chromosome_name))

chromos <- c(seq(19), c("X", "Y"))
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

nctypes <- "lncRNA
miRNA
misc_RNA
scRNA
snRNA
snoRNA
ribozyme
sRNA
scaRNA"
nctypes <- strsplit(nctypes, "\n")[[1]]

biotypes <- c("protein_coding", igvtypes, nctypes)
# biotypes <- c("protein_coding", igvtypes)
# biotypes <- igvtypes

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

# add window that is basically start and end position
dat.win <- gos.filt %>%
  dplyr::rename(gene = external_gene_name, seqnames = chromosome_name, tss = transcription_start_site, gstart = transcript_start, gend = transcript_end) %>%
  dplyr::select(gene, seqnames, tss, gstart, gend, isoform) %>%
  # define start and end, also extend upstream and downstream
  mutate(start = gstart - upstream, 
         end = gend + downstream,
         seqnames = paste("chr", seqnames, sep = "")) %>%
  dplyr::select(-tss, -gstart, -gend) %>%
  arrange(seqnames, start) %>%
  mutate(peak.name = paste(paste(seqnames, paste(start, end, sep = "-"), sep = ":"), gene, isoform, sep = ";")) %>%
  ungroup() %>%
  dplyr::select(seqnames, start, end, peak.name, gene) %>%
  rowwise() %>%
  mutate(gene.length = end - start)
  # mutate(dist = 0)

# take median size?
dat.win.filt <- dat.win %>%
  group_by(gene) %>%
  filter(gene.length == quantile(gene.length, p = pfilt, type = 3))


# # collapse identical start and end sites
# rows.keep <- which(!duplicated(subset(dat.win, select = c(seqnames, start, end))))
# dat.win.filt <- dat.win[rows.keep, ]

#

# Write -------------------------------------------------------------------

# write bedfile to table


if (withchr){
  data.table::fwrite(dat.win.filt, 
                     file = outf,
                     sep = "\t", col.names = FALSE, row.names = FALSE)
} else {
  data.table::fwrite(dat.win.filt %>% mutate(seqnames = gsub("chr", "", seqnames)),
                     file = outf,
                     sep = "\t", col.names = FALSE, row.names = FALSE)
}



