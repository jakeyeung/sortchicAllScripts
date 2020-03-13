# Jake Yeung
# Date of Creation: 2020-03-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/2-explore_Giladi_et_al_write_gene_tables_from_celltypefilt_DE_toponly.R
# 

rm(list=ls())

library(topicmodels)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(scchicFuncs)
library(ggforce)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)

library(ggforce)

library(biomaRt)
library(dplyr)
library(GenomicRanges)

# States ------------------------------------------------------------------

keeptop <- 200
logfc.min <- 0
pvalmax <- 0.01
winsize <- 2
jsign <- "low"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/giladi_filtered.seurat"
dir.create(outdir)
inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat.celltypes_filt.rds"
jbase <- basename(inf)
jbase <- ClipLast(jbase, jsep = "\\.", jsep.out = ".")

outf <- file.path(outdir, paste0("gene_tss_winsize.", winsize, ".", jbase, ".pval_", pvalmax, ".logfc_", logfc.min, ".", jsign, ".bed"))
de.output <- readRDS(inf)

head(de.output)

if (jsign == "high"){
  jsub <- de.output %>%
    group_by(cluster) %>%
    filter(p_val_adj <= pvalmax & avg_logFC > logfc.min) %>%
    arrange(desc(abs(avg_logFC)))
} else if (jsign == "low") {
  jsub <- de.output %>%
    group_by(cluster) %>%
    filter(p_val_adj <= pvalmax & avg_logFC < logfc.min) %>%
    arrange(desc(abs(avg_logFC)))
} else {
  warning("jsign must be high or low: found", jsign)
}

jsub.count <- jsub %>%
  group_by(cluster) %>%
  summarise(ngene = length(gene))
print(jsub.count)

genes.keep <- unique(jsub$gene)


# Biomart init ------------------------------------------------------------

mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl')

gos <- getBM(
  attributes=c("external_gene_name", "chromosome_name", "transcription_start_site", "start_position", "end_position", "strand", "gene_biotype"),
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

biotypes <- c("protein_coding", igvtypes)

gos.filt <- gos %>% 
  filter(chromosome_name %in% chromos) %>%
  filter(gene_biotype %in% biotypes) %>%
  filter(!grepl("^Gm", external_gene_name)) %>%
  filter(external_gene_name %in% genes.keep) %>%
  group_by(chromosome_name) %>%
  arrange(transcription_start_site) %>%
  group_by(external_gene_name) %>%
  mutate(isoform = seq(length(external_gene_name)))

# check how far TSSs are apart from each other

# add isoform number based on TSS

# add window

dat.win <- gos.filt %>%
  dplyr::rename(gene = external_gene_name, seqnames = chromosome_name, tss = transcription_start_site, gstart = start_position, gend = end_position) %>%
  dplyr::select(gene, seqnames, tss, isoform) %>%
  mutate(start = tss - winsize / 2,
         end = tss + winsize / 2,
         seqnames = paste("chr", seqnames, sep = "")) %>%
  dplyr::select(-tss) %>%
  arrange(seqnames, start) %>%
  mutate(peak.name = paste(paste(seqnames, paste(start, end, sep = "-"), sep = ":"), gene, isoform, sep = ";")) %>%
  ungroup() %>%
  dplyr::select(seqnames, start, end, peak.name)

# Write -------------------------------------------------------------------

# write bedfile to table
data.table::fwrite(dat.win, file = outf, sep = "\t", col.names = FALSE, row.names = FALSE)



