# Jake Yeung
# Date of Creation: 2021-02-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_deeper/3-run_GREAT_on_regions.compare_GC.compare_dists.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)

library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rGREAT)
library(JFuncs)



hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks



# Load bins  --------------------------------------------------------------

bsize <- 50000

inf.lost.lst <- lapply(jmarks, function(jmark){
  # inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/tables_top_6085_four_marks_dynamic_bins/lost_bins.", jmark, ".2021-02-18.bed"))
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/tables_top_6085_four_marks_dynamic_bins/lost_bins.", jmark, ".minlog2fc_0.5.2021-02-19.bed"))
  return(inf)
})

inf.high.lst <- lapply(jmarks, function(jmark){
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/High_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.", jmark, ".2021-01-30.txt"))
  
})

dat.lost.lst <- lapply(inf.lost.lst, function(inf){
  fread(inf)
})

dat.high.lst <- lapply(inf.high.lst, function(inf){
  print(inf)
  fread(inf)
})


# Run GREAT  --------------------------------------------------------------

regions.lost.lst <- lapply(dat.lost.lst, function(dat.lost){
  dat.lost %>%
    dplyr::rename(seqnames = V1,
                  start = V2,
                  end = V3,
                  region_coord = V4)
})

regions.high.lst <- lapply(jmarks, function(jmark){
  rbind(regions.lost.lst[[jmark]], dat.high.lst[[jmark]] %>% dplyr::select(seqnames, start, end, region_coord))
})


regions.range.lost.lst <- lapply(regions.lost.lst, function(regions.lost){
  makeGRangesFromDataFrame(as.data.frame(regions.lost)) 
})

regions.range.high.lst <- lapply(regions.high.lst, function(regions.high){
  makeGRangesFromDataFrame(as.data.frame(regions.high)) 
})


out.great.lost.k27 <- submitGreatJob(gr = regions.range.lost.lst$H3K27me3, bg = regions.range.high.lst$H3K27me3, species = "mm10", request_interval = 30)
out.great.lost.k9 <- submitGreatJob(gr = regions.range.lost.lst$H3K9me3, bg = regions.range.high.lst$H3K9me3, species = "mm10", request_interval = 30)

out.tb.lost.k27 <- getEnrichmentTables(out.great.lost.k27, ontology=availableOntologies(out.great.lost.k27))
out.tb.lost.k9 <- getEnrichmentTables(out.great.lost.k9, ontology=availableOntologies(out.great.lost.k9))


# Check lost bins  --------------------------------------------------------

out.tb.lost.k27.ordered <- lapply(out.tb.lost.k27, function(x){
  x %>% arrange(Hyper_Adjp_BH)
})

out.tb.lost.k9.ordered <- lapply(out.tb.lost.k9, function(x){
  x %>% arrange(Hyper_Adjp_BH)
})

print(head(out.tb.lost.k27.ordered$`GO Biological Process`))
print(head(out.tb.lost.k27.ordered$`GO Molecular Function`))
print(head(out.tb.lost.k27.ordered$`GO Cellular Component`))

print(head(out.tb.lost.k9.ordered$`GO Biological Process`))
print(head(out.tb.lost.k9.ordered$`GO Molecular Function`))
print(head(out.tb.lost.k9.ordered$`GO Cellular Component`))

# Write tables ------------------------------------------------------------





