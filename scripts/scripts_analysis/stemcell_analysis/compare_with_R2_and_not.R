# Jake Yeung
# Date of Creation: 2020-06-01
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/7-calculate_high_dimensional_fold_changes.R
# Summarize fold changes: find the different movements of genes in a genome-wide manner? 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(DropletUtils)

library(JFuncs)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

hubprefix <- "/home/jyeung/hub_oudenaarden"
jorg <- "org.Mm.eg.db"
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
jinf.tss <- file.path(hubprefix, "jyeung/data/databases/gene_tss/first_transcript_tss/gene_tss_winsize.10000.first_transcript.bed")


jmarks <- c("H3K4me1"); names(jmarks) <- jmarks

jmark <- "H3K4me3"

# Constants ---------------------------------------------------------------

    
# indir.bins <- file.path(hubprefix, "jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.count_tables_bl_filt")
# indir.bins <- file.path(hubprefix, "jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS.first_transcript")

indir.bins.orig <- file.path(hubprefix, "jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_SlidingWindow")
indir.bins.new <- file.path(hubprefix, "jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_SlidingWindow.slurm.noR2")
indir.bins.new.AddOne <- file.path(hubprefix, "jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_SlidingWindow.slurm.noR2.AddOne")
indir.bins.buys <- file.path(hubprefix, "jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_SlidingWindow.slurm.noR2.Buys")
indir.bins.buys.withR2 <- file.path(hubprefix, "jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_SlidingWindow.slurm.withR2.Buys")

inf.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.stringentDE.forAvO/filtered_cells_for_pseudobulk.2020-05-31.rds")

assertthat::assert_that(dir.exists(indir.bins.orig))
assertthat::assert_that(dir.exists(indir.bins.new))
assertthat::assert_that(dir.exists(indir.bins.buys))
assertthat::assert_that(dir.exists(indir.bins.buys.withR2))



# Load meta ---------------------------------------------------------------


dat.meta <- readRDS(inf.meta)


# Load bins matrix  -------------------------------------------------------

inf.old <- file.path(indir.bins.orig, paste0(jmark, ".countTableTSS.mapq_40.SlidingWindow_10000.blfiltered.csv"))
inf.new <- file.path(indir.bins.new, paste0(jmark, ".mapq_40.SlidingWindow_10000.blfiltered.csv"))
inf.new.AddOne <- file.path(indir.bins.new.AddOne, paste0(jmark, ".mapq_40.SlidingWindow_10000.blfiltered.csv"))
inf.buys <- file.path(indir.bins.buys, paste0(jmark, ".mapq_40.SlidingWindow_10000.blfiltered.csv"))
inf.buys.withR2 <- file.path(indir.bins.buys.withR2, paste0(jmark, ".mapq_40.SlidingWindow_10000.blfiltered.csv"))

mat.old <- ReadMatSlideWinFormat(inf.old, add.chromo = FALSE)
mat.new <- ReadMatSlideWinFormat(inf.new, add.chromo = FALSE)
mat.new.AddOne <- ReadMatSlideWinFormat(inf.new.AddOne, add.chromo = FALSE)

mat.buys <- ReadMatSlideWinFormat(inf.buys, add.chromo = FALSE)
mat.buys.withR2 <- ReadMatSlideWinFormat(inf.buys.withR2, add.chromo = FALSE)

mats.lst <- list(old = mat.old, new = mat.new, new.AddOne = mat.new.AddOne, buys.noR2 = mat.buys, buys.withR2 = mat.buys.withR2)

ncuts.lst <- lapply(mats.lst, colSums)

# Get cells ---------------------------------------------------------------

cells.keep <- unique(sapply(mats.lst, colnames))

bins.keep <- Reduce(f = intersect, x = lapply(mats.lst, rownames))

jcell <- cells.keep

plot(log2(ncuts.lst$new[cells.keep]), log2(ncuts.lst$buys.withR2[cells.keep]))

plot(ncuts.lst$new[cells.keep], ncuts.lst$new.AddOne[cells.keep])




