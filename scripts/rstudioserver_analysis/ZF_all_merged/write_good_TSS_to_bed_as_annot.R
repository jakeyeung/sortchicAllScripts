# Jake Yeung
# Date of Creation: 2020-06-18
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/write_good_TSS_to_bed_as_annot.R
# These TSS are H3K4me3 reference with the highest signal 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(JFuncs)

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/zebrafish.SelectedTSSfromWKM"
dir.create(outdir)

out.annot <- file.path(outdir, "gene_tss.SelectedFromWKM.winsize_50000.species_drerio.bed")

# oad bins ----------------------------------------------------------------

jdate <- "2020-06-09"
indir.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/setup_matrix_for_poisson_regression.likeBM.redo_count_tables"
assertthat::assert_that(dir.exists(indir.tss))

jprefix <- file.path(indir.tss, paste0("integrated_analysis.", jdate, ".UseTSSfromH3K4me3.likeBM"))

infrdata <- paste0(jprefix, ".RData")
assertthat::assert_that(file.exists(infrdata))

load(infrdata, v=T)


# Write rownames to output  -----------------------------------------------

rnames.all <- lapply(tss.mats.sc.filt.zf, function(jmat){
  rownames(jmat)
})

rnames.common <- Reduce(intersect, rnames.all)

tss.coord <- sapply(rnames.common, function(x) strsplit(x, ";")[[1]][[1]])


# Write to output ---------------------------------------------------------

dat.bed <- data.frame(chromo = sapply(tss.coord, GetChromo), 
                      Start = sapply(tss.coord, GetStart, returnAsInt=TRUE),
                      End = sapply(tss.coord, GetEnd, returnAsInt=TRUE),
                      coord = rnames.common, 
                      stringsAsFactors = FALSE)

fwrite(dat.bed, file = out.annot, quote = FALSE, sep = "\t", col.names = FALSE)




