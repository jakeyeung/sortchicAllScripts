# Jake Yeung
# Date of Creation: 2021-09-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/5-check_Ku_2021_new_tables.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Load bed  ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_total_reads_with_remapped_Ku2021.", Sys.Date(), ".NewTables.pdf")
# pdf(outpdf, useDingbats = FALSE)

inmain <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Ku_et_al_2021")

inf.k4 <- file.path(inmain, "GSE139857_scH3K4me3_7798_cell_counts_at_WBCpeaks.txt")
inf.k4.rnames <- file.path(inmain, "GSE139857_peakname3.txt")
inf.k4.cnames <- file.path(inmain, "from_email_KejiZhao_WaiLim/wc_uniq_H3K4me3.txt")
inf.k27 <- file.path(inmain, "GSE139857_scK27_4lanes_rc_at_79100peaks_width10k.txt")
inf.k27.rnames <- file.path(inmain, "from_email_KejiZhao_WaiLim/combined_sicer_E100_peaks.sort.merge.txt")
inf.k27.cnames <- file.path(inmain, "from_email_KejiZhao_WaiLim/wc_uniq_H3K27me3.txt")

mat.k4 <- fread(inf.k4)
bed.k4 <- fread(inf.k4.rnames)
rnames.k4 <- paste(paste(bed.k4$V1, bed.k4$V2, sep = ":"), bed.k4$V3, sep = "-")

mat.k27 <- fread(inf.k27)
bed.k27 <- fread(inf.k27.rnames)
rnames.k27 <- paste(paste(bed.k27$V1, bed.k27$V2, sep = ":"), bed.k27$V3, sep = "-")

dat.total.k27 <- fread(inf.k27.cnames)
dat.total.k4 <- fread(inf.k4.cnames)

cnames.k27 <- dat.total.k27$V2


mat.k4.sparse <- Matrix(as.matrix(mat.k4), sparse = TRUE)
mat.k27.sparse <- Matrix(as.matrix(mat.k27), sparse = TRUE)

rm(mat.k4)
rm(mat.k27)

rownames(mat.k4.sparse) <- rnames.k4
rownames(mat.k27.sparse) <- rnames.k27
colnames(mat.k27.sparse) <- cnames.k27


# Count unique reads  -----------------------------------------------------

jmarks <- c("H3K4me3", "H3K27me3")
mat.lst <- list(mat.k4.sparse, mat.k27.sparse)
names(mat.lst) <- jmarks

dat.total.lst <- list(dat.total.k4, dat.total.k27)
names(dat.total.lst) <- jmarks

dat.peak.sums.lst <- lapply(mat.lst, function(jmat){
  data.frame(cell = colnames(jmat), cuts_in_peaks = colSums(jmat), stringsAsFactors = FALSE)
})

dat.total.sums.lst <- lapply(dat.total.lst, function(dat.total){
  data.frame(cell = dat.total$V2, cuts_total = dat.total$V1, stringsAsFactors = FALSE)
})

dat.merge.sums.long <- lapply(jmarks, function(jmark){
  left_join(dat.peak.sums.lst[[jmark]], dat.total.sums.lst[[jmark]]) %>%
    mutate(mark = jmark)
})  %>%
  bind_rows() %>%
  mutate(dataset = "Ku2021_scChICseq")

ggplot(dat.merge.sums.long, aes(x = mark, y = cuts_in_peaks / cuts_total)) + 
  geom_boxplot() + 
  geom_point() +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# save to file 
outf <- file.path(hubprefix, paste0("jyeung/data/scChiC/public_data/Ku_et_al_2021/processed_tables/total_cuts_peak_cuts_Ku_2021.", Sys.Date(), ".txt"))
fwrite(dat.merge.sums.long, file = outf, sep = "\t")






