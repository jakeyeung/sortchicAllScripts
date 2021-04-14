# Jake Yeung
# 2-filter_mat_by_cells.R
# 2020-08-23
# DESCRIPTION
# 
#     Filters csv file by cells from an .rds mat
# 
# FOR HELP
# 
#     Rscript 2-filter_mat_by_cells.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2020-08-23
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

library(data.table)
library(scchicFuncs)
library(dplyr)

inf.csv <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/mouse/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters/mouse-BM-EtOH-H3K4me3-1.sorted.tagged.countTable.TSS.csv"
inf.rds <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda/H3K4me3_BM.rds"
outf.rds <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/mouse/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters/mouse-BM-EtOH-H3K4me3-1.sorted.tagged.countTable.TSS.good_cells_filt.rds"

dat.new <- ReadMatTSSFormat(inf = inf.csv, as.sparse = TRUE, add.coord = FALSE, sort.rnames = FALSE)

dat.old <- readRDS(inf.rds)

cnames.keep <- colnames(dat.old)

cnames.keep.i <- colnames(dat.new) %in% cnames.keep

dim(dat.new)
dat.new.filt <- dat.new[, cnames.keep.i]
dim(dat.new.filt)
saveRDS(dat.new.filt, file = outf.rds)
