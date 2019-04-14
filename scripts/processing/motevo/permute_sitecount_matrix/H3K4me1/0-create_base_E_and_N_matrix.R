# Jake Yeung
# 0-create_base_E_and_N_matrix.R
# 2019-04-01
# DESCRIPTION
# 
#     Create E and N matrix before permuting the rows in N
# 
# FOR HELP
# 
#     Rscript 0-create_base_E_and_N_matrix.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-04-01
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)


inf <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/H3K4me1/mara_output/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me1.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K50/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me1.filt_0.99.center_TRUE_K50/r.Robj"

# write Emat and Nmat here
outdir <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_H3K4me1/mara_input"

load(inf, v=T)

saveRDS(E, file.path(outdir, "exprs_mat.rds"))
saveRDS(N, file.path(outdir, "sitecounts_mat.rds"))

