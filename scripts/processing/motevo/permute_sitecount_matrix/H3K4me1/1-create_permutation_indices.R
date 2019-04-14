# Jake Yeung
# 1-create_permutation_indices.R
# 2019-04-01
# DESCRIPTION
# 
#     Create permutation index which will be read later row by row
# 
# FOR HELP
# 
#     Rscript 1-create_permutation_indices.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-04-01
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

set.seed(0)
inf <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_H3K4me1/mara_input/sitecounts_mat.rds"
outf <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_H3K4me1/mara_input/row_indices_shuffled.rds"
N <- readRDS(inf)  # get number of rows

n.shuffles <- 10000

# permute from row.seq

mat.out <- matrix(nrow = n.shuffles, ncol = nrow(N))
print(mat.out[1:5, 1:5])
for (i in seq(n.shuffles)){
  row.i.shuffled <- sample.int(nrow(N), nrow(N), replace=FALSE) # shuffle it up
  mat.out[i, ] <- row.i.shuffled
}
print(mat.out[1:5, 1:5])
# write output
saveRDS(mat.out, file = outf)
