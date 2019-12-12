# Jake Yeung
# Date of Creation: 2019-12-09
# File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/retag_analysis/LDA_analysis.R
# Analysis 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

# Constants ---------------------------------------------------------------

inmain <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZFWKM_allmarks_mergedtagged_dedupfixed_redo"


# Load umaps  -------------------------------------------------------------



