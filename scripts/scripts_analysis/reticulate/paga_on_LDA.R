# Jake Yeung
# Date of Creation: 2019-01-22
# File: ~/projects/scChiC/scripts/scripts_analysis/reticulate/paga_on_LDA.R
# PAGA on LDA


library(reticulate)
# use_python("/Users/yeung/anaconda3/bin/python")
use_condaenv("base", "/Users/yeung/anaconda3/bin/conda")

# need to set up ~/.Renviron RETICULATE_PYTHON=mypath https://stackoverflow.com/questions/50145643/unable-to-change-python-path-in-reticulate-r
reticulate::repl_python()  



import scanpy.api as sc

# try PAGA
adata_lda = sc.read("/Users/yeung/projects/scChiC/outputs_R/lda_output/chip.H3K27me3.normalized_counts.txt")

# root_cell = 50
# adata_lda.uns['iroot'] = root_cell

sc.tl.pca(adata_lda, svd_solver='arpack')
sc.pp.neighbors(adata_lda, n_neighbors=5, n_pcs=0)  # dont do PCA it's already reduced

sc.tl.draw_graph(adata_lda)
sc.tl.louvain(adata_lda, resolution=1.01)
sc.tl.paga(adata_lda, groups='louvain')
sc.pl.paga(adata_lda, threshold = 0.1)

# find marker genes
sc.tl.rank_genes_groups(adata_lda, 'louvain', method='t-test')
sc.pl.rank_genes_groups(adata_lda, n_genes=25, sharey=False)

# sc.tl.rank_genes_groups(adata_lda, 'louvain', method='logreg')
# sc.pl.rank_genes_groups(adata_lda, n_genes=25, sharey=False)

# sc.tl.draw_graph(adata_lda, init_pos='paga')
# sc.pl.draw_graph(adata_lda, color='louvain_anno', legend_loc='louvain')

