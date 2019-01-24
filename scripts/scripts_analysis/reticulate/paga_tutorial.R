# Jake Yeung
# Date of Creation: 2019-01-21
# File: ~/projects/scChiC/scripts/scripts_analysis/reticulate/paga_tutorial.R
# PAGA tutorial


library(reticulate)

repl_python()

import scanpy.api as sc

adata = sc.datasets.paul15()
sc.pp.recipe_zheng17(adata)

sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.draw_graph(adata)

sc.tl.louvain(adata, resolution=1.0)
sc.tl.paga(adata, groups='louvain')
sc.pl.paga(adata, color=['louvain'])

# try PAGA
adata_lda = sc.read("/Users/yeung/projects/scChiC/outputs_R/lda_output/chip.H3K27me3.normalized_counts.txt")

# root_cell = 50
# adata_lda.uns['iroot'] = root_cell

sc.tl.pca(adata_lda, svd_solver='arpack')
sc.pp.neighbors(adata_lda, n_neighbors=5, n_pcs=0)

sc.tl.draw_graph(adata_lda)
sc.tl.louvain(adata_lda, resolution=1.01)
sc.tl.paga(adata_lda, groups='louvain')
sc.pl.paga(adata_lda, color=['louvain'])


sc.tl.draw_graph(adata_lda, init_pos='paga')
sc.pl.draw_graph(adata_lda, color='louvain_anno', legend_loc='louvain')


# try with another dataset

adata = sc.read('/Users/yeung/Downloads/nestorowa_corrected_log2_transformed_counts.txt', cache=True)
# sc.pp.recipe_weinreb17(adata, log=False)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.louvain(adata, resolution=1)
sc.tl.draw_graph(adata, layout='fa', random_state=1)
# sc.pl.draw_graph(adata, color=['louvain'], legend_loc='on data')

sc.tl.paga(adata, groups='louvain')
sc.pl.paga(adata, color=['louvain'])
