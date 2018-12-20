# Jake Yeung
# Date of Creation: 2018-12-20
# File: ~/projects/scChiC/scripts/scripts_analysis/metacell/metacell_tutorial.R
# Run tutorial from 
# https://tanaylab.bitbucket.io/metacell-r/articles/basic_pbmc8k.html


if(!dir.exists("testdb")) dir.create("testdb/")
scdb_init("testdb/", force_reinit=T)

mcell_import_scmat_10x("test", base_dir="http://www.wisdom.weizmann.ac.il/~atanay/metac_data/pbmc_8k/")
mat = scdb_mat("test")
print(dim(mat@mat))

if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")


mcell_plot_umis_per_cell("test")

mat = scdb_mat("test")
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
ig_genes = c(grep("^IGJ", nms, v=T), 
             grep("^IGH",nms,v=T),
             grep("^IGK", nms, v=T), 
             grep("^IGL", nms, v=T))

bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))


mcell_mat_ignore_genes(new_mat_id="test", mat_id="test", bad_genes, reverse=F) 


mcell_mat_ignore_small_cells("test", "test", 800)


mcell_add_gene_stat(gstat_id="test", mat_id="test", force=T)


mcell_gset_filter_varmean(gset_id="test_feats", gstat_id="test", T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id = "test_feats", gstat_id="test", T_tot=100, T_top3=2)


mcell_plot_gstats(gstat_id="test", gset_id="test_feats")


gset = scdb_gset("test_feats")
feat2 = gset_get_feat_mat(gset_id = "test_feats", mat_id = "test", downsamp=T, add_non_dsamp=T)
mcell_add_cgraph_from_mat_bknn(mat_id="test", 
                               gset_id = "test_feats", 
                               graph_id="test_graph",
                               K=100,
                               dsamp=T)