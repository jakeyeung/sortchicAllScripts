# Jake Yeung
# Date of Creation: 2018-12-22
# File: ~/projects/scChiC/scripts/scripts_analysis/metacell/run_metacell_on_merged_peaks.R
# Run metacell on merged peaks 

library(metacell)

indat <- "/tmp/metacell_inputs/PZ-BM-H3K4me1.merged.25000kbmerge.mat"
outdir <- "/tmp/metacell_outputs.merged"
# varthres <- 0.01
dir.create(outdir)

# varthres <- 0.001
# if (is.na(varthres)){
#   stop(paste("Threshold must be numeric, found:"), args[[3]])
# }

scdb_init(outdir, force_reinit=T)

matname <- "H3K4me1_merged"
matname.filt.cells <- "H3K4me1_filt"
matname.filt.bins <- "H3K4me1_noemptbin"
jgset <- "H3K4me1_feats"
jgraph <- "H3K4me1_graph"
jcoc <- "H3K4me1_coc500"
jmc <- "H3K4me1_mc"
newjmc <- "H3K4me1_nemc"
colorizepath <- "/private/tmp/metacell_inputs/pbmc_mc_colorize.txt"
jproj <- "H3K4me1_2dproj"
jdset <- "H3K4me1"

# minUMIs <- 6000
# minUMIs <- 0  # we take subset of rows so this is meaningless here
# maxUMIs <- 10000000
# maxbinsums <- 10000

figsdir <- file.path(outdir, "figs")
scfigs_init(figsdir)

mcell_import_scmat_tsv(matname, fn=indat, dset_nm=jdset)
mat = scdb_mat(matname)
print(dim(mat@mat))

# remove bad peaks
binsums <- Matrix::rowSums(as.matrix(mat@mat))
# emptybins <- names(which(binsums <= 1 | binsums > maxbinsums))
emptybins <- vector('integer')  # no bad peaks
mcell_mat_ignore_genes(new_mat_id=matname.filt.bins, mat_id=matname, ig_genes=emptybins)

mat <- scdb_mat(matname.filt.bins)

# ignore bad cells
# cell_sizes <- Matrix::colSums(mat@mat)
# large_cells <- names(which(cell_sizes>maxUMIs))
# small_cells <- names(which(cell_sizes<minUMIs))
bad_cells <- c(small_cells, large_cells)
bad_cells <- vector('integer')

print("Plot UMIs after")
mcell_plot_umis_per_cell(matname.filt.bins)

mcell_mat_ignore_cells(new_mat_id=matname.filt.cells, mat_id=matname.filt.bins, ig_cells = bad_cells)
mat <- scdb_mat(matname.filt.cells)

print("Dimensions after filtering")
print(dim(mat@mat))


print("Plot UMIs after")
mcell_plot_umis_per_cell(matname.filt.cells)

print("Add gene stat")
mcell_add_gene_stat(gstat_id=matname.filt.cells, mat_id=matname.filt.cells, force=T)


print("Filter varmean")
mcell_gset_filter_varmean(gset_id=jgset, gstat_id=matname.filt.cells, T_vm=0, force_new=T)
# print("Filter cov")
# mcell_gset_filter_cov(gset_id = jgset, gstat_id=matname.filt.cells, T_tot=99, T_top3=Inf)
# mcell_plot_gstats(gstat_id=matname.filt.cells, gset_id=jgset)

print("Plot gstats")
mcell_plot_gstats(gstat_id=matname.filt.cells, gset_id=jgset)

print("Add cgraph")

gset1 = scdb_gset(jgset)
feat1 = gset_get_feat_mat(jgset, matname.filt.cells, downsamp=F, add_non_dsamp=F)
mcell_add_cgraph_from_mat_bknn(mat_id=matname.filt.cells, 
                               gset_id = jgset, 
                               graph_id=jgraph,
                               K=25,  # tuning param
                               dsamp=F)


print("Coclust")
mcell_coclust_from_graph_resamp(
  coc_id=jcoc, 
  graph_id=jgraph,
  min_mc_size=10, 
  p_resamp=0.75, n_resamp=500)

print("MC from coclust")
mcell_mc_from_coclust_balanced(
  coc_id=jcoc, 
  mat_id=matname.filt.cells,
  mc_id=jmc, 
  K=15, min_mc_size=15, alpha=2)

print("Plot outlier")
mcell_plot_outlier_heatmap(mc_id=jmc, mat_id = matname.filt.cells, T_lfc=3)

print("MC split filt")
mcell_mc_split_filt(new_mc_id=newjmc, 
                    mc_id=jmc, 
                    mat_id=matname.filt.cells,
                    T_lfc=3, plot_mats=F)

print("Gset from markers")
mcell_gset_from_mc_markers(gset_id=jgset, mc_id=jmc)

marks_colors = read.table(colorizepath, sep="\t", h=T, stringsAsFactors=F)
mc_colorize(newjmc, marker_colors=marks_colors)
mc = scdb_mc(newjmc)
mcell_mc_plot_marks(mc_id=newjmc, gset_id=jgset, mat_id=matname.filt.cells)

lfp = log2(mc@mc_fp)

mcell_mc2d_force_knn(mc2d_id=jproj,mc_id=,jmc, graph_id=jgraph)
tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id=jproj)

mc_hc = mcell_mc_hclust_confu(mc_id=jmc, graph_id=jgraph)

mc_sup = mcell_mc_hierarchy(mc_id=jmc, mc_hc=mc_hc, T_gap=0.04)
mcell_mc_plot_hierarchy(mc_id=jmc, 
                        graph_id=jgraph, 
                        mc_order=mc_hc$order, 
                        sup_mc = mc_sup, 
                        width=2800, heigh=2000, min_nmc=2)

# mcell_mc2d_force_knn(mc2d_id=jproj,mc_id=jmc, graph_id=jgraph)
# mcell_mc2d_plot(mc2d_id=jproj)




# Are cells in cluster 5 different from others? ---------------------------

# load("/private/tmp/metacell_outputs.merged/mc2d.H3K4me1_2dproj.Rda", v=T)

# jcol <- mc@mc

jmc2d = scdb_mc2d(jproj)


clstr.assigns <- mc@mc[names(jmc2d@sc_x)]
clstr5.vs.all <- sapply(clstr.assigns, function(x) ifelse(x == 5, "blue", "red"))
# clst5 <- 
plot(jmc2d@sc_x, jmc2d@sc_y, pch = 20, col = clstr5.vs.all)

# get average clstr 5 versus all and find differentially expressed peaks
# mat = scdb_mat(matname.filt.cells)
clstr.fg <- 5
cells.fg <- which(mc@mc == clstr.fg)
print(length(cells.fg))

cells.bg <- which(mc@mc != clstr.fg)
print(length(cells.bg))

# top hits are chrY, chr6, chr7
gset.obj@gene_set


# Write metacell assignments to output ------------------------------------

outdir <- "outputs_R/metacell_output"
dir.create(outdir)
outf <- file.path(outdir, "PZ-BM-H3K4me1.merged.25000kbmerge.cluster_assign.txt")
outdf <- data.frame(bamname = names(mc@mc),
                    clusterid = mc@mc)
write.table(outdf, file = outf, quote = FALSE, sep = "\t", row.names = FALSE)



