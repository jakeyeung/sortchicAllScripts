# Jake Yeung
# Date of Creation: 2022-04-27
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/20-plot_dynamic_bins_across_pseudotimes_four_marks_cleaned_batch_corrected_fewer_ctypes_TSS_marker_genes_genesets.R
# =

rm(list=ls())

library(scchicFuncs)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)

jratio <- 0.66

jmarks <- c("k4me1", "k4me3", "k27me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarksold) <- jmarks


# Load new colors ---------------------------------------------------------

inf.colors.fixed <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/ctypes_on_umap_batch_corrected_colors_fixed/dat_colors_DC_monocyte_fixed.2022-05-17.txt"
dat.colors.fixed <- fread(inf.colors.fixed)


# Load meta ----------------------------------------------------------------

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt")
  dat.meta <- fread(inf.meta) %>%
    left_join(., dat.colors.fixed) %>%
    rowwise() %>%
    mutate(colcode = colcodenew)
  # replace colcode with colcodenew
})

dat.meta.colors <- subset(dat.meta.lst$k4me1, select = c(ctype.from.LL, colcode))
ctype2col <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)





# Load impute -------------------------------------------------------------


dat.impute.lst <- lapply(jmarks, function(jmark){
  inf.impute <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/primetime_objects/dat_impute_tss_", jmark, ".2022-04-24.rds")
  readRDS(inf.impute)
})

# Load trajs --------------------------------------------------------------

indir.traj <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned"
dat.trajs <- lapply(jmarks, function(jmark){
  inf.trajs <- file.path(indir.traj, paste0("trajs_outputs.", jmark, ".rds"))
  if (jmark == "k27me3"){
    inf.trajs <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned2_batch_corrected_eryth_fix/trajs_outputs_batch_corrected.k27me3.2022-04-21.rds")
  }
  print(inf.trajs)
  readRDS(inf.trajs)
})



# Get marker genes --------------------------------------------------------

dat.markers <- fread("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/genesets/geneset_on_umap.binskeep_0.niter_500.2020-12-09.from_LDA_topics.condensed.heatmap.famousgenes.keepn_400.refmark_H3K4me3.2020-12-09.txt")

print(unique(dat.markers$jset))

jsets <- c("Bcells", "Eryths", "Granulocytes"); names(jsets) <- jsets
# jsets <- c("HSPCs"); names(jsets) <- jsets

# markergenes <- c("Ly6c2", "Ly6g", "S100a8", "S100a2", "Chil3", "Sox6", "Tal1", "Gata1", "Ebf1", "Cd79a", "Cd79b", "Hoxa9", "Meis1", "Runx2", "Kit", "Hlf", "Erg", "Cd34", "Stat4", "Tcf7", "Cebpe", "Pax5", "Cd180", "Cd72", "Blnk")
# markergenes.grep <- paste(markergenes, collapse = "$|")

# rnames <- rownames(dat.impute.lst$k4me1)
# rnames.markers <- rnames[grepl(markergenes.grep, rnames)]



# check 
jmarktest <- "k9me3"
jmarktest <- "k27me3"
jctypetest <- "Granulocytes"


outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umap_and_trajs/k27me3_batch_corrected_fewer_ctypes_TSS_genesets"
dir.create(outdir)

# jctypes <- c("Granulocytes", "Bcells", "Eryths")
jctypes <- c("Granulocytes")
names(jctypes) <- jctypes

# jrname <- grep("Cd72", rows.keep, value = TRUE)
# jrname <- grep("Ly6c2", rows.keep, value = TRUE)

parallel::mclapply(jmarks, function(jmark){
  
  print(jmark)
  
  # cell2ctype <- 
  
  for (jctype in jctypes){
    
    # pdfout <- file.path(outdir, paste0("umap_on_trajs_fix_color_TSS_common_root_fewer_plots_no_alpha.", jmark, ".", jctype, ".", Sys.Date(), ".pdf"))
    pdfout <- file.path(outdir, paste0("umap_on_trajs_fix_color_TSS_common_root_fewer_plots_colors_fixed_HSPCs.", jmark, ".", jctype, ".", Sys.Date(), ".pdf"))
    pdf(pdfout, useDingbats = FALSE)
    print(jctype)
    
    # rnames.markers <- unique(subset(dat.markers, jset == jctype)$gene)
    rnames.markers <- unique(subset(dat.markers, jset == "HSPCs")$gene)
    
    rows.keep.i <- which(rownames(dat.impute.lst[[jmark]]) %in% rnames.markers)
    print(length(rows.keep.i))
    
    rows.keep <- rownames(dat.impute.lst[[jmark]])[rows.keep.i]
    
    print(length(rows.keep))
    
      jtitle <- paste(jmark, jctype, "nbins:", length(rows.keep))
      
      dat.signal <- data.frame(cell = colnames(dat.impute.lst[[jmark]]), 
                               signal = colMeans(dat.impute.lst[[jmark]][rows.keep, ]), 
                               stringsAsFactors = FALSE)
      dat.exprs <- dat.signal %>%
        left_join(dat.meta.lst[[jmark]], ., by = "cell")
      
      m.signal <- ggplot(dat.exprs, aes(x = umap1, y = umap2, color = signal)) + 
        geom_point() + 
        scale_color_viridis_c() + 
        ggtitle(jtitle) + 
        theme_bw() + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      m.signal.batch <- ggplot(dat.exprs, aes(x = umap1, y = umap2, color = signal)) + 
        geom_point() + 
        scale_color_viridis_c() + 
        facet_wrap(~batch) + 
        ggtitle(jtitle) + 
        theme_bw() + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      print(m.signal)
      print(m.signal.batch)
      
      # show trajectory
      
      traj.ctypes <- names(dat.trajs[[jmark]])
      names(traj.ctypes) <- traj.ctypes
      dat.trajs.long <- lapply(traj.ctypes, function(traj.ctype){
        dat.trajs.sub <- dat.trajs[[jmark]][[traj.ctype]] %>%
          filter(is.ctype) %>%
          dplyr::select(cell, ctype.ordered, ptime) %>%
          mutate(traj = traj.ctype)
        # colcode = ctype2col[[traj.ctype]])
      }) %>%
        bind_rows() %>%
        left_join(., dat.signal, by = "cell") %>%
        left_join(., subset(dat.meta.lst[[jmark]], select = c(cell, batch, ctype.from.LL, colcode)), by = "cell")
      
      # order so Granulocyte is first? 
      dat.trajs.long <- dat.trajs.long %>%
        arrange(traj == jctype) %>%
        rowwise() %>%
        # mutate(jalpha = ifelse(traj == jctype, 0.9, 0.5),
        mutate(jalpha = ifelse(traj == jctype, 0.9, 0.9),
               colcodetraj = ctype2col[[traj]])
      
      m.traj <- ggplot(dat.trajs.long, aes(x = ptime, y = signal, color = colcode, alpha = jalpha)) + 
        geom_point() + 
        scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                              guide = "legend") + 
        theme_bw() + 
        ggtitle(paste(jmark, jctype)) + 
        theme(aspect.ratio=jratio, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
      print(m.traj)
      
      # fit with common root
      dat.hscs <- subset(dat.trajs.long, ctype.from.LL %in% c("HSCs", "LT", "ST", "MPPs"))
      dat.nonhscs <- subset(dat.trajs.long, !ctype.from.LL %in% c("HSCs", "LT", "ST", "MPPs"))
      dat.commonroot <- lapply(split(dat.nonhscs, f = dat.nonhscs$traj), function(jdat){
        jtraj <- unique(jdat$traj)
        jcolcode <- unique(jdat$colcode)
        jcolcodetrj <- unique(jdat$colcodetraj)
        jdat <- rbind(jdat, dat.hscs)
        jdat$traj <- jtraj
        jdat$colcodetraj <- jdat$colcodetraj
        jdat$colcode <- jdat$colcode
        return(jdat)
      }) %>%
        bind_rows()
      
      ptime.cutoff <- 0.1
      
      dat.hscs <- subset(dat.trajs.long, ptime < ptime.cutoff)
      dat.nonhscs <- subset(dat.trajs.long, ptime > ptime.cutoff)
      
      dat.commonroot <- lapply(split(dat.nonhscs, f = dat.nonhscs$traj), function(jdat){
        jtraj <- unique(jdat$traj)
        jcolcodetraj <- unique(jdat$colcodetraj)
        
        # print(jtraj)
        # print(jcolcodetraj)
        
        jdat <- rbind(jdat, dat.hscs)
        jdat$traj <- jtraj
        jdat$colcodetraj <- jcolcodetraj
        # jdat$colcode <- jcolcode
        return(jdat)
      }) %>%
        bind_rows()
      
      m.traj.fit.smooth <- ggplot(mapping = aes(x = ptime, y = signal, alpha = jalpha, group = traj)) + 
        geom_point(mapping = aes(color = colcode), data = dat.trajs.long) + 
        # stat_smooth() + 
        # geom_smooth(method = "lm", se = FALSE, formula = y ~ 1 + x + x ^ 2) +
        geom_smooth(mapping = aes(color = colcodetraj), 
                    method = "gam", se = FALSE, formula = y ~ s(x, bs = "cs", k = 4),
                    data = dat.commonroot) +
        scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                              guide = "legend") + 
        theme_bw() + 
        ggtitle(jtitle) + 
        theme(aspect.ratio=jratio, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
      print(m.traj.fit.smooth)
      
      m.traj.fit <- ggplot(mapping = aes(x = ptime, y = signal, alpha = jalpha, group = traj)) + 
        geom_point(mapping = aes(color = colcode), data = dat.trajs.long) + 
        # stat_smooth() + 
        # geom_smooth(method = "lm", se = FALSE, formula = y ~ 1 + x + x ^ 2) +
        geom_smooth(mapping = aes(color = colcodetraj), 
                    method = "gam", se = FALSE, formula = y ~ s(x, bs = "cs", k = 4),
                    data = dat.commonroot) +
        scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                              guide = "legend") + 
        theme_bw() + 
        ggtitle(paste(jmark, jctype)) + 
        theme(aspect.ratio=jratio, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
      print(m.traj.fit)
      
      m.traj.batch <- ggplot(dat.trajs.long, aes(x = ptime, y = signal, color = colcode, alpha = jalpha)) + 
        geom_point() + 
        ggtitle(paste(jmark, jctype)) + 
        scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                              guide = "legend") + 
        theme_bw() + 
        facet_wrap(~batch) + 
        theme(aspect.ratio=jratio, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
      print(m.traj.batch)
      
      m.traj.batch2 <- ggplot(dat.trajs.long, aes(x = ptime, y = signal, color = batch, alpha = jalpha)) + 
        geom_point(alpha = 0.25) + 
        geom_point() + 
        theme_bw() + 
        ggtitle(paste(jmark, jctype)) + 
        facet_wrap(~traj) + 
        theme(aspect.ratio=jratio, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
      print(m.traj.batch2)
      
    dev.off()
    }
}, mc.cores = length(jmarks))


