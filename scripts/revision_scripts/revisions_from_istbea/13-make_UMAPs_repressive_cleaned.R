# Jake Yeung
# Date of Creation: 2022-02-17
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/13-make_UMAPs_repressive_cleaned.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


jmarksnew <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarksnew) <- jmarksnew
jmarksrepress <- c("k27me3", "k9me3"); names(jmarksrepress) <- jmarksrepress

# Load metas for cell types --------------------------------------------------------------

indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up"

dat.meta.lst <- lapply(jmarksnew, function(jmark){
  fread(file.path(indir.meta, paste0("metadata.", jmark, ".txt")))
})


# Update UMAPs with re-LDA'd repressive marks -----------------------------

indir.cleaned.meta <- "metadata_plate_experi_batch.cleaned.dynamicbins.k27me3.2022-02-17.txt"

dat.cleaned.meta.lst <- lapply(jmarksnew, function(jmark){
  if (!jmark %in% jmarksrepress){
    dat.meta <- dat.meta.lst[[jmark]]
  } else {
    inf.meta.cleaned <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata/repressive_cleaned/metadata_plate_experi_batch.cleaned.dynamicbins.", jmark, ".2022-02-17.txt")
    dat.meta.cleaned <- fread(inf.meta.cleaned) %>%
      dplyr::select(cell, umap1, umap2) 
    assertthat::assert_that(nrow(dat.meta.cleaned) > 0)
    dat.meta.old <- dat.meta.lst[[jmark]] %>%
      dplyr::select(-umap1, -umap2)
    dat.meta <- left_join(dat.meta.old, dat.meta.cleaned)
  }
  return(dat.meta)
})



# Show UMAP  --------------------------------------------------------------

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/repressive_cleaned_up"
outf <- file.path(outdir, paste0("repressive_cleaned_up.", Sys.Date(), ".pdf"))
pdf(outf, useDingbats = FALSE)

jalpha <- 0.5
m.lst <- lapply(jmarksnew, function(jmark){
  jdat <- dat.cleaned.meta.lst[[jmark]]
  ncells <- nrow(jdat)
  jsub <- jdat %>%
    dplyr::select(ctype.from.LL, colcode) 
  jsub <- jsub[!duplicated(jsub), ]
  
  ggplot(jdat, aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point(alpha = jalpha) + 
    theme_bw() + 
    ggtitle(jmark, paste0("ncells: ", ncells)) + 
    scale_color_identity( labels = jsub$ctype.from.LL, breaks = jsub$colcode,
                          guide = "legend")  + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

print(m.lst)


m2.lst <- lapply(jmarksnew, function(jmark){
  jdat <- dat.cleaned.meta.lst[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = ctype)) + 
    geom_point(alpha = jalpha) + 
    facet_wrap(~ctype.from.LL) + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m2.lst)




# Show cell enrichment ----------------------------------------------------

# DCs
# monocytes
# pDCs

jmark <- "k4me1"

for (jmark in jmarksnew){
  print(jmark)
  jdat <- dat.cleaned.meta.lst[[jmark]]
  jctypes <- sort(unique(jdat$ctype.from.LL))
  for (jctype in jctypes){
    print(jctype)
    jdat.highlight <- jdat %>%
      rowwise() %>%
      mutate(highlight.LL = ctype.from.LL == jctype,
             highlight = ctype == jctype) %>%
      arrange(highlight.LL)
    
    jcolcode <- unique(subset(jdat.highlight, ctype.from.LL == jctype)$colcode)
    
    m <- ggplot(jdat.highlight, aes(x = umap1, y = umap2, color = highlight.LL)) + 
      geom_point() + 
      theme_bw() + 
      scale_color_manual(values = c("grey80", jcolcode)) + 
      ggtitle(paste(jmark, jctype)) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
    
    m <- ggplot(jdat.highlight %>% arrange(highlight), aes(x = umap1, y = umap2, color = highlight)) + 
      geom_point() + 
      theme_bw() + 
      scale_color_manual(values = c("grey80", jcolcode)) + 
      ggtitle(paste(jmark, jctype)) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
  }
}
dev.off()




