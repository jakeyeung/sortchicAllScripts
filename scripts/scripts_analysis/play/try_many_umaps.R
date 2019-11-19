# Jake Yeung
# Date of Creation: 2019-03-24
# File: ~/projects/scchic/scripts/scripts_analysis/play/try_many_umaps.R
# Test different neighbors and umaps


nn.change <- seq(17, 23)
jmark <- "H3K4me3"
# clstr.hash.new <- DoLouvain(topics.mat.new.lst[[1]], custom.settings.louv.new.lst, dat.umap.long = NULL)
custom.settings.louv.new.lst <- lapply(nn.louv.new, function(x) GetUmapSettings(x, jmetric.louv, jmindist.louv, jseed.louv))

pdf(paste0("~/data/scchic/many_umaps.", jmark, ".pdf"))

for (nn in nn.change){
  nn.new <- c("H3K4me1" = nn, "H3K4me3" = nn, "H3K27me3" = nn, "H3K9me3" = nn)
  print(nn.new)
  custom.settings.new.lst <- lapply(nn.new, function(x) GetUmapSettings(x, jmetric.louv, 0.01, jseed.louv))
  jsets <- GetUmapSettings(nn, jmetric.louv, jmindist.louv, jseed.louv)
  
  # custom.settings.new <- GetUmapSettings(nn=nn.louv[[1]], jmetric=jmetric.louv, jmindist=jmindist.louv, seed=jseed.louv)
  
  # do UMAP on new settings
  dat.umap.new <- umap(topics.mat.new.lst[[jmark]], custom.settings.new.lst[[jmark]])
  
  #  assign cluster to umap
    # dat.umap.new <- dat.umap.new.lst[[jmark]]
    dat.umap.long.new <- data.frame(umap1 = dat.umap.new$layout[, 1], umap2 = dat.umap.new$layout[, 2], cell = rownames(dat.umap.new$layout), stringsAsFactors = FALSE)
    dat.umap.long.new$louvain <- as.character(sapply(dat.umap.long.new$cell, function(x) clstr.hash.new.lst[[jmark]][[x]]))
    dat.umap.long.new$mark <- jmark
  
  m1 <- ggplot(dat.umap.long.new, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(nn.new[[1]])
  print(m1)
}
dev.off()