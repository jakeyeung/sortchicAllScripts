# Jake Yeung
# Date of Creation: 2019-12-29
# File: ~/projects/scchic/scripts/rstudioserver_analysis/intestinal_all_merged/0-merge_scraped_with_noscraping.R
# Find the enterocytes by merging scraped with noscraping 

# mix k27me3 with k27me3

outdir <- "/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/qc_scraped_unscraped_merged"

# jmark <- "k27me3"
jmarks <- c("k4me1", "k36me3")

for (jmark in jmarks){
  print(jmark)
  inf1 <- paste0("/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/quality_controls_intestines/mat.Scraped.Unenriched.", jmark, ".TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.2019-12-23.rds")
  inf2 <- paste0("/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/quality_controls_intestines/mat.NoScraping.AllMerged.", jmark, ".TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.2019-12-23.rds")
  assertthat::assert_that(file.exists(inf1))
  assertthat::assert_that(file.exists(inf2))
  
  
  mat1 <- readRDS(inf1)
  mat2 <- readRDS(inf2)
  
  # get common rows
  rnames.common <- intersect(rownames(mat1), rownames(mat2))
  mat.merged <- cbind(mat1[rnames.common, ], mat2[rnames.common, ])
  print(lapply(list(mat1, mat2, mat.merged), dim))
  
  # save to output
  outf.base <- file.path(outdir, paste0("mat.ScrapeUnscrapeMerged.", jmark))
  if (!file.exists(paste0(outf.base, ".rds"))){
    saveRDS(mat.merged, paste0(outf.base, ".rds"))
  }
  writeMM(mat.merged, sparse = TRUE, file = paste0(outf.base, ".mm"))
  write.table(rownames(mat.merged), file = paste0(outf.base, ".rownames"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
  write.table(colnames(mat.merged), file = paste0(outf.base, ".colnames"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
  
}


# mix k9me3 with k9m3 
