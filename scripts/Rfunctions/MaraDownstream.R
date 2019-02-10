# Jake Yeung
# Date of Creation: 2019-02-10
# File: ~/projects/scchic/scripts/Rfunctions/MaraDownstream.R
# MARA downstream functions

LoadMARA <- function(mdir, bc.inf = "data/barcode_summaries/barcodes/maya_384NLA.bc"){
  act.mat <- fread(file.path(mdir, "Activities"), header = FALSE)
  se.mat <- fread(file.path(mdir, "StandardError"), header = FALSE)
  cnames <- unlist(fread(file.path(mdir, "Colnames"), header = FALSE), use.names = FALSE)
  zscores <- fread(file.path(mdir, "Zscores"), header = FALSE)
  colnames(zscores) <- c("motif", "zscore")
  zscores <- zscores %>% arrange(desc(zscore))
  
  barcodes <- read.table(bc.inf, col.names = "Barcodes", stringsAsFactors = FALSE)
  cellhash <- hash(rownames(barcodes), unlist(barcodes))
  cellhash.bc <- hash(unlist(barcodes), paste("cell", rownames(barcodes), sep = ""))
  
  cnames.new <- sapply(cnames, function(x){
    jmark <- strsplit(x, "\\.")[[1]][[1]]
    jtiss <- strsplit(x, "\\.")[[1]][[2]]
    jmouse <- strsplit(x, "\\.")[[1]][[3]]
    jrep <- paste0("rep", strsplit(x, "\\.")[[1]][[4]])
    jcell <- cellhash.bc[[strsplit(x, "\\.")[[1]][[6]]]]
    xnew <- paste(jtiss, jmark, jmouse, jrep, jcell, sep = "_")
  })
  
  colnames(act.mat) <- c("motif", cnames.new)
  colnames(se.mat) <- c("motif", cnames.new)
  
  act.long <- tidyr::gather(act.mat, -motif, key = "cell", value = "activity")
  
  return(list(act.mat = act.mat, se.mat = se.mat, zscores = zscores, act.long = act.long))
}