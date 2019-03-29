# Jake Yeung
# Date of Creation: 2019-03-27
# File: ~/projects/scchic/scripts/Rfunctions/MatchCellNameToSample.R
# Match samples properly 

# go to /hpc/hub_oudenaarden/avo/scChiC and annotate combine_bin_files_BM_H3K4me1.R for all marks
# find match between cell name to bam name

SwapRepNameInCell <- function(x, rep.switch.hash){
  # BM_H3K4me1_m1_S9_cell152 -> rep1 or whatever
  cellname <- strsplit(x, "_")[[1]][[5]]
  jkey <- paste(strsplit(x, "_")[[1]][1:4], collapse = "_")
  jval <- rep.switch.hash[[jkey]]
  cellname.new <- paste(jval, cellname, sep = "_")
  return(cellname.new)
}

GetRepSwitchHash <- function(experihash){
  jkeys <- hash::keys(experihash)
  jvals <- hash::values(experihash)
  # cname is in form:
  # xnew <- paste(jtiss, jmark, jmouse, jrep, jcell, sep = "_")
  rep.switch.lst <- list()
  for (i in seq(length(jvals))){
    jval <- jvals[[i]]
    jkey <- jkeys[[i]]
    jmark <- strsplit(jkey, "-")[[1]][[4]]
    jmouse <- strsplit(jkey, "-")[[1]][[3]]
    jrepS <- strsplit(jval, "_")[[1]][[2]]
    jrep <- paste0("rep", strsplit(jkey, "-")[[1]][[5]])
    jtiss <- strsplit(jkey, "-")[[1]][[2]]
    
    experi.old <- paste(jtiss, jmark, jmouse, jrepS, sep = "_")  # ignoring the cell
    experi.new <- paste(jtiss, jmark, jmouse, jrep, sep = "_")
    rep.switch.lst[[experi.old]] <- experi.new
  }
  rep.swich.hash <- hash::hash(rep.switch.lst)
}



GetCellToSampName <- function(){
  # /hpc/hub_oudenaarden/avo/scChiC/combine_bin_files_BM_H3K9me3.R
  # 
  cell2samp <- list("BM_H3K9me3_m1_rep1" = "PZ-BM-m1-H3K9me3-2_AH3VGVBGX9_S1",
                    "BM_H3K9me3_m1_rep2" = "PZ-BM-m1-H3K9me3-1_H2GV2BGX9_S13",
                    "BM_H3K9me3_m2_rep1" = "PZ-BM-m2-H3K9me3-1_H2GV2BGX9_S17",
                    "BM_H3K9me3_m2_rep2" = "PZ-BM-m2-H3K9me3-2_AH3VGVBGX9_S5",
                    "BM_H3K27me3_m1_rep1" = "PZ-BM-m1-H3K27me3-2_AH3VGVBGX9_S2",
                    "BM_H3K27me3_m1_rep2" = "",
                    "BM_H3K27me3_m2_rep1" = "",
                    "BM_H3K27me3_m2_rep2" = "")
  
       
}

