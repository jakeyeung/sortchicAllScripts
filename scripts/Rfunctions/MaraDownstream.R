# Jake Yeung
# Date of Creation: 2019-02-10
# File: ~/projects/scchic/scripts/Rfunctions/MaraDownstream.R
# MARA downstream functions

GetCellHashBC <- function(inf.experinames = "data/outputs/experinames.txt", bc.inf = "data/barcode_summaries/barcodes/maya_384NLA.bc"){
  experis <- read.table(inf.experinames, stringsAsFactors = FALSE)
  experihash <- hash::hash(sapply(experis$V1, function(x) strsplit(x, "_")[[1]][[1]]), sapply(experis$V1, function(x) paste(strsplit(x, "_")[[1]][2:3], collapse="_")))
  switch.rep.hash <- GetRepSwitchHash(experihash)
  
  barcodes <- read.table(bc.inf, col.names = "Barcodes", stringsAsFactors = FALSE)
  cellhash <- hash::hash(rownames(barcodes), unlist(barcodes))
  cellhash.bc <- hash::hash(unlist(barcodes), paste("cell", rownames(barcodes), sep = "")) 
}

SwitchColnames <- function(cnames, inf.experinames = "data/outputs/experinames.txt", bc.inf = "data/barcode_summaries/barcodes/maya_384NLA.bc", rep.prefix = "", jsplit = "\\."){
  source("scripts/Rfunctions/MatchCellNameToSample.R")
  experis <- read.table(inf.experinames, stringsAsFactors = FALSE)
  experihash <- hash::hash(sapply(experis$V1, function(x) strsplit(x, "_")[[1]][[1]]), sapply(experis$V1, function(x) paste(strsplit(x, "_")[[1]][2:3], collapse="_")))
  switch.rep.hash <- GetRepSwitchHash(experihash)
  
  barcodes <- read.table(bc.inf, col.names = "Barcodes", stringsAsFactors = FALSE)
  cellhash <- hash::hash(rownames(barcodes), unlist(barcodes))
  cellhash.bc <- hash::hash(unlist(barcodes), paste("cell", rownames(barcodes), sep = ""))
  
  cnames.new <- sapply(cnames, function(x){
    jmark <- strsplit(x, jsplit)[[1]][[1]]
    jtiss <- strsplit(x, jsplit)[[1]][[2]]
    jmouse <- strsplit(x, jsplit)[[1]][[3]]
    jrep <- paste0(rep.prefix, strsplit(x, jsplit)[[1]][[4]])
    jcell <- cellhash.bc[[strsplit(x, jsplit)[[1]][[6]]]]
    if (is.null(jcell)){
      warning(paste(strsplit(x, jsplit)[[1]][[6]], "not in hash, returning NA"))
      xnew <- NA
    } else {
      xnew <- paste(jtiss, jmark, jmouse, jrep, jcell, sep = "_")
    }
  })
  # need to source
  # source("scripts/Rfunctions/MatchCellNameToSample.R")
  cnames.new <- sapply(cnames.new, function(x) SwapRepNameInCell(x, switch.rep.hash))
  return(cnames.new)
}

LoadMARA <- function(mdir, bc.inf = "data/barcode_summaries/barcodes/maya_384NLA.bc", fix.tech.rep = FALSE, rep.prefix = "rep", swap.tech.rep = NULL, filesuffix = "", make.cnames = TRUE){
  act.mat <- fread(file.path(mdir, paste0("Activities", filesuffix)), header = FALSE)
  se.mat <- fread(file.path(mdir, paste0("StandardError", filesuffix)), header = FALSE)
  cnames <- unlist(fread(file.path(mdir, paste0("Colnames", filesuffix)), header = FALSE), use.names = FALSE)
  
  zscores <- fread(file.path(mdir, "Zscores"), header = FALSE)
  colnames(zscores) <- c("motif", "zscore")
  zscores <- zscores %>% arrange(desc(zscore))

  if (make.cnames){
    barcodes <- read.table(bc.inf, col.names = "Barcodes", stringsAsFactors = FALSE)
    cellhash <- hash(rownames(barcodes), unlist(barcodes))
    cellhash.bc <- hash(unlist(barcodes), paste("cell", rownames(barcodes), sep = ""))
    cnames.new <- sapply(cnames, function(x){
      jmark <- strsplit(x, "\\.")[[1]][[1]]
      jtiss <- strsplit(x, "\\.")[[1]][[2]]
      jmouse <- strsplit(x, "\\.")[[1]][[3]]
      jrep <- paste0(rep.prefix, strsplit(x, "\\.")[[1]][[4]])
      jcell <- cellhash.bc[[strsplit(x, "\\.")[[1]][[6]]]]
      xnew <- paste(jtiss, jmark, jmouse, jrep, jcell, sep = "_")
    })
    if (fix.tech.rep){
      for (mouse in c("m1", "m2")){
        # get tech reps for each mouse
        cnames.i <- grep(mouse, cnames.new)
        cnames.sub <- cnames.new[cnames.i]
        techreps <- unique(sapply(cnames.sub, function(x) strsplit(x, "_")[[1]][[4]]))  # repS9 and repS11
        techreps.number <- sapply(techreps, function(techrep) readr::parse_number(techrep))  # 9, 11
        # convert number to rank
        techreps.rank <- rank(techreps.number)
        techreps.new <- paste("rep", techreps.rank, sep = "")
        cnames.sub.new <- cnames.sub  # init
        for (i in seq(length(techreps.new))){
          techrep <- techreps[[i]]
          techrep.new <- techreps.new[[i]]
          cnames.sub.new <- gsub(paste0("_", techrep, "_"), paste0("_", techrep.new, "_"), cnames.sub.new)
        }
        # replace repS11 with rep1
        # cnames.sub.new <- mapply(function(techrep, techrep.new, cname) return(gsub(techrep, techrep.new, cname)), techreps, techreps.new, cnames.sub)
        cnames.new[cnames.i] <- cnames.sub.new
      }
    }
    
    if (!is.null(swap.tech.rep)){
      # need to source
      # source("scripts/Rfunctions/MatchCellNameToSample.R")
      cnames.new <- sapply(cnames.new, function(x) SwapRepNameInCell(x, swap.tech.rep))
    }
  } else {
    cnames.new <- cnames
  }
  
  
  colnames.lst <- data.frame(old.colnames = cnames, new.colnames = cnames.new)
  
  colnames(act.mat) <- c("motif", cnames.new)
  colnames(se.mat) <- c("motif", cnames.new)
  
  act.long <- tidyr::gather(act.mat, -motif, key = "cell", value = "activity")
  
  return(list(act.mat = act.mat, se.mat = se.mat, zscores = zscores, act.long = act.long, colnames.lst = colnames.lst))
}

SwapTechRep <- function(x){
  # swap tech rep 
}

LoadLDA <- function(inf){
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  out.lda <- ChooseBestLDA(out.lda)
  (kchoose <- out.lda@k)
  return(out.lda)
}

DoUmapFromLDA <- function(out.lda, custom.settings){
  tm.result <- topicmodels::posterior(out.lda)
  topics.mat <- tm.result$topics
  # terms.mat <- tm.result$terms
  
  # custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)
  # custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)
  
  dat.umap <- umap(topics.mat, config = custom.settings)
  rownames(dat.umap$layout) <- rownames(topics.mat)
  dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout))
  return(dat.umap.long)
}

PlotMotifInUmap <- function(jmotif, dat.merged, zscores, jmark, jsize = 1, colvec = c("blue", "gray95", "red")){
  dat.sub <- subset(dat.merged, motif == jmotif)
  jzscore <- signif(subset(zscores, motif == jmotif)$zscore, digits = 2)
  jtitle <- paste(jmotif, "zscore:", jzscore, jmark)
  # add ordering statistic
  # https://stackoverflow.com/questions/11805354/geom-point-put-overlapping-points-with-highest-values-on-top-of-others
  dat.sub <- dat.sub %>%
    group_by(motif) %>%
    mutate(orderrank = rank(activity, ties.method = "first")) %>%
    arrange(orderrank)
  m <- ggplot(dat.sub, aes(x = umap1, y = umap2, color = activity, order = orderrank)) + geom_point(size = jsize) +
  # m <- ggplot(dat.sub, aes(x = umap1, y = umap2, color = activity)) + geom_point(size = jsize) + 
    ggtitle(jtitle) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  if (length(colvec) == 3){
    midpt <- min(dat.sub$activity) + (max(dat.sub$activity) - min(dat.sub$activity)) / 2
    m <- m + 
      scale_color_gradient2(low = colvec[[1]], mid = colvec[[2]], high = muted(colvec[[3]]), midpoint = midpt)
  } else if (length(colvec) == 2){
    m <- m + 
      scale_color_gradient(low = colvec[[1]], high = muted(colvec[[2]]), space = "Lab")
  } else {
    warning(paste("Colvec must be length 2 or 3"))
    print(colvec)
  }
  return(m)
}

