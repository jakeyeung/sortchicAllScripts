GetTrajMixed <- function(inf.springtraj = "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-15.RData", 
                         inf.traj =  "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"){
  assertthat::assert_that(file.exists(inf.springtraj))
  load(inf.springtraj, v=T)
  assertthat::assert_that(file.exists(inf.traj))
  load(inf.traj, v=T)
  
  trajs.mixed <- c(trajs["H3K4me1"], trajs.spring[c("H3K4me3", "H3K27me3", "H3K9me3")])
  # rename bcell trajectory as lymphoid
  names(trajs.mixed$H3K4me1)[which(names(trajs.mixed$H3K4me1) == "lymphoid")] <- "lymphod.naive"
  names(trajs.mixed$H3K4me1)[which(names(trajs.mixed$H3K4me1) == "bcell")] <- "lymphoid"
  
  # handle dat.umap.long
  dat.umap.mixed <- list(H3K4me1 = dat.umap.long.trajs[["H3K4me1"]], 
                         H3K4me3 = dat.trajs.long %>% filter(mark == "H3K4me3"),
                         H3K27me3 = dat.trajs.long %>% filter(mark == "H3K27me3"), 
                         H3K9me3 = dat.trajs.long %>% filter(mark == "H3K9me3"))
  return(list(trajs.mixed = trajs.mixed, dat.umap.mixed = dat.umap.mixed))
}

SwapTrajName <- function(traj.vec){
  # for H3K4me1 there is bcell and tcell and lymphoid
  # rename bcell as lymphoid, lymphoid as lymphoid.naive
  # as variance_over_pseudotime_standarderror.R
  
  # first rename lymphoid as lymphoid.naive
  traj.vec <- gsub("lymphoid", "lymphoid.naive", traj.vec)
  # rename bcell as lymphoid
  traj.vec <- gsub("bcell", "lymphoid", traj.vec)
  return(traj.vec)
}

RenameTraj <- function(x){
  if (x == "eryth"){
    return("Erythryoid")
  } else if (x == "granu") {
    return("Myeloid")
  } else if (x == "lymphoid"){
    return("Lymphoid") 
  }
}

RenameTraj.rev <- function(x){
  if (x == "Erythryoid"){
    return("eryth")
  } else if (x == "Myeloid") {
    return("granu")
  } else if (x == "Lymphoid"){
    return("lymphoid") 
  }
}


GetTopGenes <- function(dat.long, jcelltype, jtopic, topn = 20){
  jsub <- dat.long %>% 
    rowwise() %>%
    mutate(CellType2 = ifelse(CellType == "nucleate_erythrocyte", TRUE, FALSE)) %>% 
    group_by(CellType2, Gene_Name) %>%
    summarise(zscore = mean(logFPKM)) %>%
    group_by(Gene_Name) %>% 
    summarise(zscore.diff = diff(range(zscore)))
  # jsub.hash <- hash(jsub$Gene_Name, jsub$zscore.diff)
  # get top genes
  jsub.top <- top.peaks.annotated %>% 
    filter(topic == jtopic & rnk <= 100) %>%
    group_by(topic, SYMBOL) %>%
    filter(beta == max(beta)) %>%
    arrange(rnk)
  jsub.top <- left_join(jsub.top, jsub, by = c("SYMBOL" = "Gene_Name")) %>% 
    arrange(desc(zscore.diff))
  jsub.sub <- jsub.top[1:topn, ]
  # genes <- jsub.top$SYMBOL[1:topn]
  return(jsub.sub)
}

# do Principal Curves on filtered genes
InferTrajOnUmap <- function(dat, cname = "granu", init.on = "umap2", flip.lambda = FALSE, return.obj = FALSE){
  inmat <- subset(dat[which(dat[[cname]] == TRUE), ], select = c(umap1, umap2, cell))
  rownames(inmat) <- inmat$cell
  inmat$cell <- NULL
  # order from origin outwards
  inmat <- inmat[order(inmat$umap1 ^ 2 + inmat$umap2 ^ 2), ]
  # create initial ordering
  if (init.on == "umap1"){
    start.x <- inmat[["umap1"]]
    start.y <- rep(mean(inmat[["umap2"]]), nrow(inmat))
  } else if (init.on == "umap2"){
    start.x <- rep(mean(inmat[["umap1"]], nrow(inmat)))
    start.y <- inmat[["umap2"]]
  } else {
    warning(paste("init on must be umap1 or umap2", init.on))
  }
  startvals <- as.matrix(data.frame(umap1 = start.x, umap2 = start.y))
  proj <- principal_curve(as.matrix(inmat))
  dat.pca.proj <- as.data.frame(proj$s)
  dat.pca.proj$cell <- rownames(proj$s)
  # scale lambda from 0 to 1
  dat.pca.proj$lambda <- proj$lambda / max(proj$lambda)
  dat.pca.proj <- dat.pca.proj %>% arrange(lambda)
  if (flip.lambda){
    dat.pca.proj$lambda <- 1 - dat.pca.proj$lambda
  }
  # if end of lambda is near origin, then switch
  # endpoint <- dat.pca.proj[, 1]$umap1 *
  if (!return.obj){
    return(dat.pca.proj)
  } else {
    return(list(dat.pca.proj = dat.pca.proj, pc.obj = proj))
  }
}


CatActivityCountsTraj <- function(dat.umap.long, mat.imput, act.long.merged, traj, jpeak, jscale = TRUE, jscale.all = FALSE){
  peak.counts <- data.frame(peak = jpeak, 
                            cell = colnames(mat.imput),
                            counts = log10(mat.imput[jpeak, ]))
  if (jscale.all){
    peak.counts$counts <- scale(peak.counts$counts, center = TRUE, scale = TRUE)
  }
  act.counts <- subset(act.long.merged, motif == jmotif)
  if (jscale.all){
    act.counts$activity <- scale(act.counts$activity, center = TRUE, scale = TRUE)
  }
  dat.counts <- left_join(peak.counts, act.counts)
  act.counts$counts <- act.counts$activity
  act.counts <- left_join(act.counts, traj)
  peak.counts <- left_join(peak.counts, traj)
  # show umap, ad cebpb activity 
  # make dat.merge for activity and counts plot 
  dat.umap.long.cat <- left_join(dat.umap.long, act.counts %>% dplyr::select(cell, counts))
  dat.umap.long.cat <- left_join(dat.umap.long.cat, peak.counts %>% dplyr::rename(peak.counts = counts) %>% dplyr::select(cell, peak.counts))
  dat.merge <- bind_rows(act.counts %>% filter(!is.na(lambda)) %>% mutate(assay = "activity"), 
                         peak.counts %>% filter(!is.na(lambda)) %>% mutate(assay = "counts"))
  if (jscale){
    dat.merge <- dat.merge %>%
      group_by(assay) %>%
      mutate(counts = scale(counts))
  }
  return(list(dat.umap.long.cat=dat.umap.long.cat, dat.merge=dat.merge))
}