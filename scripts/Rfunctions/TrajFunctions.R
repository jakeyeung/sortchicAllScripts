
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
InferTrajOnUmap <- function(dat, cname = "granu", init.on = "umap1", flip.lambda = FALSE){
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
  return(dat.pca.proj)
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