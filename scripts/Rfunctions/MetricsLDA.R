# Jake Yeung
# Date of Creation: 2019-01-04
# File: ~/projects/scChiC/scripts/Rfunctions/MetricsLDA.R
# Metrics for LDA

CalculateMetrics <- function(models, topics, verbose=TRUE){
  if (verbose) cat("calculate metrics:\n")
  metrics <- c("CaoJuan2009", "Deveaud2014")
  result <- data.frame(topics)
  for(m in metrics) {
    if (verbose) cat(sprintf("  %s...", m))
    if (! m %in% c("Griffiths2004", "CaoJuan2009", "Arun2010", "Deveaud2014")) {
      cat(" unknown!\n")
    } else {
      result[m] <- switch(m,
                          "Griffiths2004" = Griffiths2004(models, control),
                          "CaoJuan2009"   = CaoJuan2009(models),
                          "Arun2010"      = Arun2010(models, dtm),
                          "Deveaud2014"   = Deveaud2014(models),
                          NaN
      )
      if (verbose) cat(" done.\n")
    }
  }
  return(result)
}

#' @keywords internal
Griffiths2004 <- function(models, control) {
  # log-likelihoods (remove first burning stage)
  burnin  <- ifelse("burnin" %in% names(control), control$burnin, 0)
  logLiks <- lapply(models, function(model) {
    utils::tail(model@logLiks, n = length(model@logLiks) - burnin/control$keep)
    # model@logLiks[-(1 : (control$burnin/control$keep))]
  })
  # harmonic means for every model
  metrics <- sapply(logLiks, function(x) {
    # code is a little tricky, see explanation in [Ponweiser2012 p. 36]
    # ToDo: add variant without "Rmpfr"
    llMed <- stats::median(x)
    metric <- as.double(
      llMed - log( Rmpfr::mean( exp( -Rmpfr::mpfr(x, prec=2000L) + llMed )))
    )
    return(metric)
  })
  return(metrics)
}

#' @keywords internal
CaoJuan2009 <- function(models) {
  metrics <- sapply(models, function(model) {
    # topic-word matrix
    m1 <- exp(model@beta)
    # pair-wise cosine distance
    pairs <- utils::combn(nrow(m1), 2)
    cos.dist <- apply(pairs, 2, function(pair) {
      x <- m1[pair[1], ]
      y <- m1[pair[2], ]
      # dist <- lsa::cosine(x, y)
      dist <- crossprod(x, y) / sqrt(crossprod(x) * crossprod(y))
      return(dist)
    })
    # metric
    metric <- sum(cos.dist) / (model@k*(model@k-1)/2)
    return(metric)
  })
  return(metrics)
}

#' @keywords internal
Arun2010 <- function(models, dtm) {
  # length of documents (count of words)
  len <- slam::row_sums(dtm)
  # evaluate metrics
  metrics <- sapply(models, FUN = function(model) {
    # matrix M1 topic-word
    m1 <- exp(model@beta) # rowSums(m1) == 1
    m1.svd <- svd(m1)
    cm1 <- as.matrix(m1.svd$d)
    # matrix M2 document-topic
    m2   <- model@gamma   # rowSums(m2) == 1
    cm2  <- len %*% m2    # crossprod(len, m2)
    norm <- norm(as.matrix(len), type="m")
    cm2  <- as.vector(cm2 / norm)
    # symmetric Kullback-Leibler divergence
    divergence <- sum(cm1*log(cm1/cm2)) + sum(cm2*log(cm2/cm1))
    return ( divergence )
  })
  return(metrics)
}

#' Deveaud2014
#' @keywords internal
Deveaud2014 <- function(models) {
  metrics <- sapply(models, function(model) {
    ### original version
    # topic-word matrix
    m1 <- exp(model@beta)
    # prevent NaN
    if (any(m1 == 0)) { m1 <- m1 + .Machine$double.xmin }
    # pair-wise Jensen-Shannon divergence
    pairs  <- utils::combn(nrow(m1), 2)
    jsd <- apply(pairs, 2, function(pair) {
      x <- m1[pair[1], ]
      y <- m1[pair[2], ]
      ### standard Jensen-Shannon divergence
      # m <- (x + y) / 2
      # jsd <- 0.5 * sum(x*log(x/m)) + 0.5 * sum(y*log(y/m))
      ### divergence by Deveaud2014
      jsd <- 0.5 * sum(x*log(x/y)) + 0.5 * sum(y*log(y/x))
      return(jsd)
    })
    
    #     ### optimized version
    #     m1   <- model@beta
    #     m1.e <- exp(model@beta)
    #     pairs  <- utils::combn(nrow(m1), 2)
    #     jsd <- apply(pairs, 2, function(pair) {
    #       x   <- m1[pair[1], ]
    #       y   <- m1[pair[2], ]
    #       x.e <- m1.e[pair[1], ]
    #       y.e <- m1.e[pair[2], ]
    #       jsd <- ( sum(x.e*(x-y)) + sum(y.e*(y-x)) ) / 2
    #       return(jsd)
    #     })
    
    # metric
    metric <- sum(jsd) / (model@k*(model@k-1))
    return(metric)
  })
  return(metrics)
}


#' FindTopicsNumber_plot
#'
#' Support function to analyze optimal topic number. Use output of the
#' \code{\link{FindTopicsNumber}} function.
#'
#' @param values Data-frame with first column named `topics` and other columns
#'   are values of metrics.
#'
#' @examples
#' \dontrun{
#'
#' library(topicmodels)
#' data("AssociatedPress", package="topicmodels")
#' dtm <- AssociatedPress[1:10, ]
#' optimal.topics <- FindTopicsNumber(dtm, topics = 2:10,
#'   metrics = c("Arun2010", "CaoJuan2009", "Griffiths2004")
#' )
#' FindTopicsNumber_plot(optimal.topics)
#' }
#'
#' @export
#' @import ggplot2
FindTopicsNumber_plot <- function(values) {
  # normalize to [0,1]
  columns <- base::subset(values, select = 2:ncol(values))
  values <- base::data.frame(
    values["topics"],
    base::apply(columns, 2, function(column) {
      scales::rescale(column, to = c(0, 1), from = range(column))
    })
  )
  
  # melt
  values <- reshape2::melt(values, id.vars = "topics", na.rm = TRUE)
  
  # separate max-arg & min-arg metrics
  values$group <- values$variable %in% c("Griffiths2004", "Deveaud2014")
  values$group <- base::factor(
    values$group,
    levels = c(FALSE, TRUE),
    labels = c("minimize", "maximize")
  )
  
  # standart plot
  p <- ggplot(values, aes_string(x = "topics", y = "value", group = "variable"))
  p <- p + geom_line()
  p <- p + geom_point(aes_string(shape = "variable"), size = 3)
  p <- p + guides(size = FALSE, shape = guide_legend(title = "metrics:"))
  p <- p + scale_x_continuous(breaks = values$topics)
  p <- p + labs(x = "number of topics", y = NULL)
  
  # separate in two parts
  p <- p + facet_grid(group ~ .)
  
  # style
  # p <- p + theme_bw(base_size = 14, base_family = "") %+replace% theme(
  p <- p + theme_bw() %+replace% theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey70"),
    panel.grid.minor.x = element_blank(),
    legend.key = element_blank(),
    strip.text.y = element_text(angle = 90)
  )
  
  # move strip block to left side
  g <- ggplotGrob(p)
  g$layout[g$layout$name == "strip-right", c("l", "r")] <- 3
  grid::grid.newpage()
  grid::grid.draw(g)
  
  # return(p)
}