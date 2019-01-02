
FitGlmRow <- function(row, pseudo, size, returnobj=FALSE){
  # use Offset by size of library
  # https://stats.stackexchange.com/questions/66791/where-does-the-offset-go-in-poisson-negative-binomial-regression
  # fit GLM for a row of a sparse matrix, should save some space?
  if (!is.null(nrow(row))){
    # probably a matrix of many rows, sum them up
    print(paste("Merging", nrow(row), "rows"))
    row <- Matrix::colSums(row)
  }
  dat <- data.frame(pseudo = pseudo,
                    counts = row,
                    size = size)
  # m1.pois <- glm(counts ~ 1 + pseudo, data = dat, family = "poisson")
  m1.pois <- glm(counts ~ 1 + pseudo + offset(log(size)), data = dat, family = "poisson")
  # m1.pois.null <- glm(counts ~ 1, data = dat, family = "poisson")
  # compare <- 2*logLik(m1.pois) - 2*logLik(m1.pois.null)
  pval <- pchisq(m1.pois$deviance, df=m1.pois$df.residual, lower.tail=FALSE)
  if (!returnobj){
    return(list(int = coef(m1.pois)[[1]], pseudo = coef(m1.pois)[[2]], pval = pval))
  } else {
    return(m1.pois)  # takes more memory??
  }
}