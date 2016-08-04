weights.meta <- function(object,
                         comb.fixed = object$comb.fixed,
                         comb.random = object$comb.random,
                         ...) {
  
  ##
  ##
  ## (1) Check for meta object
  ##
  ##
  chkclass(object, "meta")
  
  
  if (!(comb.fixed | comb.random)) {
    warning("Information missing which weights should be calculated (see arguments 'comb.fixed' and 'comb.random')")
    return(NULL)
  }

  w.fixed  <- object$w.fixed
  w.random <- object$w.random
  ##
  p.fixed  <- 100 * w.fixed  / sum(w.fixed,  na.rm = TRUE)
  p.random <- 100 * w.random / sum(w.random, na.rm = TRUE)

  res <- data.frame(w.fixed, p.fixed, w.random, p.random)
  ##
  rownames(res) <- object$studlab
  ##
  if (!comb.fixed) {
    res$w.fixed <- NULL
    res$p.fixed <- NULL
  }
  ##
  if (!comb.random) {
    res$w.random <- NULL
    res$p.random <- NULL
  }
  
  res
}
