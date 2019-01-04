#' Calculate absolute and percentage weights for meta-analysis
#' 
#' @description
#' This function returns a data frame containing information on
#' absolute and percentage weights of individual studies contributing
#' to fixed effect and random effects meta-analysis.
#' 
#' @param object An object of class \code{meta}.
#' @param comb.fixed A logical indicating whether absolute and
#'   percentage weights from the fixed effect model should be
#'   calculated.
#' @param comb.random A logical indicating whether absolute and
#'   percentage weights from the random effects model should be
#'   calculated.
#' @param \dots Additional arguments (ignored at the moment).
#' 
#' @return
#' A data frame with the following variables is returned:
#' \tabular{lll}{
#' \bold{Variable} \tab \bold{Definition} \tab \bold{Condition} \cr
#' w.fixed \tab absolute weights in fixed effect model \tab (if
#'   \code{comb.fixed = TRUE}) \cr
#' p.fixed \tab percentage weights in fixed effect model \tab (if
#'   \code{comb.fixed = TRUE}) \cr
#' w.random \tab absolute weights in random effects model \tab (if
#'   \code{comb.random = TRUE}) \cr
#' p.random \tab percentage weights in random effects model \tab (if
#'   \code{comb.random = TRUE})
#' }
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{metagen}}
#' 
#' @examples
#' data(Fleiss93cont)
#' # Do meta-analysis (fixed effect and random effects model)
#' #
#' meta1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, study,
#'                   data = Fleiss93cont, sm = "SMD")
#' 
#' # Print weights for fixed effect and random effects meta-analysis
#' #
#' weights(meta1)
#' 
#' # Do meta-analysis (only random effects model)
#' #
#' meta2 <- update(meta1, comb.fixed = FALSE)
#' 
#' # Print weights for random effects meta-analysis
#' #
#' weights(meta2)
#' 
#' # Print weights for fixed effect and random effects meta-analysis
#' #
#' weights(meta2, comb.fixed = TRUE)
#' 
#' @export
#' @export weights.meta


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
