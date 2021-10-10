#' Calculate absolute and percentage weights for meta-analysis
#' 
#' @description
#' This function returns a data frame containing information on
#' absolute and percentage weights of individual studies contributing
#' to fixed effect and random effects meta-analysis.
#' 
#' @param object An object of class \code{meta}.
#' @param fixed A logical indicating whether absolute and
#'   percentage weights from the fixed effect model should be
#'   calculated.
#' @param random A logical indicating whether absolute and
#'   percentage weights from the random effects model should be
#'   calculated.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#' 
#' @return
#' A data frame with the following variables is returned:
#' \tabular{lll}{
#' \bold{Variable} \tab \bold{Definition} \tab \bold{Condition} \cr
#' w.fixed \tab absolute weights in fixed effect model \tab (if
#'   \code{fixed = TRUE}) \cr
#' p.fixed \tab percentage weights in fixed effect model \tab (if
#'   \code{fixed = TRUE}) \cr
#' w.random \tab absolute weights in random effects model \tab (if
#'   \code{random = TRUE}) \cr
#' p.random \tab percentage weights in random effects model \tab (if
#'   \code{random = TRUE})
#' }
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{metagen}}
#' 
#' @examples
#' data(Fleiss1993cont)
#' # Do meta-analysis (fixed effect and random effects model)
#' #
#' meta1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'                   data = Fleiss1993cont, sm = "SMD",
#'                   studlab = paste(study, year))
#' 
#' # Print weights for fixed effect and random effects meta-analysis
#' #
#' weights(meta1)
#' 
#' # Do meta-analysis (only random effects model)
#' #
#' meta2 <- update(meta1, fixed = FALSE)
#' 
#' # Print weights for random effects meta-analysis
#' #
#' weights(meta2)
#' 
#' # Print weights for fixed effect and random effects meta-analysis
#' #
#' weights(meta2, fixed = TRUE)
#' 
#' @method weights meta
#' @export


weights.meta <- function(object,
                         fixed = object$fixed,
                         random = object$random,
                         warn.deprecated = gs("warn.deprecated"),
                         ...) {
  
  ##
  ##
  ## (1) Check for meta object
  ##
  ##
  chkclass(object, "meta")
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  fixed <-
    deprecated(fixed, missing(fixed), args, "comb.fixed",
               warn.deprecated)
  ##
  random <-
    deprecated(random, missing(random), args, "comb.random",
               warn.deprecated)
  
  
  if (!(fixed | random)) {
    warning("Information missing which weights should be ",
            "calculated (see arguments 'fixed' and 'random')")
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
  if (!fixed) {
    res$w.fixed <- NULL
    res$p.fixed <- NULL
  }
  ##
  if (!random) {
    res$w.random <- NULL
    res$p.random <- NULL
  }
  
  res
}
