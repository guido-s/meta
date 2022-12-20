#' Calculate absolute and percentage weights for meta-analysis
#' 
#' @description
#' This function returns a data frame containing information on
#' absolute and percentage weights of individual studies contributing
#' to common effect and random effects meta-analysis.
#' 
#' @param object An object of class \code{meta}.
#' @param common A logical indicating whether absolute and
#'   percentage weights from the common effect model should be
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
#' w.common \tab absolute weights in common effect model \tab (if
#'   \code{common = TRUE}) \cr
#' p.common \tab percentage weights in common effect model \tab (if
#'   \code{common = TRUE}) \cr
#' w.random \tab absolute weights in random effects model \tab (if
#'   \code{random = TRUE}) \cr
#' p.random \tab percentage weights in random effects model \tab (if
#'   \code{random = TRUE})
#' }
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{metagen}}
#' 
#' @examples
#' data(Fleiss1993cont)
#' # Do meta-analysis (common effect and random effects model)
#' #
#' meta1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, studlab = paste(study, year), sm = "SMD")
#' 
#' # Print weights for common effect and random effects meta-analysis
#' #
#' weights(meta1)
#' 
#' # Do meta-analysis (only random effects model)
#' #
#' meta2 <- update(meta1, common = FALSE)
#' 
#' # Print weights for random effects meta-analysis
#' #
#' weights(meta2)
#' 
#' # Print weights for common effect and random effects meta-analysis
#' #
#' weights(meta2, common = TRUE)
#' 
#' @method weights meta
#' @export


weights.meta <- function(object,
                         common = object$common,
                         random = object$random,
                         warn.deprecated = gs("warn.deprecated"),
                         ...) {
  
  ##
  ##
  ## (1) Check for meta object
  ##
  ##
  chkclass(object, "meta")
  object <- update(object)
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  missing.common <- missing(common)
  common <-
    deprecated(common, missing.common, args, "comb.fixed",
               warn.deprecated)
  common <-
    deprecated(common, missing.common, args, "fixed",
               warn.deprecated)
  ##
  random <-
    deprecated(random, missing(random), args, "comb.random",
               warn.deprecated)
  
  
  if (!(common | random)) {
    warning("Information missing which weights should be ",
            "calculated (see arguments 'common' and 'random')")
    return(NULL)
  }

  w.common  <- object$w.common
  w.random <- object$w.random
  ##
  p.common <- 100 * w.common / sum(w.common, na.rm = TRUE)
  p.random <- 100 * w.random / sum(w.random, na.rm = TRUE)

  res <- data.frame(w.common, p.common, w.random, p.random)
  ##
  rownames(res) <- object$studlab
  ##
  if (!common) {
    res$w.common <- NULL
    res$p.common <- NULL
  }
  ##
  if (!random) {
    res$w.random <- NULL
    res$p.random <- NULL
  }
  
  res
}
