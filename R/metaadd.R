#' Add pooled results from external analysis to meta-analysis
#' 
#' @description
#' Add pooled results from external analysis to an existing
#' meta-analysis object. This is useful, for example, to add results
#' from a Bayesian meta-analysis which is not implemented in R package
#' \bold{meta}.
#' 
#' @param x Meta-analysis object.
#' @param type A character string or vector indicating whether added
#'   results are from common effect, random effects model or
#'   prediction interval. Either \code{"common"}, \code{"random"} or
#'   \code{"prediction"}, can be abbreviated.
#' @param TE Pooled estimate(s).
#' @param lower Lower limit(s) of confidence or prediction interval.
#' @param upper Upper limit(s) of confidence or prediction interval.
#' @param statistic Test statistic(s).
#' @param pval P-value(s).
#' @param text A character string or vector used in printouts and
#'   forest plot to label the added results.
#' @param data An optional data frame containing the new results.
#' 
#' @details
#' In R package \bold{meta}, objects of class \code{"meta"} contain
#' results of both common effect and random effects
#' meta-analyses. This function enables the user to add the pooled
#' results of an additional analysis to an existing meta-analysis
#' object. This is useful, for example, to add the result of a
#' Bayesian meta-analysis.
#'
#' Arguments \code{TE}, \code{lower} and \code{upper} have to be
#' provided if \code{type = "common"} or \code{type = "random"}. For
#' \code{type = "prediction"}, only arguments \code{lower} and
#' \code{upper} have to be provided and other arguments are ignored.
#'
#' Note, R function \code{\link{metamerge}} can be used to add
#' meta-analysis results of another meta-analysis object (see
#' \code{\link{meta-object}}).
#' 
#' @return
#' An object of class \code{"meta"} with corresponding generic
#' functions (see \code{\link{meta-object}}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metagen}}, \code{\link{metamerge}}
#' 
#' @examples
#' data(Fleiss1993bin)
#' 
#' # Use REML estimator of tau2 (default)
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac, data = Fleiss1993bin,
#'   studlab = paste(study, year), sm = "OR",
#'   text.random = "Random effects model (REML)", overall.hetstat = FALSE)
#' 
#' # Use DerSimonian-Laird estimator of tau2
#' m2 <- update(m1, method.tau = "DL")
#' 
#' # Add random effects results of second meta-analysis
#' m12 <- metaadd(m1, type = "random", data = m2,
#'   TE = TE.random,
#'   lower = lower.random, upper = upper.random,
#'   statistic = statistic.random, pval = pval.random,
#'   text = "Random effects model (DL)")
#' m12
#' 
#' forest(m12)
#' 
#' @export metaadd


metaadd <- function(x, type,
                    TE, lower, upper, statistic = NA, pval = NA,
                    text, data = NULL) {
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  
  chkclass(x, "meta")
  res <- updateversion(x)
  ##
  missing.type <- missing(type)
  missing.TE <- missing(TE)
  missing.lower <- missing(lower)
  missing.upper <- missing(upper)
  missing.statistic <- missing(statistic)
  missing.pval <- missing(pval)
  missing.text <- missing(text)
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  if (nulldata)
    data <- sfsp
  ##
  ## Catch 'type', 'TE', 'lower', 'upper', 'statistic', and 'pval' from data:
  ##
  type <- catch("type", mc, data, sfsp)
  ##
  if (any(type %in% c("common", "random")))
    if (missing.TE)
      stop("Argument 'TE' must be provided for common effect or ",
           "random effects model.",
           call. = FALSE)
  ##
  if (missing.lower)
    stop("Argument 'lower' must be provided.",
         call. = FALSE)
  ##
  if (missing.upper)
    stop("Argument 'upper' must be provided.",
         call. = FALSE)
  ##
  TE <- catch("TE", mc, data, sfsp)
  lower <- catch("lower", mc, data, sfsp)
  upper <- catch("upper", mc, data, sfsp)
  statistic <- catch("statistic", mc, data, sfsp)
  pval <- catch("pval", mc, data, sfsp)
  text <- catch("text", mc, data, sfsp)
  ##
  type <- setchar(type, c("common", "random", "prediction"))
  if (length(type) == 1 & length(lower) > 1)
    type <- rep_len(type, length(lower))
  ##
  k.all <- length(type)
  ##
  if (!missing.TE) {
    chknumeric(TE)
    chklength(TE, k.all, name = "type")
  }
  ##
  chknumeric(lower)
  chklength(lower, k.all, name = "type")
  ##
  chknumeric(upper)
  chklength(upper, k.all, name = "type")
  ##
  if (!missing.statistic) {
    chknumeric(statistic)
    chklength(statistic, k.all, name = "type")
  }
  else
    statistic <- rep_len(NA, k.all)
  ##
  if (!missing.pval) {
    chknumeric(pval)
    chklength(pval, k.all, name = "type")
  }
  else
    pval <- rep_len(NA, k.all)
  ##
  if (!missing.text) {
    if (length(text) == 1 & k.all > 1)
      text <- rep_len(text, k.all)
    else
      chklength(text, k.all)
  }
  else
    text <- rep_len("Added result", k.all)
  ##
  for (i in seq_along(type)) {
    if (type[i] == "common") {
      res$common <- TRUE
      res$TE.common <- c(res$TE.common, TE[i])
      res$lower.common <- c(res$lower.common, lower[i])
      res$upper.common <- c(res$upper.common, upper[i])
      res$statistic.common <- c(res$statistic.common, statistic[i])
      res$pval.common <- c(res$pval.common, pval[i])
      res$text.common <- c(res$text.common, text[i])
      res$method <- c(res$method, "")
    }
    ##
    if (type[i] == "random") {
      res$random <- TRUE
      res$TE.random <- c(res$TE.random, TE[i])
      res$seTE.random <- c(res$seTE.random, NA)
      res$lower.random <- c(res$lower.random, lower[i])
      res$upper.random <- c(res$upper.random, upper[i])
      res$df.random <- c(res$df.random, NA)
      res$statistic.random <- c(res$statistic.random, statistic[i])
      res$pval.random <- c(res$pval.random, pval[i])
      res$text.random <- c(res$text.random, text[i])
      res$method.tau <- c(res$method.tau, "")
      res$method.random <- c(res$method.random, "")
      res$method.random.ci <- c(res$method.random.ci, "")
      res$adhoc.hakn.ci <- c(res$adhoc.hakn.ci, "")
    }
    ##
    if (type[i] == "prediction") {
      res$prediction <- TRUE
      res$lower.predict <- c(res$lower.predict, lower[i])
      res$upper.predict <- c(res$upper.predict, upper[i])
      res$text.predict <- c(res$text.predict, text[i])
      res$method.predict <- c(res$method.predict, "")
      res$adhoc.hakn.pi <- c(res$adhoc.hakn.pi, "")
    }
  }
  
  res
}
