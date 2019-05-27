#' Calculate the number needed to treat (NNT)
#' 
#' @description
#' Calculate the number needed to treat (NNT) from estimated risk
#' difference, risk ratio, or odds ratio, and a baseline risk.
#' 
#' @aliases nnt nnt.default nnt.meta
#' 
#' @param x An object of class \code{meta}, or estimated treatment
#'   effect, i.e., risk difference(s), risk ratio(s), or odds ratio(s).
#' @param p.c Baseline risk (control group event probability).
#' @param sm Summary measure.
#' @param lower Lower confidence interval limit.
#' @param upper Upper confidence interval limit.
#' @param comb.fixed A logical indicating whether NNTs should be
#'   calculated based on fixed effect estimate.
#' @param comb.random A logical indicating whether NNTs should be
#'   calculated based on random effects estimate.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.prop Minimal number of significant digits for
#'   proportions, see \code{print.default}.
#' @param big.mark A character used as thousands separator.
#' @param \dots Additional arguments (ignored at the moment).
#' 
#' @details
#' The number needed to treat (NNT) can be easily computed from an
#' estimated risk difference (RD), risk ratio (RR), or odds ratio (OR)
#' and a given baseline risk (Higgins & Green, 2011, section 12.5).
#'
#' Accordlingly, this function can be used to calculate NNTs for
#' meta-analyses generated with \code{\link{metabin}} or
#' \code{\link{metagen}} if argument \code{sm} was equal to
#' \code{"RD"}, \code{"RR"}, or \code{"OR"}. It is also possible to
#' directly provide estimated treatment effects without conducting a
#' meta-analysis (see Examples).
#'
#' The baseline risk can be specified using argument \code{p.c}. If
#' this argument is missing, the minimum, mean, and maximum of the
#' control event probabilities in the meta-analysis are used for
#' \code{\link{metabin}}; otherwise the control event probabilities
#' 0.1, 0.2, \dots, 0.9 are used.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{metagen}}
#' 
#' @references
#' Higgins, J.P.T and S. Green (2011):
#' \emph{Cochrane Handbook for Systematic Reviews of Interventions
#'   Version 5.1.0 [Updated March 2011]}.
#' The Cochrane Library: http://www.cochrane-handbook.org
#' 
#' @examples
#' # Calculate NNT for RD = -0.21
#' # (Cochrane Handbook, version 5.1, subsection 12.5.4.1)
#' nnt(-0.21, sm = "RD")
#' 
#' # Calculate NNT for RR = 0.92 and baseline risk p.c = 0.3
#' # (Cochrane Handbook, version 5.1, subsection 12.5.4.2)
#' nnt(0.92, p.c = 0.3, sm = "RR")
#'
#' # Calculate NNT for OR = 0.73 and baseline risk p.c = 0.3
#' # (Cochrane Handbook, version 5.1, subsection 12.5.4.3)
#' nnt(0.73, p.c = 0.3, sm = "OR")
#'
#' # Use Mantel-Haenszel odds ratio to calculate NNTs
#' data(Olkin95)
#' m1 <- metabin(event.e, n.e, event.c, n.c, data = Olkin95,
#'               comb.random = FALSE)
#' nnt(m1, comb.random = TRUE)
#' 
#' @rdname nnt
#' @export nnt


nnt <- function(x, ...) 
  UseMethod("nnt")





#' @rdname nnt
#' @method nnt meta
#' @export
#' @export nnt.meta


nnt.meta <- function(x, p.c,
                     comb.fixed = x$comb.fixed,
                     comb.random = x$comb.random, ...) {
  
  
  ##
  ## (1) Check for meta object and summary measure
  ##
  chkclass(x, "meta")
  if (!(x$sm %in% c("RD", "RR", "OR")))
    stop("Calculation of NNTs only possible for risk difference, ",
         "risk ratio, or odds ratio as summary measure (argument 'sm').")


  ##
  ## (2) Check / set baseline risks
  ##
  if (missing(p.c)) {
    if (x$sm == "RD")
      p.c <- NA
    else {
      if (!is.null(x$event.c) & !is.null(x$n.c))
        p.c <- as.vector(summary(x$event.c / x$n.c,
                                 na.rm = TRUE)[c("Min.", "Mean", "Max.")])
      else
        p.c <- 1:9 / 10
    }
  }
  ##
  chknumeric(p.c, 0, 1)
  
  
  ##
  ## (3) Check other arguments
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  
  
  ##
  ## (4) Calculate NNTs
  ##
  res <- list()
  ##
  if (comb.fixed) {
    res$fixed <- data.frame(p.c = p.c, NNT = NA,
                            lower.NNT = NA, upper.NNT = NA)
    ##
    res$TE.fixed    <- x$TE.fixed
    res$lower.fixed <- x$lower.fixed
    res$upper.fixed <- x$upper.fixed
  }
  ##
  if (comb.random) {
    res$random <- data.frame(p.c = p.c, NNT = NA,
                             lower.NNT = NA, upper.NNT = NA)
    ##
    res$TE.random    <- x$TE.random
    res$lower.random <- x$lower.random
    res$upper.random <- x$upper.random
  }
  ##
  if (x$sm == "RD") {
    if (comb.fixed) {
      res$fixed$NNT       <- -1 / res$TE.fixed
      res$fixed$lower.NNT <- -1 / res$lower.fixed
      res$fixed$upper.NNT <- -1 / res$upper.fixed
    }
    ##
    if (comb.random) {
      res$random$NNT       <- -1 / res$TE.random
      res$random$lower.NNT <- -1 / res$lower.random
      res$random$upper.NNT <- -1 / res$upper.random
    }
  }
  ##
  else if (x$sm == "RR") {
    if (comb.fixed) {
      res$fixed$NNT       <- logRR2NNT(res$TE.fixed, p.c)
      res$fixed$lower.NNT <- logRR2NNT(res$lower.fixed, p.c)
      res$fixed$upper.NNT <- logRR2NNT(res$upper.fixed, p.c)
    }
    ##
    if (comb.random) {
      res$random$NNT       <- logRR2NNT(res$TE.random, p.c)
      res$random$lower.NNT <- logRR2NNT(res$lower.random, p.c)
      res$random$upper.NNT <- logRR2NNT(res$upper.random, p.c)
    }
  }
  ##
  else if (x$sm == "OR") {
    if (comb.fixed) {
      res$fixed$NNT       <- logOR2NNT(res$TE.fixed, p.c)
      res$fixed$lower.NNT <- logOR2NNT(res$lower.fixed, p.c)
      res$fixed$upper.NNT <- logOR2NNT(res$upper.fixed, p.c)
    }
    ##
    if (comb.random) {
      res$random$NNT       <- logOR2NNT(res$TE.random, p.c)
      res$random$lower.NNT <- logOR2NNT(res$lower.random, p.c)
      res$random$upper.NNT <- logOR2NNT(res$upper.random, p.c)
    }
  }
  ##
  res$sm <- x$sm
  res$comb.fixed <- comb.fixed
  res$comb.random <- comb.random
  ##
  res$call <- match.call()
  res$version <- packageDescription("meta")$Version
  ##
  class(res) <- "nnt.meta"
  
  
  res
}





#' @rdname nnt
#' @method nnt default
#' @export
#' @export nnt.default


nnt.default <- function(x, p.c, sm, lower, upper, ...) {
  
  ##
  ## (1) Check arguments
  ##
  if (missing(x))
    stop("Argument 'x' is mandatory.")
  ##
  if (missing(sm))
    stop("Argument 'sm' is mandatory.")
  ##
  sm <- setchar(sm, c("RD", "RR", "OR"))
  ##
  if (missing(p.c)) {
    if (sm == "RD")
      p.c <- NA
    else
      stop("Argument 'p.c' is mandatory.")
  }
  ##
  chknumeric(p.c, 0, 1)
  ##
  missing.lower <- missing(lower)
  if (!missing.lower)
      chklength(lower, length(x), "nnt.default")
  ##
  missing.upper <- missing(upper)
  if (!missing.upper)
    chklength(upper, length(x), "nnt.default")
  ##
  ## Check length of arguments 'x' and 'p.c'
  ##
  if (length(x) != length(p.c)) {
    if (length(x) > length(p.c))
      bad <- length(x) %% length(p.c) != 0
    else
      bad <- length(p.c) %% length(x) != 0
    ##
    if (bad)
      stop("Lengths of arguments 'x' and 'p.c' do not match.",
             call. = FALSE)
  }
  
  
  ##
  ## (2) Calculate NNTs
  ##
  res <- data.frame(p.c = p.c,
                    NNT = NA, lower.NNT = NA, upper.NNT = NA,
                    x = x)
  ##
  if (missing.lower)
    res$lower.NNT <- NULL
  if (missing.upper)
    res$upper.NNT <- NULL
  ##
  if (!missing.lower)
    res$lower.x <- lower
  if (!missing.upper)
    res$upper.x <- upper
  ##
  if (sm == "RD") {
    res$NNT <- -1 / res$x
    if (!missing.lower)
      res$lower.NNT <- -1 / res$lower.x
    if (!missing.upper)
      res$upper.NNT <- -1 / res$upper.x
  }
  ##
  else if (sm == "RR") {
    res$NNT <- logRR2NNT(log(res$x), res$p.c)
    if (!missing.lower)
      res$lower.NNT <- logRR2NNT(log(res$lower.x), res$p.c)
    if (!missing.upper)
      res$upper.NNT <- logRR2NNT(log(res$upper.x), res$p.c)
  }
  ##
  else if (sm == "OR") {
    res$NNT <- logOR2NNT(log(res$x), res$p.c)
    if (!missing.lower)
      res$lower.NNT <- logOR2NNT(log(res$lower.x), res$p.c)
    if (!missing.upper)
      res$upper.NNT <- logOR2NNT(log(res$upper.x), res$p.c)
  }
  ##
  names(res)[names(res) == "x"] <- sm
  names(res)[names(res) == "lower.x"] <- paste0("lower.", sm)
  names(res)[names(res) == "upper.x"] <- paste0("upper.", sm)

  
  res
}





#' @rdname nnt
#' @method print nnt.meta
#' @export
#' @export print.nnt.meta


print.nnt.meta <- function(x,
                           comb.fixed = x$comb.fixed,
                           comb.random = x$comb.random,
                           digits = gs("digits"),
                           digits.prop = gs("digits.prop"),
                       big.mark = gs("big.mark"),
                           ...) {
  
  
  ##
  ## (1) Check for nnt.meta object
  ##
  chkclass(x, "nnt.meta")
  
  
  ##
  ## (2) Check arguments
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  ##
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.prop, min = 0, single = TRUE)


  if (comb.fixed) {
    cat("Fixed effect model: \n\n")
    x$fixed$p.c <- formatN(round(x$fixed$p.c, digits.prop),
                           digits.prop, big.mark = big.mark)
    ##
    x$fixed$NNT <- formatN(round(x$fixed$NNT, digits),
                           digits, big.mark = big.mark)
    ##
    x$fixed$lower.NNT <- formatN(round(x$fixed$lower.NNT, digits),
                                 digits, big.mark = big.mark)
    ##
    x$fixed$upper.NNT <- formatN(round(x$fixed$upper.NNT, digits),
                                 digits, big.mark = big.mark)
    ##
    prmatrix(x$fixed, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(x$fixed)))
    if (comb.random)
      cat("\n")
  }
  ##
  if (comb.random) {
    cat("Random effects model: \n\n")
    x$random$p.c <- formatN(round(x$random$p.c, digits.prop),
                           digits.prop, big.mark = big.mark)
    ##
    x$random$NNT <- formatN(round(x$random$NNT, digits),
                           digits, big.mark = big.mark)
    ##
    x$random$lower.NNT <- formatN(round(x$random$lower.NNT, digits),
                                 digits, big.mark = big.mark)
    ##
    x$random$upper.NNT <- formatN(round(x$random$upper.NNT, digits),
                                 digits, big.mark = big.mark)
    ##
    prmatrix(x$random, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(x$fixed)))
  }
  
  
  invisible(NULL)
}





##
## Auxillary R functions
##
logRR2NNT <- function(logRR, p.c)
  1 / (p.c * (1 - exp(logRR)))
##
logOR2NNT <- function(logOR, p.c)
  1 / (p.c - exp(logOR) * p.c / (1 - p.c + exp(logOR) * p.c))
