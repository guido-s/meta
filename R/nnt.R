#' Calculate the number needed to treat (NNT)
#' 
#' @description
#' Calculate the number needed to treat (NNT) from estimated risk
#' difference, risk ratio, or odds ratio, and a baseline risk.
#' 
#' @aliases nnt nnt.default nnt.meta
#' 
#' @param x An object of class \code{meta}, or estimated treatment
#'   effect, i.e., risk difference(s), risk ratio(s), or odds
#'   ratio(s).
#' @param p.c Baseline risk (control group event probability).
#' @param sm Summary measure.
#' @param lower Lower confidence interval limit.
#' @param upper Upper confidence interval limit.
#' @param common A logical indicating whether NNTs should be
#'   calculated based on common effect estimate.
#' @param random A logical indicating whether NNTs should be
#'   calculated based on random effects estimate.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.prop Minimal number of significant digits for
#'   proportions, see \code{print.default}.
#' @param big.mark A character used as thousands separator.
#' @param \dots Additional arguments (to catch deprecated arguments).
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
#' data(Olkin1995)
#' m1 <- metabin(ev.exp, n.exp, ev.cont, n.cont, data = Olkin1995,
#'               random = FALSE)
#' nnt(m1, random = TRUE)
#' 
#' @rdname nnt
#' @export nnt


nnt <- function(x, ...) 
  UseMethod("nnt")





#' @rdname nnt
#' @method nnt meta
#' @export


nnt.meta <- function(x, p.c,
                     common = x$common,
                     random = x$random, ...) {
  
  
  ##
  ## (1) Check for meta object and summary measure
  ##
  chkclass(x, "meta")
  x <- updateversion(x)
  ##
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
  args  <- list(...)
  common <- deprecated(common, missing(common), args, "fixed", FALSE)
  chklogical(common)
  chklogical(random)
  
  
  ##
  ## (4) Calculate NNTs
  ##
  res <- list()
  ##
  if (common) {
    res$nnt.common <- data.frame(p.c = p.c, NNT = NA,
                                lower.NNT = NA, upper.NNT = NA)
    ##
    res$TE.common    <- x$TE.common
    res$lower.common <- x$lower.common
    res$upper.common <- x$upper.common
  }
  ##
  if (random) {
    res$nnt.random <- data.frame(p.c = p.c, NNT = NA,
                                 lower.NNT = NA, upper.NNT = NA)
    ##
    res$TE.random    <- x$TE.random
    res$lower.random <- x$lower.random
    res$upper.random <- x$upper.random
  }
  ##
  if (x$sm == "RD") {
    if (common) {
      res$nnt.common$NNT       <- -1 / res$TE.common
      res$nnt.common$lower.NNT <- -1 / res$lower.common
      res$nnt.common$upper.NNT <- -1 / res$upper.common
    }
    ##
    if (random) {
      res$nnt.random$NNT       <- -1 / res$TE.random
      res$nnt.random$lower.NNT <- -1 / res$lower.random
      res$nnt.random$upper.NNT <- -1 / res$upper.random
    }
  }
  ##
  else if (x$sm == "RR") {
    if (common) {
      res$nnt.common$NNT       <- logRR2NNT(res$TE.common, p.c)
      res$nnt.common$lower.NNT <- logRR2NNT(res$lower.common, p.c)
      res$nnt.common$upper.NNT <- logRR2NNT(res$upper.common, p.c)
    }
    ##
    if (random) {
      res$nnt.random$NNT       <- logRR2NNT(res$TE.random, p.c)
      res$nnt.random$lower.NNT <- logRR2NNT(res$lower.random, p.c)
      res$nnt.random$upper.NNT <- logRR2NNT(res$upper.random, p.c)
    }
  }
  ##
  else if (x$sm == "OR") {
    if (common) {
      res$nnt.common$NNT       <- logOR2NNT(res$TE.common, p.c)
      res$nnt.common$lower.NNT <- logOR2NNT(res$lower.common, p.c)
      res$nnt.common$upper.NNT <- logOR2NNT(res$upper.common, p.c)
    }
    ##
    if (random) {
      res$nnt.random$NNT       <- logOR2NNT(res$TE.random, p.c)
      res$nnt.random$lower.NNT <- logOR2NNT(res$lower.random, p.c)
      res$nnt.random$upper.NNT <- logOR2NNT(res$upper.random, p.c)
    }
  }
  ##
  res$sm <- x$sm
  res$common <- common
  res$random <- random
  ##
  res$call <- match.call()
  res$version <- packageDescription("meta")$Version
  ##
  res$fixed <- common
  ##
  class(res) <- "nnt.meta"
  
  
  res
}





#' @rdname nnt
#' @method nnt default
#' @export


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


print.nnt.meta <- function(x,
                           common = x$common,
                           random = x$random,
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
  args  <- list(...)
  common <- deprecated(common, missing(common), args, "fixed", FALSE)
  if (is.null(common) & !is.null(x$fixed))
    common <- x$fixed
  chklogical(common)
  chklogical(random)
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.prop, min = 0, length = 1)


  if (common) {
    cat(paste0(gs("text.common"), ": \n\n"))
    x$nnt.common$p.c <- formatN(round(x$nnt.common$p.c, digits.prop),
                               digits.prop, big.mark = big.mark)
    ##
    x$nnt.common$NNT <- formatN(round(x$nnt.common$NNT, digits),
                               digits, big.mark = big.mark)
    ##
    x$nnt.common$lower.NNT <- formatN(round(x$nnt.common$lower.NNT, digits),
                                     digits, big.mark = big.mark)
    ##
    x$nnt.common$upper.NNT <- formatN(round(x$nnt.common$upper.NNT, digits),
                                     digits, big.mark = big.mark)
    ##
    prmatrix(x$nnt.common, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(x$nnt.common)))
    if (random)
      cat("\n")
  }
  ##
  if (random) {
    cat(paste0(gs("text.random"), ": \n\n"))
    x$nnt.random$p.c <- formatN(round(x$nnt.random$p.c, digits.prop),
                                digits.prop, big.mark = big.mark)
    ##
    x$nnt.random$NNT <- formatN(round(x$nnt.random$NNT, digits),
                                digits, big.mark = big.mark)
    ##
    x$nnt.random$lower.NNT <- formatN(round(x$nnt.random$lower.NNT, digits),
                                      digits, big.mark = big.mark)
    ##
    x$nnt.random$upper.NNT <- formatN(round(x$nnt.random$upper.NNT, digits),
                                      digits, big.mark = big.mark)
    ##
    prmatrix(x$nnt.random, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(x$nnt.common)))
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
