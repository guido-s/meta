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
#' @param small.values A character string specifying whether small
#'   treatment effects indicate a beneficial (\code{"good"}) or
#'   harmful (\code{"bad"}) effect, can be abbreviated.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.prop Minimal number of significant digits for
#'   proportions, see \code{print.default}.
#' @param big.mark A character used as thousands separator.
#' @param \dots Additional arguments (to catch deprecated arguments).
#' 
#' @details
#' The number needed to treat (NNT) is the estimated number of
#' patients who need to be treated with a new treatment instead of a
#' standard for one additional patient to benefit (Laupacis et al.,
#' 1988; Cook & Sackett, 1995). This definition of the NNT implies
#' that the new treatment is more beneficial than the standard. If the
#' new treatment is indeed less beneficial than the standard, the NNT
#' gives the number of patients treated with the new treatment to
#' observe an additional harmful event. Accordingly, the abbreviations
#' NNTB and NNTH can be used to distinguish between beneficial and
#' harmful NNTs (Altman, 1998).
#'
#' NNTs can be easily computed from an estimated risk difference (RD),
#' risk ratio (RR), or odds ratio (OR) and a given baseline risk
#' (Higgins et al., 2022, section 15.4.4). Accordlingly, NNTs can be
#' calculated for meta-analyses generated with \code{\link{metabin}}
#' or \code{\link{metagen}} if argument \code{sm} was equal to
#' \code{"RD"}, \code{"RR"}, or \code{"OR"}. It is also possible to
#' provide only estimated treatment effects and baseline risks (see
#' Examples).
#'
#' The baseline risk can be specified using argument \code{p.c}. If
#' this argument is missing, the minimum, mean, and maximum of the
#' control event probabilities in the meta-analysis are used for
#' \code{\link{metabin}} and control event probabilities of 0.1, 0.2,
#' \dots, 0.9 are used for \code{\link{metagen}}.
#'
#' Argument \code{small.values} can be used to specify whether small
#' treatment effects indicate a beneficial (\code{"good"}) or harmful
#' (\code{"bad"}) effect. For \code{small.values = "small"}, odds and
#' risk ratios below 1 and risk differences below 0 indicate that the
#' new treatment is beneficial. For \code{small.values = "bad"}, odds
#' and risk ratios above 1 and risk differences above 0 indicate that
#' the new treatment is beneficial.
#'
#' \subsection{Interpretation of (positive and negative) NNTs}{
#' A positive value for the estimated NNT indicates that the new
#' treatment is beneficial, i.e., the NNT is actually an NNTB. On the
#' other hand, a negative value for the estimated NNT indicates that
#' the new treatment is harmful, i.e., the NNT is actually an NNTH.
#' 
#' The minimal value for the NNTB is 1. In this extreme case the new
#' treatment is 100\% effective and the standard treatment is 0\%
#' effective, i.e., only one patient has to be treated with the new
#' treatment for one additional patient to benefit. The NNTB increases
#' with decreasing difference between the two risks. If both risks are
#' equal, the NNTB is infinite.
#'
#' The other extreme is an NNT of -1 if the new treatment is 0\%
#' effective and the standard is 100\% effective. Here, one additional
#' harmful event is observed for each patient treated with the new
#' treatment. The NNT approaches minus infinity if the difference
#' between the two risks decreases to 0. Finally, an NNT of -1
#' translates to an NNTH of 1 with possible values from 1 to infinity.
#' }
#' 
#' \subsection{Confidence interval for the NNT}{
#' Confidence limits for an NNT are derived from the lower and upper
#' confidence limits of the summary measure using the same formulae as
#' for the NNT (Higgins et al., 2022, section 15.4.4).
#'
#' A peculiar problem arises if the confidence interval for the
#' summary measure includes the null effect (i.e., RR = 1, OR = 1, or
#' RD = 0). In this case the confidence interval for the NNT contains
#' both NNTB and NNTH values and it seemingly does not include the
#' estimated NNT.
#'
#' As described above, a positive NNT value corresponds to an NNTB and
#' the absolute value of a negative NNT is equal to an
#' NNTH. Accordingly, a confidence interval for the NNT from 20 to -5
#' translates to NNTB values between 20 and infinity and NNTH values
#' between 5 and infinity (Altman, 1998).
#' }
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{metagen}}
#' 
#' @references
#' Altman DG (1998):
#' Confidence intervals for the number needed to treat.
#' \emph{British Medical Journal},
#' \bold{317}, 1309--12
#'
#' Cook RJ, Sackett DL (1995):
#' The Number Needed to Treat: A Clinically Useful Measure of
#' Treatment Effect.
#' \emph{British Medical Journal},
#' \bold{310}, 452--54
#' 
#' Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch
#' VA (editors) (2022):
#' \emph{Cochrane Handbook for Systematic Reviews of Interventions
#'   Version 6.3 (updated February 2022)}.
#' Available from www.training.cochrane.org/handbook
#'
#' Laupacis A, Sackett DL, Roberts RS (1988):
#' An Assessment of Clinically Useful Measures of the Consequences of
#' Treatment.
#' \emph{New England Journal of Medicine},
#' \bold{318}, 1728--33
#' 
#' @examples
#' # Calculate NNT for RD = -0.12 and -0.22
#' # (Cochrane Handbook, version 6.3, subsection 15.4.4.1)
#' nnt(c(-0.12, -0.22), sm = "RD")
#' 
#' # Calculate NNT for RR = 0.92 and baseline risk p.c = 0.3
#' # (Cochrane Handbook, version 6.3, subsection 15.4.4.2)
#' nnt(0.92, p.c = 0.3, sm = "RR")
#'
#' # Calculate NNT for OR = 0.73 and baseline risk p.c = 0.3
#' # (Cochrane Handbook, version 6.3, subsection 15.4.4.3)
#' nnt(0.73, p.c = 0.3, sm = "OR")
#'
#' # Use Mantel-Haenszel odds ratio to calculate NNTs
#' data(Olkin1995)
#' m1 <- metabin(ev.exp, n.exp, ev.cont, n.cont, data = Olkin1995,
#'               random = FALSE)
#' nnt(m1)
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
                     random = x$random,
                     small.values = "good",
                     ...) {
  
  
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
      if (!is.null(x$event.c) & !is.null(x$n.c)) {
        p.c <-
          unique(as.vector(
            summary(x$event.c / x$n.c,
                    na.rm = TRUE)[c("Min.", "Mean", "Max.")]))
        if (!all(p.c == 0))
          p.c <- p.c[p.c != 0]
      }
      else
        p.c <- seq(0.1, 0.9, by = 0.1)
    }
  }
  ##
  chknumeric(p.c, 0, 1)
  
  
  ##
  ## (3) Check other arguments
  ##
  small.values <- setchar(small.values, c("good", "bad"))
  ##
  args  <- list(...)
  missing.common <- missing(common)
  common <- deprecated(common, missing.common, args, "fixed", FALSE)
  chklogical(common)
  chklogical(random)
  ##
  if (missing.common & !common & !random & x$k == 1)
    common <- TRUE
  
  
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
      res$nnt.common$NNT <- -1 / res$TE.common
      res$nnt.common$lower.NNT <- -1 / res$lower.common
      res$nnt.common$upper.NNT <- -1 / res$upper.common
    }
    ##
    if (random) {
      res$nnt.random$NNT <- -1 / res$TE.random
      res$nnt.random$lower.NNT <- -1 / res$lower.random
      res$nnt.random$upper.NNT <- -1 / res$upper.random
    }
  }
  ##
  else if (x$sm == "RR") {
    if (common) {
      res$nnt.common$NNT <- logRR2nnt(res$TE.common, p.c)
      res$nnt.common$lower.NNT <- logRR2nnt(res$lower.common, p.c)
      res$nnt.common$upper.NNT <- logRR2nnt(res$upper.common, p.c)
    }
    ##
    if (random) {
      res$nnt.random$NNT <- logRR2nnt(res$TE.random, p.c)
      res$nnt.random$lower.NNT <- logRR2nnt(res$lower.random, p.c)
      res$nnt.random$upper.NNT <- logRR2nnt(res$upper.random, p.c)
    }
  }
  ##
  else if (x$sm == "OR") {
    if (common) {
      res$nnt.common$NNT <- logOR2nnt(res$TE.common, p.c)
      res$nnt.common$lower.NNT <- logOR2nnt(res$lower.common, p.c)
      res$nnt.common$upper.NNT <- logOR2nnt(res$upper.common, p.c)
    }
    ##
    if (random) {
      res$nnt.random$NNT <- logOR2nnt(res$TE.random, p.c)
      res$nnt.random$lower.NNT <- logOR2nnt(res$lower.random, p.c)
      res$nnt.random$upper.NNT <- logOR2nnt(res$upper.random, p.c)
    }
  }
  ##
  ## Switch direction and lower and upper limits
  ##
  if (small.values == "bad") {
    if (common) {
      res$nnt.common$NNT <- -res$nnt.common$NNT
      tmp.lower <- res$nnt.common$lower.NNT
      res$nnt.common$lower.NNT <- -res$nnt.common$upper.NNT
      res$nnt.common$upper.NNT <- -tmp.lower
    }
    ##
    if (random) {
      res$nnt.random$NNT <- -res$nnt.random$NNT
      tmp.lower <- res$nnt.random$lower.NNT
      res$nnt.random$lower.NNT <- -res$nnt.random$upper.NNT
      res$nnt.random$upper.NNT <- -tmp.lower
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


nnt.default <- function(x, p.c, sm, lower, upper,
                        small.values = "good",
                        ...) {
  
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
  small.values <- setchar(small.values, c("good", "bad"))
  
  
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
    res$NNT <- logRR2nnt(log(res$x), res$p.c)
    if (!missing.lower)
      res$lower.NNT <- logRR2nnt(log(res$lower.x), res$p.c)
    if (!missing.upper)
      res$upper.NNT <- logRR2nnt(log(res$upper.x), res$p.c)
  }
  ##
  else if (sm == "OR") {
    res$NNT <- logOR2nnt(log(res$x), res$p.c)
    if (!missing.lower)
      res$lower.NNT <- logOR2nnt(log(res$lower.x), res$p.c)
    if (!missing.upper)
      res$upper.NNT <- logOR2nnt(log(res$upper.x), res$p.c)
  }
  ##
  names(res)[names(res) == "x"] <- sm
  names(res)[names(res) == "lower.x"] <- paste0("lower.", sm)
  names(res)[names(res) == "upper.x"] <- paste0("upper.", sm)
  ##
  ## Switch direction and lower and upper limits
  ##
  if (small.values == "bad") {
    res$NNT <- -res$NNT
    ##
    tmp.lower <- res$lower.NNT
    if (!is.null(res$upper.NNT))
      res$lower.NNT <- -res$upper.NNT
    if (!is.null(tmp.lower))
      res$upper.NNT <- -tmp.lower
  }
  
  
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
logRR2nnt <- function(x, p.c)
    1 / (p.c * (1 - exp(x)))
##
logOR2nnt <- function(x, p.c)
  1 / (p.c - exp(x) * p.c / (1 - p.c + exp(x) * p.c))
