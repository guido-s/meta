#' Calculate the number needed to treat (NNT)
#' 
#' @description
#' Calculate the number needed to treat (NNT) from estimated risk
#' difference, risk ratio, odds ratio, or hazard ratio, and a baseline
#' probability, i.e., control group event probability for binary outcomes or
#' survival probability for hazard ratios.
#' 
#' @aliases nnt nnt.default nnt.meta
#' 
#' @param x An object of class \code{meta}, or estimated treatment
#'   effect(s), i.e., risk difference(s), risk ratio(s), odds
#'   ratio(s), or hazard ratio(s).
#' @param p.c Baseline probability, i.e., control group event probability
#'   for binary outcomes or survival probability in control group for
#'   hazard ratios.
#' @param sm Summary measure.
#' @param lower Lower confidence interval limit.
#' @param upper Upper confidence interval limit.
#' @param common A logical indicating whether NNTs should be
#'   calculated based on common effect estimate.
#' @param random A logical indicating whether NNTs should be
#'   calculated based on random effects estimate.
#' @param small.values A character string specifying whether small
#'   treatment effects indicate a beneficial (\code{"desirable"}) or
#'   harmful (\code{"undesirable"}) effect, can be abbreviated.
#' @param transf A logical indicating whether treatment estimates and
#'   confidence limits are transformed or on the original scale.
#'   If \code{transf = TRUE}, inputs are expected to be log odds ratios instead
#'   of odds ratios for \code{sm = "OR"}, for example.
#' @param digits Minimal number of significant digits to print NNT and its
#'   confidence interval, see \code{print.default}.
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
#' risk ratio (RR), or odds ratio (OR) and a given baseline probability
#' (Higgins et al., 2023, section 15.4.4). It is also possible to
#' calculate NNTs from hazard ratios (HR) (Altman & Andersen,
#' 1999). Accordingly, NNTs can be calculated for meta-analyses
#' generated with \code{\link{metabin}} or \code{\link{metagen}} if
#' argument \code{sm} was equal to \code{"RD"}, \code{"RR"},
#' \code{"OR"}, or \code{"HR"}. It is also possible to provide only
#' estimated treatment effects and baseline probabilities (see Examples).
#'
#' The baseline probability can be specified using argument \code{p.c}. If
#' this argument is missing, the minimum, mean, and maximum of the
#' control event probabilities in the meta-analysis are used for
#' \code{\link{metabin}} and control event probabilities of 0.1, 0.2,
#' \dots, 0.9 are used for \code{\link{metagen}}. Note, the survival instead of
#' mortality probability must be provided for hazard ratios.
#'
#' Argument \code{small.values} can be used to specify whether small
#' treatment effects indicate a beneficial (\code{"desirable"}) or harmful
#' (\code{"undesirable"}) effect. For \code{small.values = "desirable"}, odds,
#' risk and hazard ratios below 1 and risk differences below 0 indicate that the
#' new treatment is beneficial. For \code{small.values = "undesirable"}, odds,
#' risk and hazard ratios above 1 and risk differences above 0 indicate that
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
#' for the NNT (Higgins et al., 2023, section 15.4.4).
#'
#' A peculiar problem arises if the confidence interval for the
#' summary measure includes the null effect (i.e., RR = 1, OR = 1, HR
#' = 1, or RD = 0). In this case the confidence interval for the NNT
#' contains both NNTB and NNTH values and it seemingly does not
#' include the estimated NNT.
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
#' Altman DG, Andersen PK (1999):
#' Calculating the number needed to treat for trials where the outcome
#' is time to an event.
#' \emph{British Medical Journal},
#' \bold{319}, 1492--95
#' 
#' Cook RJ, Sackett DL (1995):
#' The number needed to treat: a clinically useful measure of
#' treatment effect.
#' \emph{British Medical Journal},
#' \bold{310}, 452--54
#' 
#' Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch
#' VA (editors) (2023):
#' \emph{Cochrane Handbook for Systematic Reviews of Interventions
#'   Version 6.4 (updated August 2023)}.
#' Available from \url{https://www.cochrane.org/authors/handbooks-and-manuals/handbook}
#'
#' Laupacis A, Sackett DL, Roberts RS (1988):
#' An assessment of clinically useful measures of the consequences of
#' treatment.
#' \emph{New England Journal of Medicine},
#' \bold{318}, 1728--33
#' 
#' @examples
#' # Calculate NNTs for risk differences of -0.12 and -0.22
#' # (Cochrane Handbook, version 6.3, subsection 15.4.4.1)
#' nnt(c(-0.12, -0.22), sm = "RD")
#' 
#' # Calculate NNT for risk ratio of 0.92 and baseline risk of 0.3
#' # (Cochrane Handbook, version 6.3, subsection 15.4.4.2)
#' nnt(0.92, p.c = 0.3, sm = "RR")
#'
#' # Calculate NNT for odds ratio of 0.73 and baseline risk of 0.3
#' # (Cochrane Handbook, version 6.3, subsection 15.4.4.3)
#' nnt(0.73, p.c = 0.3, sm = "OR")
#'
#' # Calculate NNTs for Mantel-Haenszel odds ratio
#' data(Olkin1995)
#' m1 <-
#'   metabin(ev.exp, n.exp, ev.cont, n.cont, data = Olkin1995, random = FALSE)
#' nnt(m1)
#'
#' # Calculate NNTs from hazard ratio at two and four years (example from
#' # Altman & Andersen, 1999). Note, argument 'p.c' must provide survival
#' # probabilities instead of mortality rates for the control group.
#' nnt(0.72, lower = 0.55, upper = 0.92, sm = "HR", p.c = 1 - c(0.33, 0.49))
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
                     small.values = "desirable",
                     ...) {
  
  
  ##
  ## (1) Check for meta object and summary measure
  ##
  chkclass(x, "meta")
  x <- updateversion(x)
  ##
  if (!(x$sm %in% c("RD", "RR", "OR", "HR")))
    stop("Calculation of NNTs only possible for risk difference, ",
         "risk ratio, odds ratio, or hazard ratio as summary measure ",
         "(argument 'sm').")
  
  
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
      else {
        p.c <- seq(0.1, 0.9, by = 0.1)
        if (x$sm == "HR")
          p.c <- rev(p.c)
      }
    }
  }
  ##
  chknumeric(p.c, 0, 1)
  
  
  ##
  ## (3) Check other arguments
  ##
  small.values <- setsv(small.values)
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
  if (missing.common & x$k == 0 & x$k.all == 1 &
      all(is.na(x$TE.common)) & any(!is.na(x$TE))) {
    common <- TRUE
    x$TE.common <- x$TE
  }
  
  
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
  else if (x$sm == "HR") {
    if (common) {
      res$nnt.common$NNT <- logHR2nnt(res$TE.common, p.c)
      res$nnt.common$lower.NNT <- logHR2nnt(res$lower.common, p.c)
      res$nnt.common$upper.NNT <- logHR2nnt(res$upper.common, p.c)
    }
    ##
    if (random) {
      res$nnt.random$NNT <- logHR2nnt(res$TE.random, p.c)
      res$nnt.random$lower.NNT <- logHR2nnt(res$lower.random, p.c)
      res$nnt.random$upper.NNT <- logHR2nnt(res$upper.random, p.c)
    }
  }
  ##
  ## Switch direction and lower and upper limits
  ##
  if (small.values == "undesirable") {
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
  res$level.ma <- x$level.ma
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
                        small.values = "desirable",
                        transf = FALSE,
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
  sm <- setchar(sm, c("RD", "RR", "OR", "HR"))
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
  #
  small.values <- setsv(small.values)
  #
  chklogical(transf)
  
  
  ##
  ## (2) Calculate NNTs
  ##
  res <- data.frame(x = x, p.c = p.c,
                    NNT = NA, lower.NNT = NA, upper.NNT = NA)
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
  #
  else if (sm == "RR") {
    res$NNT <- logRR2nnt(if (!transf) log(res$x) else res$x, res$p.c)
    if (!missing.lower)
      res$lower.NNT <-
        logRR2nnt(if (!transf) log(res$lower.x) else res$lower.x, res$p.c)
    if (!missing.upper)
      res$upper.NNT <-
        logRR2nnt(if (!transf) log(res$upper.x) else res$upper.x, res$p.c)
  }
  #
  else if (sm == "OR") {
    res$NNT <- logOR2nnt(if (!transf) log(res$x) else res$x, res$p.c)
    if (!missing.lower)
      res$lower.NNT <-
        logOR2nnt(if (!transf) log(res$lower.x) else res$lower.x, res$p.c)
    if (!missing.upper)
      res$upper.NNT <-
        logOR2nnt(if (!transf) log(res$upper.x) else res$upper.x, res$p.c)
  }
  #
  else if (sm == "HR") {
    res$NNT <- logHR2nnt(if (!transf) log(res$x) else res$x, res$p.c)
    if (!missing.lower)
      res$lower.NNT <-
        logHR2nnt(if (!transf) log(res$lower.x) else res$lower.x, res$p.c)
    if (!missing.upper)
      res$upper.NNT <-
        logHR2nnt(if (!transf) log(res$upper.x) else res$upper.x, res$p.c)
  }
  ##
  names(res)[names(res) == "x"] <- sm
  names(res)[names(res) == "lower.x"] <- paste0("lower.", sm)
  names(res)[names(res) == "upper.x"] <- paste0("upper.", sm)
  ##
  ## Remove baseline probability for risk difference
  ##
  if (sm == "RD")
    res[["p.c"]] <- NULL
  else if (sm == "HR")
    names(res)[names(res) == "p.c"] <- "Surv.c"
  ##
  ## Switch direction and lower and upper limits
  ##
  if (small.values == "undesirable") {
    res$NNT <- -res$NNT
    ##
    tmp.lower <- res$lower.NNT
    if (!is.null(res$upper.NNT))
      res$lower.NNT <- -res$upper.NNT
    if (!is.null(tmp.lower))
      res$upper.NNT <- -tmp.lower
  }
  #
  attr(res, "small.values") <- small.values
  #
  class(res) <- c("nnt.default", class(res))
  
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
  
  
  ci.lab <- paste0(round(100 * replaceNULL(x$level.ma), 1), "%-CI")
  p.lab <- if (x$sm == "HR") "Surv.c" else "p.c"
  
  
  if (common) {
    cat(paste0("Number needed to treat (",
               tolower(gs("text.common")), "): \n\n"))
    ##
    nnt.common <- cbind(TE = x$TE.common, x$nnt.common)
    if (x$sm != "RD")
      nnt.common$TE <- exp(nnt.common$TE)
    ##
    ci.nnt <-
      sum(!is.na(x$lower.common)) > 0 & sum(!is.na(x$upper.common)) > 0
    ##
    sign <-
      all(x$lower.common[!is.na(x$lower.common)] > 0) |
      all(x$upper.common[!is.na(x$upper.common)] < 0)
    ##
    nnt.common$p.c <- formatN(round(nnt.common$p.c, digits.prop),
                              digits.prop, big.mark = big.mark)
    ##
    if (all(nnt.common$NNT < 0)) {
      lab.nnt <- "NNTH"
      lab.nnt2 <- "NNTB"
      ##
      nnt.common$NNT <- -nnt.common$NNT
      tmp.u <- -nnt.common$upper.NNT
      nnt.common$upper.NNT <- -nnt.common$lower.NNT
      nnt.common$lower.NNT <- tmp.u
    }
    else if (all(nnt.common$NNT > 0)) {
      lab.nnt <- "NNTB"
      lab.nnt2 <- "NNTH"
    }
    else {
      lab.nnt <- "NNT"
      lab.nnt2 <- "NNT2"
    }
    ##
    nnt.common$TE <-
      formatN(round(nnt.common$TE, digits),
              digits, big.mark = big.mark)
    ##
    nnt.common$NNT <- formatN(round(nnt.common$NNT, digits),
                              digits, big.mark = big.mark)
    ##
    if (sign) {
      nnt.common$lower.NNT <-
        formatCI(formatN(round(nnt.common$lower.NNT, digits),
                         digits, big.mark = big.mark),
                 formatN(round(nnt.common$upper.NNT, digits),
                         digits, big.mark = big.mark))
      nnt.common$upper.NNT <- NULL
      ##
      colnames(nnt.common) <- c(x$sm, p.lab, lab.nnt, ci.lab)
      ##
      if (!ci.nnt)
        nnt.common <- nnt.common[, -4]
      ##
      prmatrix(nnt.common, quote = FALSE, right = TRUE,
               rowlab = rep("", nrow(nnt.common)))
    }
    else {
      bracktype <- setchar(gs("CIbracket"), c("[", "(", "{", ""))
      if (bracktype == "[") {
        bracketLeft <- paste0("[", lab.nnt, " ")
        bracketRight <- "]"
      }
      else if (bracktype == "(") {
        bracketLeft <- paste0("(", lab.nnt, " ")
        bracketRight <- ")"
      }
      else if (bracktype == "{") {
        bracketLeft <- paste0("{", lab.nnt, " ")
        bracketRight <- "}"
      }
      else if (bracktype == "") {
        bracketLeft <- paste0("", lab.nnt, " ")
        bracketRight <- ""
      }
      ##
      nnt.common$lower.NNT <-
        formatCI(formatN(round(nnt.common$lower.NNT, digits),
                         digits, big.mark = big.mark),
                 formatN(round(-nnt.common$upper.NNT, digits),
                         digits, big.mark = big.mark),
                 bracket.left = bracketLeft,
                 bracket.right = bracketRight,
                 separator = paste0(" to Inf to ", lab.nnt2, " "))
      nnt.common$upper.NNT <- NULL
      ##
      colnames(nnt.common) <- c(x$sm, p.lab, lab.nnt, ci.lab)
      ##
      if (!ci.nnt)
        nnt.common <- nnt.common[, -4]
      ##
      prmatrix(nnt.common, quote = FALSE, right = TRUE,
               rowlab = rep("", nrow(nnt.common)))
    }
    ##
    if (random)
      cat("\n")
  }
  ##
  if (random) {
    cat(paste0("Number needed to treat (",
               tolower(gs("text.random")), "): \n\n"))
    ##
    nnt.random <- cbind(TE = x$TE.random, x$nnt.random)
    if (x$sm != "RD")
      nnt.random$TE <- exp(nnt.random$TE)
    ##
    ci.nnt <-
      sum(!is.na(x$lower.random)) > 0 & sum(!is.na(x$upper.random)) > 0
    ##
    sign <-
      all(x$lower.random[!is.na(x$lower.random)] > 0) |
      all(x$upper.random[!is.na(x$upper.random)] < 0)
    ##
    nnt.random$p.c <- formatN(round(nnt.random$p.c, digits.prop),
                              digits.prop, big.mark = big.mark)
    ##
    if (all(nnt.random$NNT < 0)) {
      lab.nnt <- "NNTH"
      lab.nnt2 <- "NNTB"
      ##
      nnt.random$NNT <- -nnt.random$NNT
      tmp.u <- -nnt.random$upper.NNT
      nnt.random$upper.NNT <- -nnt.random$lower.NNT
      nnt.random$lower.NNT <- tmp.u
    }
    else if (all(nnt.random$NNT > 0)) {
      lab.nnt <- "NNTB"
      lab.nnt2 <- "NNTH"
    }
    else {
      lab.nnt <- "NNT"
      lab.nnt2 <- "NNT2"
    }
    ##
    nnt.random$TE <-
      formatN(round(nnt.random$TE, digits),
              digits, big.mark = big.mark)
    ##
    nnt.random$NNT <- formatN(round(nnt.random$NNT, digits),
                              digits, big.mark = big.mark)
    ##
    if (sign) {
      nnt.random$lower.NNT <-
        formatCI(formatN(round(nnt.random$lower.NNT, digits),
                         digits, big.mark = big.mark),
                 formatN(round(nnt.random$upper.NNT, digits),
                         digits, big.mark = big.mark))
      nnt.random$upper.NNT <- NULL
      ##
      colnames(nnt.random) <- c(x$sm, p.lab, lab.nnt, ci.lab)
      ##
      if (!ci.nnt)
        nnt.random <- nnt.random[, -4]
      ##
      prmatrix(nnt.random, quote = FALSE, right = TRUE,
               rowlab = rep("", nrow(nnt.random)))
    }
    else {
      bracktype <- setchar(gs("CIbracket"), c("[", "(", "{", ""))
      if (bracktype == "[") {
        bracketLeft <- paste0("[", lab.nnt, " ")
        bracketRight <- "]"
      }
      else if (bracktype == "(") {
        bracketLeft <- paste0("(", lab.nnt, " ")
        bracketRight <- ")"
      }
      else if (bracktype == "{") {
        bracketLeft <- paste0("{", lab.nnt, " ")
        bracketRight <- "}"
      }
      else if (bracktype == "") {
        bracketLeft <- paste0("", lab.nnt, " ")
        bracketRight <- ""
      }
      ##
      nnt.random$lower.NNT <-
        formatCI(formatN(round(nnt.random$lower.NNT, digits),
                         digits, big.mark = big.mark),
                 formatN(round(-nnt.random$upper.NNT, digits),
                         digits, big.mark = big.mark),
                 bracket.left = bracketLeft,
                 bracket.right = bracketRight,
                 separator = paste0(" to Inf to ", lab.nnt2, " "))
      nnt.random$upper.NNT <- NULL
      ##
      colnames(nnt.random) <- c(x$sm, p.lab, lab.nnt, ci.lab)
      ##
      if (!ci.nnt)
        nnt.random <- nnt.random[, -4]
      ##
      prmatrix(nnt.random, quote = FALSE, right = TRUE,
               rowlab = rep("", nrow(nnt.random)))
    }
  }
  
  
  invisible(NULL)
}





#' @rdname nnt
#' @method print nnt.default
#' @export


print.nnt.default <- function(x,
                              digits = gs("digits"),
                              digits.prop = gs("digits.prop"),
                              big.mark = gs("big.mark"),
                              ...) {
  
  #
  # (1) Check for nnt.default object
  #
  chkclass(x, "nnt.default")
  
  #
  # (2) Check arguments
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.prop, min = 0, length = 1)
  
  #
  # (3) Round results
  #
  x$NNT <- formatN(round(x$NNT, digits), digits, big.mark = big.mark)
  if (isCol(x, "lower.NNT"))
    x$lower.NNT <-
      formatN(round(x$lower.NNT, digits), digits, big.mark = big.mark)
  if (isCol(x, "upper.NNT"))
    x$upper.NNT <-
      formatN(round(x$upper.NNT, digits), digits, big.mark = big.mark)
  #
  if (isCol(x, "p.c"))
    x$p.e <- formatN(round(x$p.c, digits.prop), digits.prop)
  if (isCol(x, "Surv.c"))
    x$p.e <- formatN(round(x$Surv.c, digits.prop), digits.prop)
  
  #
  # (4) Print results
  #
  class(x) <- class(x)[class(x) != "nnt.default"]
  print(x)
  
  invisible(NULL)
}





##
## Auxillary R functions
##
logRR2nnt <- function(x, p.c)
    1 / (p.c - exp(x) * p.c)
##
logOR2nnt <- function(x, p.c)
  1 / (p.c - exp(x) * p.c / (1 - p.c + exp(x) * p.c))
##
logHR2nnt <- function(x, Surv.c)
  1 / (Surv.c^exp(x) - Surv.c)
