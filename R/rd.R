#' Calculate the risk difference from meta-analysis results
#' 
#' @description
#' Calculate the risk difference from the risk ratio, odds ratio, or hazard
#' ratio, and a baseline probability, i.e., control group event probability for
#' binary outcomes or survival probability for hazard ratios. Results are
#' reported as absolute risk reduction/increase or absolute benefit
#' increase/reduction.
#' 
#' @aliases rd rd.default rd.meta
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
#' @param level The level used to calculate confidence intervals.
#' @param common A logical indicating whether ARRs / ABIs should be
#'   calculated based on common effect estimate.
#' @param random A logical indicating whether ARRs / ABIs should be
#'   calculated based on random effects estimate.
#' @param small.values A character string specifying whether small
#'   treatment effects indicate a beneficial (\code{"desirable"}) or
#'   harmful (\code{"undesirable"}) effect, can be abbreviated.
#' @param transf A logical indicating whether treatment estimates and
#'   confidence limits are transformed or on the original scale.
#'   If \code{transf = TRUE}, inputs are expected to be log odds ratios instead
#'   of odds ratios for \code{sm = "OR"}, for example.
#' @param pscale A numeric defining a scaling factor for printing of
#'   absolute risk reduction or absolute benefit increase.
#' @param digits Minimal number of significant digits to print ARR / ABI and its
#'   confidence interval, see \code{print.default}.
#' @param digits.prop Minimal number of significant digits for
#'   proportions, see \code{print.default}.
#' @param big.mark A character used as thousands separator.
#' @param details A logical specifying whether details should be printed.
#' @param \dots Additional arguments (ignored)
#' 
#' @details
#' Calculate the risk difference from the risk ratio, odds ratio, or hazard
#' ratio, and a baseline probability, i.e., control group event probability for
#' binary outcomes or survival probability for hazard ratios. Report results as
#' absolute risk reduction (ARR), absolute risk increase (ARI), absolute
#' benefit increase (ABI), or absolute benefit reduction (ABR).
#'
#' Absolute measures can be easily computed from an estimated risk
#' difference (RD), or a risk ratio (RR) or odds ratio (OR) and a given
#' baseline probability. It is also possible to calculate absolute
#' measures from hazard ratios (HR) (Altman & Andersen, 1999).
#' Accordingly, these measures can be calculated for meta-analyses
#' generated with \code{\link{metabin}} or \code{\link{metagen}}
#' if argument \code{sm} was equal to \code{"RD"}, \code{"RR"},
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
#' (\code{"undesirable"}) effect. For \code{small.values = "desirable"}, OR < 1,
#' RR < 1, HR < 1, or RD > 0 indicate that the new treatment is beneficial.
#' In this case, calculated risk differences are expressed as absolute risk
#' reductions (ARR). If OR > 1, RR > 1, HR > 1, or RD > 0, calculated risk
#' differences are expressed as absolute risk increases (ARI). For
#' \code{small.values = "undesirable"}, odds, OR > 1, RR > 1, HR > 1, or RD > 0
#' indicate that the new treatment is beneficial. In this case, calculated risk
#' differences are expressed as absolute benefit increases (ABI). If OR < 1,
#' RR < 1, HR < 1, or RD < 0, calculated risk differences are expressed as
#' absolute benefit reductions (ABR).
#' 
#' Argument \code{pscale} can be used to rescale ARRs, ..., ABRs, e.g.,
#' \code{pscale = 1000} means that these quantities are expressed as events
#' per 1000 observations. This is useful in situations with (very) low
#' event probabilities.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{metagen}}
#' 
#' @references
#' Altman DG, Andersen PK (1999):
#' Calculating the number needed to treat for trials where the outcome
#' is time to an event.
#' \emph{British Medical Journal},
#' \bold{319}, 1492--95
#' 
#' @examples
#' # Calculate ARRs for risk differences of -0.12 and -0.22
#' rd(c(-0.12, -0.22), sm = "RD")
#' 
#' # Calculate ARR for risk ratio of 0.92 and baseline risk of 0.3
#' rd(0.92, p.c = 0.3, sm = "RR")
#'
#' # Calculate ARR for odds ratio of 0.73 and baseline risk of 0.3
#' rd(0.73, p.c = 0.3, sm = "OR")
#'
#' # Calculate ARRs for Mantel-Haenszel odds ratio
#' data(Olkin1995)
#' ma <-
#'   metabin(ev.exp, n.exp, ev.cont, n.cont, data = Olkin1995, random = FALSE)
#' rd(ma)
#'
#' # Calculate ARRs from hazard ratio at two and four years (example from
#' # Altman & Andersen, 1999). Note, argument 'p.c' must provide survival
#' # probabilities instead of mortality rates for the control group.
#' rd(0.72, lower = 0.55, upper = 0.92, sm = "HR", p.c = 1 - c(0.33, 0.49))
#' 
#' @rdname rd
#' @export rd

rd <- function(x, ...) 
  UseMethod("rd")





#' @rdname rd
#' @method rd meta
#' @export

rd.meta <- function(x, p.c,
                    common = x$common,
                    random = x$random,
                    small.values = "desirable",
                    pscale = 1,
                    ...) {
  
  
  #
  # (1) Check for meta object and summary measure
  #
  chkclass(x, "meta")
  x <- updateversion(x)
  #
  if (!(x$sm %in% c("RD", "RR", "OR", "HR")))
    stop("Calculation of ARRs /ABIs only possible for risk difference, ",
         "risk ratio, odds ratio, or hazard ratio as summary measure ",
         "(argument 'sm').")
  #
  chknumeric(pscale, length = 1)
  
  
  #
  # (2) Check / set baseline risks
  #
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
  #
  chknumeric(p.c, 0, 1)
  
  
  #
  # (3) Check other arguments
  #
  small.values <- setsv(small.values)
  #
  missing.common <- missing(common)
  chklogical(common)
  chklogical(random)
  #
  if (missing.common & !common & !random & x$k == 1)
    common <- TRUE
  #
  if (missing.common & x$k == 0 & x$k.all == 1 &
      all(is.na(x$TE.common)) & any(!is.na(x$TE))) {
    common <- TRUE
    x$TE.common <- x$TE
  }
  
  
  #
  # (4) Calculate RDs
  #
  res <- list()
  #
  if (common)
    res$common <- data.frame(
      TE = x$TE.common, lower = x$lower.common, upper = x$upper.common,
      p.c = p.c, RD = NA, lower.RD = NA, upper.RD = NA)
  #
  if (random)
    res$random <- data.frame(
      TE = x$TE.random, lower = x$lower.random, upper = x$upper.random,
      p.c = p.c, RD = NA, lower.RD = NA, upper.RD = NA)
  #
  if (x$sm == "RD") {
    if (common) {
      res$common$RD <- x$TE.common
      res$common$lower.RD <- x$lower.common
      res$common$upper.RD <- x$upper.common
    }
    #
    if (random) {
      res$random$RD <- x$TE.random
      res$random$lower.RD <- x$lower.random
      res$random$upper.RD <- x$upper.random
    }
  }
  #
  else if (x$sm == "RR") {
    if (common) {
      res$common$RD <- logRR2rd(x$TE.common, p.c)
      res$common$lower.RD <- logRR2rd(x$lower.common, p.c)
      res$common$upper.RD <- logRR2rd(x$upper.common, p.c)
    }
    #
    if (random) {
      res$random$RD <- logRR2rd(x$TE.random, p.c)
      res$random$lower.RD <- logRR2rd(x$lower.random, p.c)
      res$random$upper.RD <- logRR2rd(x$upper.random, p.c)
    }
  }
  #
  else if (x$sm == "OR") {
    if (common) {
      res$common$RD <- logOR2rd(x$TE.common, p.c)
      res$common$lower.RD <- logOR2rd(x$lower.common, p.c)
      res$common$upper.RD <- logOR2rd(x$upper.common, p.c)
    }
    #
    if (random) {
      res$random$RD <- logOR2rd(x$TE.random, p.c)
      res$random$lower.RD <- logOR2rd(x$lower.random, p.c)
      res$random$upper.RD <- logOR2rd(x$upper.random, p.c)
    }
  }
  #
  else if (x$sm == "HR") {
    if (common) {
      res$common$RD <- logHR2rd(x$TE.common, p.c)
      res$common$lower.RD <- logHR2rd(x$lower.common, p.c)
      res$common$upper.RD <- logHR2rd(x$upper.common, p.c)
    }
    #
    if (random) {
      res$random$RD <- logHR2rd(x$TE.random, p.c)
      res$random$lower.RD <- logHR2rd(x$lower.random, p.c)
      res$random$upper.RD <- logHR2rd(x$upper.random, p.c)
    }
  }
  #
  class(res) <- "rd.meta"
  #
  attr(res, "level") <- x$level.ma
  attr(res, "sm") <- x$sm
  attr(res, "small.values") <- small.values
  #
  attr(res, "pscale") <- pscale
  #
  attr(res, "call") <- match.call()
  attr(res, "version") <- packageDescription("meta")$Version
  
  res
}





#' @rdname rd
#' @method rd default
#' @export

rd.default <- function(x, p.c, sm, lower, upper,
                       level = gs("level.ma"),
                       small.values = "desirable",
                       pscale = 1,
                       transf = FALSE,
                       ...) {
  
  #
  # (1) Check arguments
  #
  if (missing(x))
    stop("Argument 'x' is mandatory.")
  #
  if (missing(sm))
    stop("Argument 'sm' is mandatory.")
  #
  sm <- setchar(sm, c("RD", "RR", "OR", "HR"))
  #
  chklevel(level)
  #
  if (missing(p.c)) {
    if (sm == "RD")
      p.c <- NA
    else
      stop("Argument 'p.c' is mandatory.")
  }
  #
  chknumeric(p.c, 0, 1)
  #
  chknumeric(pscale, length = 1)
  #
  missing.lower <- missing(lower)
  if (!missing.lower)
    chklength(lower, length(x), "rd.default")
  #
  missing.upper <- missing(upper)
  if (!missing.upper)
    chklength(upper, length(x), "rd.default")
  #
  # Check length of arguments 'x' and 'p.c'
  #
  if (length(x) != length(p.c)) {
    if (length(x) > length(p.c))
      bad <- length(x) %% length(p.c) != 0
    else
      bad <- length(p.c) %% length(x) != 0
    #
    if (bad)
      stop("Lengths of arguments 'x' and 'p.c' do not match.",
           call. = FALSE)
  }
  #
  small.values <- setsv(small.values)
  #
  chklogical(transf)
  
  
  #
  # (2) Calculate RDs
  #
  res <- data.frame(TE = if (is_relative_effect(sm) & !transf) log(x) else x,
                    lower = NA, upper = NA,
                    p.c = p.c, RD = NA, lower.RD = NA, upper.RD = NA)
  #
  if (missing(lower)) {
    res$lower <- NULL
    res$lower.RD <- NULL
  }
  else {
    res$lower <- lower
    if (is_relative_effect(sm) & !transf)
      res$lower <- log(res$lower)
  }
  #
  if (missing(upper)) {
    res$upper <- NULL
    res$upper.RD <- NULL
  }
  else {
    res$upper <- upper
    if (is_relative_effect(sm) & !transf)
      res$upper <- log(res$upper)
  }
  #
  if (sm == "RD") {
    res$RD <- res$TE
    if (!missing.lower)
      res$lower.RD <- res$lower
    if (!missing.upper)
      res$upper.RD <- res$upper
  }
  #
  else if (sm == "RR") {
    res$RD <- logRR2rd(res$TE, res$p.c)
    if (!missing.lower)
      res$lower.RD <- logRR2rd(res$lower, res$p.c)
    if (!missing.upper)
      res$upper.RD <- logRR2rd(res$upper, res$p.c)
  }
  #
  else if (sm == "OR") {
    res$RD <- logOR2rd(res$TE, res$p.c)
    if (!missing.lower)
      res$lower.RD <- logOR2rd(res$lower, res$p.c)
    if (!missing.upper)
      res$upper.RD <- logOR2rd(res$upper, res$p.c)
  }
  #
  else if (sm == "HR") {
    res$RD <- logHR2rd(res$TE, res$p.c)
    if (!missing.lower)
      res$lower.RD <- logHR2rd(res$lower, res$p.c)
    if (!missing.upper)
      res$upper.RD <- logHR2rd(res$upper, res$p.c)
  }
  #
  class(res) <- c("rd.default", class(res))
  #
  attr(res, "level") <- level
  attr(res, "sm") <- sm
  attr(res, "small.values") <- small.values
  #
  attr(res, "pscale") <- pscale
  #
  attr(res, "call") <- match.call()
  attr(res, "version") <- packageDescription("meta")$Version
  
  res
}





#' @rdname rd
#' @method print rd.meta
#' @export

print.rd.meta <- function(x,
                          common = !is.null(x$common),
                          random = !is.null(x$random),
                          pscale = attributes(x)$pscale,
                          digits = gs("digits"),
                          digits.prop = gs("digits.prop"),
                          big.mark = gs("big.mark"),
                          details = gs("details"),
                          ...) {
  
  
  #
  # (1) Check for rd.meta object
  #
  chkclass(x, "rd.meta")
  
  
  #
  # (2) Check arguments
  #
  chklogical(common)
  chklogical(random)
  #
  chknumeric(pscale, length = 1)
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.prop, min = 0, length = 1)
  #
  chklogical(details)
  
  
  sm <- attributes(x)$sm
  sv <- attributes(x)$small.values
  #
  ci.lab <- paste0(round(100 * replaceNULL(attributes(x)$level), 1), "% CI")
  p.lab <- if (sm == "HR") "Surv.c" else "Risk.c"
  
  if (common) {
    res.c <- x$common
    #
    with.ci <- sum(!is.na(res.c$lower)) > 0 & sum(!is.na(res.c$upper)) > 0
    #
    if (sm != "RD") {
      res.c$TE <- exp(res.c$TE)
      #
      if (with.ci) {
        res.c$lower <- exp(res.c$lower)
        res.c$upper <- exp(res.c$upper)
      }
    }
    #
    if (all(res.c$RD <= 0) & any(res.c$RD < 0) & sv == "desirable") {
      lab.c <- "ARR"
      title.c <- "Absolute risk reduction"
      #
      res.c$RD <- -res.c$RD
      if (with.ci) {
        tmp.u <- -res.c$upper.RD
        res.c$upper.RD <- -res.c$lower.RD
        res.c$lower.RD <- tmp.u
      }
    }
    else if (all(res.c$RD <= 0) & any(res.c$RD < 0) & sv == "undesirable") {
      lab.c <- "ABR"
      title.c <- "Absolute benefit reduction"
      #
      res.c$RD <- -res.c$RD
      #
      if (with.ci) {
        tmp.u <- -res.c$upper.RD
        res.c$upper.RD <- -res.c$lower.RD
        res.c$lower.RD <- tmp.u
      }
    }
    else if (all(res.c$RD >= 0) & any(res.c$RD > 0) & sv == "desirable") {
      lab.c <- "ARI"
      title.c <- "Absolute risk increase"
    }
    else if (all(res.c$RD >= 0) & any(res.c$RD > 0) & sv == "undesirable") {
      lab.c <- "ABI"
      title.c <- "Absolute benefit increase"
    }
    else {
      lab.c <- "RD"
      title.c <- "Risk difference"
    }
    #
    res.c$TE <- formatN(round(res.c$TE, digits), digits, big.mark = big.mark)
    #
    if (sm != "HR")
      res.c$p.c <- formatN(round(pscale * res.c$p.c, digits.prop),
                           digits.prop, big.mark = big.mark)
    else
      res.c$p.c <- formatN(round(res.c$p.c, 2), 2, big.mark = big.mark)
    #
    res.c$RD <- formatN(round(pscale * res.c$RD, digits.prop),
                        digits.prop, big.mark = big.mark)
    #
    if (with.ci) {
      res.c$lower <-
        formatCI(
          formatN(round(res.c$lower, digits), digits, big.mark = big.mark),
          formatN(round(res.c$upper, digits), digits, big.mark = big.mark))
      res.c$upper <- NULL
      #
      res.c$lower.RD <-
        formatCI(formatN(round(pscale * res.c$lower.RD, digits.prop),
                         digits.prop, big.mark = big.mark),
                 formatN(round(pscale * res.c$upper.RD, digits.prop),
                         digits.prop, big.mark = big.mark))
      res.c$upper.RD <- NULL
    }
    else {
      res.c$lower <- NULL
      res.c$upper <- NULL
      #
      res.c$lower.RD <- NULL
      res.c$upper.RD <- NULL
    }
    #
    no_column_label <- TE <- p.c <- RD <- NULL
    #
    res.c %<>%
      mutate(no_column_label = "") %>%
      relocate(no_column_label, .after = "p.c") %>%
      rename(!!lab.c := RD) %>%
      rename(!!sm := TE) %>%
      rename(!!p.lab := p.c) %>% as.matrix()
    #
    if (with.ci) {
      colnames(res.c)[colnames(res.c) == "lower"] <- ci.lab
      colnames(res.c)[colnames(res.c) == "lower.RD"] <- ci.lab
    }
    #
    colnames(res.c)[colnames(res.c) == "no_column_label"] <- ""
    #
    if (is.null(attributes(x)$rd.default))
      title <- paste0(title.c, " (", tolower(gs("text.common")), ")")
    else
      title <- title.c
    cat(paste0(title, "\n\n"))
    #
    prmatrix(res.c, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(res.c)))
    #
    if (random)
      cat("\n")
  }
  #
  if (random) {
    res.r <- x$random
    #
    with.ci <- sum(!is.na(res.r$lower)) > 0 & sum(!is.na(res.r$upper)) > 0
    #
    if (sm != "RD") {
      res.r$TE <- exp(res.r$TE)
      #
      if (with.ci) {
        res.r$lower <- exp(res.r$lower)
        res.r$upper <- exp(res.r$upper)
      }
    }
    #
    if (all(res.r$RD <= 0) & any(res.r$RD < 0) & sv == "desirable") {
      lab.r <- "ARR"
      title.r <- "Absolute risk reduction"
      #
      res.r$RD <- -res.r$RD
      if (with.ci) {
        tmp.u <- -res.r$upper.RD
        res.r$upper.RD <- -res.r$lower.RD
        res.r$lower.RD <- tmp.u
      }
    }
    else if (all(res.r$RD <= 0) & any(res.r$RD < 0) & sv == "undesirable") {
      lab.r <- "ABR"
      title.r <- "Absolute benefit reduction"
      #
      res.r$RD <- -res.r$RD
      #
      if (with.ci) {
        tmp.u <- -res.r$upper.RD
        res.r$upper.RD <- -res.r$lower.RD
        res.r$lower.RD <- tmp.u
      }
    }
    else if (all(res.r$RD >= 0) & any(res.r$RD > 0) & sv == "desirable") {
      lab.r <- "ARI"
      title.r <- "Absolute risk increase"
    }
    else if (all(res.r$RD >= 0) & any(res.r$RD > 0) & sv == "undesirable") {
      lab.r <- "ABI"
      title.r <- "Absolute benefit increase"
    }
    else {
      lab.r <- "RD"
      title.r <- "Risk difference"
    }
    #
    res.r$TE <- formatN(round(res.r$TE, digits), digits, big.mark = big.mark)
    #
    if (sm != "HR")
      res.r$p.c <- formatN(round(pscale * res.r$p.c, digits.prop),
                           digits.prop, big.mark = big.mark)
    else
      res.r$p.c <- formatN(round(res.r$p.c, 2), 2, big.mark = big.mark)
    #
    res.r$RD <- formatN(round(pscale * res.r$RD, digits.prop),
                        digits.prop, big.mark = big.mark)
    #
    if (with.ci) {
      res.r$lower <-
        formatCI(
          formatN(round(res.r$lower, digits), digits, big.mark = big.mark),
          formatN(round(res.r$upper, digits), digits, big.mark = big.mark))
      res.r$upper <- NULL
      #
      res.r$lower.RD <-
        formatCI(formatN(round(pscale * res.r$lower.RD, digits.prop),
                         digits.prop, big.mark = big.mark),
                 formatN(round(pscale * res.r$upper.RD, digits.prop),
                         digits.prop, big.mark = big.mark))
      res.r$upper.RD <- NULL
    }
    else {
      res.r$lower <- NULL
      res.r$upper <- NULL
      #
      res.r$lower.RD <- NULL
      res.r$upper.RD <- NULL
    }
    #
    no_column_label <- TE <- p.c <- RD <- NULL
    #
    res.r %<>%
      mutate(no_column_label = "") %>%
      relocate(no_column_label, .after = "p.c") %>%
      rename(!!lab.r := RD) %>%
      rename(!!sm := TE) %>%
      rename(!!p.lab := p.c) %>% as.matrix()
    #
    if (with.ci) {
      colnames(res.r)[colnames(res.r) == "lower"] <- ci.lab
      colnames(res.r)[colnames(res.r) == "lower.RD"] <- ci.lab
    }
    #
    colnames(res.r)[colnames(res.r) == "no_column_label"] <- ""
    #
    if (is.null(attributes(x)$rd.default))
      title <- paste0(title.r, " (", tolower(gs("text.random")), ")")
    else
      title <- title.r
    cat(paste0(title, "\n\n"))
    #
    prmatrix(res.r, quote = FALSE, right = TRUE,
             rowlab = rep("", nrow(res.r)))
  }
  
  if (details) {
    if (sm == "RD")
      abs.details <- ""
    else {
      abs.details <-
        paste0("\n- Absolute effect",
               if (length(unique(c(if (common) x$common$RD,
                                   if (random) x$random$RD))) > 1) "s",
               " calculated from ",
               tolower(xlab_meta(sm, TRUE)),
               if (length(unique(c(if (common) x$common$TE,
                                   if (random) x$random$TE))) > 1) "s",
               " (argument 'sm = \"", sm, "\"')")
    }
    #
    p.details <- paste0("\n- ", p.lab, ": ",
                        if (sm != "HR") "risk " else "survival probability ",
                        "in control group")
    #
    if (pscale != 1)
      ps.details <-
        paste0("\n- Events per ",
               format(pscale, scientific = FALSE, big.mark = big.mark),
               " observations")
    else
      ps.details <- ""
    #
    if (sv == "desirable") {
      sv.details <-
        paste0("\n- ", sm, " < ",
               if (sm == "RD") "0" else "1",
               " correspond to a beneficial effect of the intervention")
    }
    else {
      sv.details <-
        paste0("\n- ", sm, " > ",
               if (sm == "RD") "0" else "1",
               " corresponds to a beneficial effect of the intervention")
    }
    #
    sv.details <-
      paste0(sv.details, "\n  (argument 'small.values = \"", sv, "\"')")
    #
    cat(paste0("\nDetails:",
               abs.details, p.details, ps.details, sv.details, "\n"))
  }
  
  invisible(NULL)
}





#' @rdname rd
#' @method print rd.default
#' @export

print.rd.default <- function(x,
                             pscale = attributes(x)$pscale,
                             digits = gs("digits"),
                             digits.prop = gs("digits.prop"),
                             big.mark = gs("big.mark"),
                             ...) {
  
  #
  # (1) Check for rd.default object
  #
  chkclass(x, "rd.default")
  
  res <- vector("list")
  res$common <- x
  class(res$common) <- "data.frame"
  class(res) <- "rd.meta"
  #
  attr(res, "level") <- attributes(x)$level
  attr(res, "sm") <- attributes(x)$sm
  attr(res, "small.values") <- attributes(x)$small.values
  #
  attr(res, "pscale") <- attributes(x)$pscale
  #
  attr(res, "rd.default") <- TRUE
  #
  print(res, pscale = pscale, digits = digits,  digits.prop = digits.prop,
        big.mark = big.mark)
  #
  invisible(NULL)
}





#
# Auxiliary R functions
#
logRR2rd <- function(x, p.c)
  exp(x) * p.c - p.c
#
logOR2rd <- function(x, p.c)
  exp(x) * p.c / (1 - p.c + exp(x) * p.c) - p.c
#
logHR2rd <- function(x, Surv.c)
  (1 - Surv.c^exp(x)) - (1 - Surv.c)
