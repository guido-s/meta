#' Meta-analysis of outcome data from Cochrane review
#' 
#' @description
#' Wrapper function to perform meta-analysis for a single outcome of a
#' Cochrane Intervention review.
#' 
#' @param x An object of class \code{rm5} created by R function
#'   \code{read.rm5}.
#' @param comp.no Comparison number.
#' @param outcome.no Outcome number.
#' @param method A character string indicating which method is to be
#'   used for pooling of studies. One of \code{"Inverse"},
#'   \code{"MH"}, or \code{"Peto"}, can be abbreviated.
#' @param sm A character string indicating which summary measure
#'   (\code{"RR"}, \code{"OR"}, \code{"RD"}, \code{"ASD"},
#'   \code{"HR"}, \code{"MD"}, or \code{"SMD"}, or \code{"ROM"}) is to
#'   be used for pooling of studies.
#' @param level The level used to calculate confidence intervals for
#'   individual studies.
#' @param level.comb The level used to calculate confidence intervals
#'   for pooled estimates.
#' @param comb.fixed A logical indicating whether a fixed effect
#'   meta-analysis should be conducted.
#' @param comb.random A logical indicating whether a random effects
#'   meta-analysis should be conducted.
#' @param hakn A logical indicating whether the method by Hartung and
#'   Knapp should be used to adjust test statistics and confidence
#'   intervals.
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2}. Either
#'   \code{"DL"}, \code{"PM"}, \code{"REML"}, \code{"ML"},
#'   \code{"HS"}, \code{"SJ"}, \code{"HE"}, or \code{"EB"}, can be
#'   abbreviated.
#' @param tau.common A logical indicating whether tau-squared should
#'   be the same across subgroups.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#' @param level.predict The level used to calculate prediction
#'   interval for a new study.
#' @param swap.events A logical indicating whether events and
#'   non-events should be interchanged.
#' @param logscale A logical indicating whether effect estimates are
#'   entered on log-scale.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratios and results
#'   for \code{sm="ZCOR"} are printed as correlations rather than
#'   Fisher's z transformed correlations, for example.
#' @param title Title of meta-analysis / systematic review.
#' @param complab Comparison label.
#' @param outclab Outcome label.
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if \code{incr} is added to studies with zero cell
#'   frequencies).
#' 
#' @details
#' Cochrane Intervention reviews are based on the comparison of two
#' interventions. Each Cochrane Intervention review can have a
#' variable number of comparisons. For each comparison, a variable
#' number of outcomes can be define. For each outcome, a seperate
#' meta-analysis is conducted. Review Manager 5 (RevMan 5) is the
#' current software used for preparing and maintaining Cochrane
#' Reviews
#' (\url{http://community.cochrane.org/tools/review-production-tools/revman-5}).
#' 
#' This wrapper function can be used to perform meta-analysis for a
#' single outcome of a Cochrane Intervention review. Internally, R
#' functions \code{metabin}, \code{metacont}, and \code{metagen} are
#' called - depending on the definition of the outcome in RevMan 5.
#' 
#' Note, it is recommended to choose the RevMan 5 settings before
#' executing \code{metacr}, i.e., \code{settings.meta("revman5")}.
#' 
#' @return
#' An object of class \code{"meta"} and \code{"metabin"},
#' \code{"metacont"}, or \code{"metagen"} depending on outcome type
#' utilised in Cochrane Intervention review for selected outcome.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{metagen}}, \code{\link{read.rm5}},
#'   \code{\link{settings.meta}}
#' 
#' @references
#' \emph{Review Manager (RevMan)}
#' [Computer program]. Version 5.3.
#' Copenhagen: The Nordic Cochrane Centre, The Cochrane Collaboration, 2014
#' 
#' @examples
#' # Locate export data file "Fleiss93_CR.csv"
#' # in sub-directory of package "meta"
#' #
#' filename <- system.file("extdata", "Fleiss93_CR.csv", package = "meta")
#' #
#' Fleiss93_CR <- read.rm5(filename)
#' 
#' # Choose RevMan 5 settings and store old settings
#' #
#' oldset <- settings.meta("revman5")
#' 
#' # Same result as R command example(Fleiss93)
#' #
#' metacr(Fleiss93_CR)
#' 
#' # Same result as R command example(Fleiss93cont)
#' #
#' metacr(Fleiss93_CR, 1, 2)
#' forest(metacr(Fleiss93_CR, 1, 2))
#' 
#' # Change summary measure to RR
#' #
#' m1 <- metacr(Fleiss93_CR)
#' update(m1, sm="RR")
#' 
#' # Use old settings
#' #
#' settings.meta(oldset)
#' 
#' @export metacr


metacr <- function(x, comp.no = 1, outcome.no = 1,
                   ##
                   method, sm,
                   ##
                   level = gs("level"), level.comb = gs("level.comb"),
                   comb.fixed, comb.random,
                   ##
                   hakn = FALSE,
                   method.tau = "DL",
                   tau.common = FALSE,
                   ##
                   prediction = gs("prediction"),
                   level.predict = gs("level.predict"),
                   ##
                   swap.events, logscale,
                   ##
                   backtransf = gs("backtransf"),
                   title, complab, outclab,
                   ##
                   keepdata = gs("keepdata"),
                   warn = FALSE) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chkclass(x, "rm5")
  ##
  if (!missing(sm))
    chknull(sm)
  chklevel(level)
  chklevel(level.comb)
  ##
  chklogical(hakn)
  method.tau <- setchar(method.tau,
                        c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB"))
  chklogical(tau.common)
  ##
  chklogical(prediction)
  chklevel(level.predict)
  ##
  if (!missing(swap.events))
    chklogical(swap.events)
  ##
  chklogical(backtransf)
  chklogical(keepdata)
  
  
  ##
  ##
  ## (2) Select data for meta-analysis
  ##
  ##
  sel <- x$comp.no == comp.no & x$outcome.no == outcome.no
  ##
  if (sum(sel) == 0) {
    warning("No data available for comp.no = ", comp.no,
            " and outcome.no = ", outcome.no, ".")
    return(NULL)
  }
  ##
  x$sel <- sel
  ##
  ## Additional checks
  ##
  if (missing(title))
    title   <- attributes(x)$title
  ##
  if (missing(complab))
    complab <- unique(x$complab[sel])
  ##
  if (missing(outclab))
    outclab <- unique(x$outclab[sel])
  ##
  label.e <- unique(x$label.e[sel])
  label.c <- unique(x$label.c[sel])
  ##
  label.left  <- unique(x$label.left[sel])
  label.right <- unique(x$label.right[sel])
  ##
  type <- unique(x$type[sel])
  ##
  if (missing(sm))
    sm <- unique(x$sm[sel])
  ##
  if (missing(method))
    method <- unique(x$method[sel])
  else
    method <- setchar(method, c("Inverse", "MH", "Peto"))
  ##
  if (missing(comb.fixed))
    comb.fixed  <- unique(x$comb.fixed[sel])
  ##
  if (missing(comb.random))
    comb.random <- unique(x$comb.random[sel])
  ##  
  if (tau.common & method == "Peto") {
    if (warn)
      warning("Argument 'tau.common' not considered for Peto method.")
    tau.common <- FALSE
  }
  ##
  if (tau.common & method == "MH") {
    if (warn)
      warning("Argument 'tau.common' not considered for Mantel-Haenszel method.")
    tau.common <- FALSE
  }
  
  
  ##
  ##
  ## (3) Calculate results for individual studies
  ##
  ##
  if (!all(is.na(x$logscale[sel]))) {
    if (!unique(x$logscale[sel]))
      x$TE[sel] <- log(x$TE[sel])
  }
  else {
    if (!missing(logscale)) {
      if (!logscale)
        x$TE[sel] <- log(x$TE[sel])
    }
    else
      if ((type == "I" & method != "Peto") & is.relative.effect(sm))
        warning("Assuming that values for 'TE' are on log scale. Please use argument 'logscale = FALSE' if values are on natural scale.")
  }
  
  
  ##
  ##
  ## (4) Do meta-analysis
  ##
  ##
  O.E <- TE <- V <- event.c <- event.e <- grplab <- mean.c <- mean.e <-
    n.c <- n.e <- sd.c <- sd.e <- seTE <- studlab <- NULL
  ##
  dropnames <- c("comp.no", "outcome.no", "group.no",
                 "lower.TE", "upper.TE", "weight",
                 "type", "method", "sm", "model", "comb.fixed", "comb.random",
                 "outclab", "k", "event.e.pooled", "n.e.pooled", "event.c.pooled",
                 "n.c.pooled", "TE.pooled", "lower.pooled", "upper.pooled",
                 "weight.pooled", "Z.pooled", "pval.TE.pooled",
                 "Q", "pval.Q", "I2", "tau2", "Q.w", "pval.Q.w", "I2.w",
                 "swap.events", "enter.n", "logscale", "label.e", "label.c",
                 "label.left", "label.right", "complab", "sel")
  ##
  varnames <- names(x)[!(names(x) %in% dropnames)]
  ##
  if (length(unique(x$group.no[sel])) > 1) {
    if (type == "D") {
      ##
      if (missing(swap.events)) {
        swap.events <- unique(x$swap.events[sel])
        swap.events <- !is.na(swap.events) && swap.events
      }
      ##
      if (swap.events)
        m1 <- metabin(n.e - event.e, n.e, n.c - event.c, n.c,
                      sm = sm, method = method, studlab = studlab,
                      data = x[sel, varnames],
                      comb.fixed = comb.fixed, comb.random = comb.random,
                      method.tau = method.tau, hakn = hakn,
                      tau.common = tau.common,
                      level = level, level.comb = level.comb,
                      prediction = prediction, level.predict = level.predict,
                      byvar = grplab,
                      bylab = "grp",
                      print.byvar = FALSE,
                      backtransf = backtransf,
                      title = title,
                      complab = complab, outclab = outclab,
                      label.e = label.e, label.c = label.c,
                      label.left = label.left, label.right = label.right,
                      RR.cochrane = TRUE, warn = warn,
                      keepdata = keepdata)
      else
        m1 <- metabin(event.e, n.e, event.c, n.c,
                      sm = sm, method = method, studlab = studlab,
                      data = x[sel, varnames],
                      comb.fixed = comb.fixed, comb.random = comb.random,
                      method.tau = method.tau, hakn = hakn,
                      tau.common = tau.common,
                      level = level, level.comb = level.comb,
                      prediction = prediction, level.predict = level.predict,
                      byvar = grplab,
                      bylab = "grp",
                      print.byvar = FALSE,
                      backtransf = backtransf,
                      title = title,
                      complab = complab, outclab = outclab,
                      label.e = label.e, label.c = label.c,
                      label.left = label.left, label.right = label.right,
                      RR.cochrane = TRUE, warn = warn,
                      keepdata = keepdata)
    }
    ##
    if (type == "C")
      m1 <- metacont(n.e, mean.e, sd.e,
                     n.c, mean.c, sd.c,
                     sm = sm, studlab = studlab,
                     data = x[sel, varnames],
                     comb.fixed = comb.fixed, comb.random = comb.random,
                     method.tau = method.tau, hakn = hakn,
                     tau.common = tau.common,
                     level = level, level.comb = level.comb,
                     prediction = prediction, level.predict = level.predict,
                     byvar = grplab,
                     bylab = "grp",
                     print.byvar = FALSE,
                     title = title,
                     complab = complab, outclab = outclab,
                     label.e = label.e, label.c = label.c,
                     label.left = label.left, label.right = label.right,
                     keepdata = keepdata)
    ##
    if (type == "P")
      m1 <- metagen(O.E / V, sqrt(1 / V),
                    sm = sm, studlab = studlab,
                    data = x[sel, varnames],
                    comb.fixed = comb.fixed, comb.random = comb.random,
                    method.tau = method.tau, hakn = hakn,
                    tau.common = tau.common,
                    level = level, level.comb = level.comb,
                    prediction = prediction, level.predict = level.predict,
                    byvar = grplab,
                    bylab = "grp",
                    print.byvar = FALSE,
                    backtransf = backtransf,
                    title = title,
                    complab = complab, outclab = outclab,
                    label.e = label.e, label.c = label.c,
                    label.left = label.left, label.right = label.right,
                    keepdata = keepdata)
    ##
    if (type == "I" & method != "Peto")
      m1 <- metagen(TE, seTE,
                    sm = sm, studlab = studlab,
                    data = x[sel, varnames],
                    comb.fixed = comb.fixed, comb.random = comb.random,
                    method.tau = method.tau, hakn = hakn,
                    tau.common = tau.common,
                    level = level, level.comb = level.comb,
                    prediction = prediction, level.predict = level.predict,
                    byvar = grplab,
                    bylab = "grp",
                    print.byvar = FALSE,
                    n.e = n.e,
                    n.c = n.c,
                    backtransf = backtransf,
                    title = title,
                    complab = complab, outclab = outclab,
                    label.e = label.e, label.c = label.c,
                    label.left = label.left, label.right = label.right,
                    keepdata = keepdata)
    ##
    if (type == "I" & method == "Peto")
      m1 <- metagen(O.E / V, sqrt(1 / V),
                    sm = sm, studlab = studlab,
                    data = x[sel, varnames],
                    comb.fixed = comb.fixed, comb.random = comb.random,
                    method.tau = method.tau, hakn = hakn,
                    tau.common = tau.common,
                    level = level, level.comb = level.comb,
                    prediction = prediction, level.predict = level.predict,
                    byvar = grplab,
                    bylab = "grp",
                    print.byvar = FALSE,
                    backtransf = backtransf,
                    title = title,
                    complab = complab, outclab = outclab,
                    label.e = label.e, label.c = label.c,
                    label.left = label.left, label.right = label.right,
                    keepdata = keepdata)
  }
  else {
    if (type == "D") {
      ##
      if (missing(swap.events)) {
        swap.events <- unique(x$swap.events[sel])
        swap.events <- !is.na(swap.events) && swap.events
      }
      ##
      if (swap.events)
        m1 <- metabin(n.e - event.e, n.e, n.c - event.c, n.c,
                      sm = sm, method = method, studlab = studlab,
                      data = x[sel, varnames],
                      comb.fixed = comb.fixed, comb.random = comb.random,
                      method.tau = method.tau, hakn = hakn,
                      tau.common = tau.common,
                      level = level, level.comb = level.comb,
                      prediction = prediction, level.predict = level.predict,
                      backtransf = backtransf,
                      title = title,
                      complab = complab, outclab = outclab,
                      label.e = label.e, label.c = label.c,
                      label.left = label.left, label.right = label.right,
                      RR.cochrane = TRUE, warn = warn,
                      keepdata = keepdata)
      else
        m1 <- metabin(event.e, n.e, event.c, n.c,
                      sm = sm, method = method, studlab = studlab,
                      data = x[sel, varnames],
                      comb.fixed = comb.fixed, comb.random = comb.random,
                      method.tau = method.tau, hakn = hakn,
                      tau.common = tau.common,
                      level = level, level.comb = level.comb,
                      prediction = prediction, level.predict = level.predict,
                      backtransf = backtransf,
                      title = title,
                      complab = complab, outclab = outclab,
                      label.e = label.e, label.c = label.c,
                      label.left = label.left, label.right = label.right,
                      RR.cochrane = TRUE, warn = warn,
                      keepdata = keepdata)
    }
    ##
    if (type == "C")
      m1 <- metacont(n.e, mean.e, sd.e,
                     n.c, mean.c, sd.c,
                     sm = sm, studlab = studlab,
                     data = x[sel, varnames],
                     comb.fixed = comb.fixed, comb.random = comb.random,
                     method.tau = method.tau, hakn = hakn,
                     tau.common = tau.common,
                     level = level, level.comb = level.comb,
                     prediction = prediction, level.predict = level.predict,
                     title = title,
                     complab = complab, outclab = outclab,
                     label.e = label.e, label.c = label.c,
                     label.left = label.left, label.right = label.right,
                     keepdata = keepdata)
    ##
    if (type == "P")
      m1 <- metagen(O.E / V, sqrt(1 / V),
                    sm = sm, studlab = studlab,
                    data = x[sel, varnames],
                    comb.fixed = comb.fixed, comb.random = comb.random,
                    method.tau = method.tau, hakn = hakn,
                    tau.common = tau.common,
                    level = level, level.comb = level.comb,
                    prediction = prediction, level.predict = level.predict,
                    backtransf = backtransf,
                    title = title,
                    complab = complab, outclab = outclab,
                    label.e = label.e, label.c = label.c,
                    label.left = label.left, label.right = label.right,
                    keepdata = keepdata)
    ##
    if (type == "I" & method != "Peto")
      m1 <- metagen(TE, seTE,
                    sm = sm, studlab = studlab,
                    data = x[sel, varnames],
                    comb.fixed = comb.fixed, comb.random = comb.random,
                    method.tau = method.tau, hakn = hakn,
                    tau.common = tau.common,
                    level = level, level.comb = level.comb,
                    prediction = prediction, level.predict = level.predict,
                    n.e = n.e,
                    n.c = n.c,
                    backtransf = backtransf,
                    title = title,
                    complab = complab, outclab = outclab,
                    label.e = label.e, label.c = label.c,
                    label.left = label.left, label.right = label.right,
                    keepdata = keepdata)
    ##
    if (type == "I" & method == "Peto")
      m1 <- metagen(O.E / V, sqrt(1 / V),
                    sm = sm, studlab = studlab,
                    data = x[sel, varnames],
                    comb.fixed = comb.fixed, comb.random = comb.random,
                    method.tau = method.tau, hakn = hakn,
                    tau.common = tau.common,
                    level = level, level.comb = level.comb,
                    prediction = prediction, level.predict = level.predict,
                    backtransf = backtransf,
                    title = title,
                    complab = complab, outclab = outclab,
                    label.e = label.e, label.c = label.c,
                    label.left = label.left, label.right = label.right,
                    keepdata = keepdata)
  }
  
  
  if (sm == "OTHER") {
    warning('Meta-analysis not possible for sm = "OTHER".')
    res <- NULL
  }
  else
    res <- m1
  
  
  res
}
