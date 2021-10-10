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
#' @param level.ma The level used to calculate confidence intervals
#'   for pooled estimates.
#' @param fixed A logical indicating whether a fixed effect
#'   meta-analysis should be conducted.
#' @param random A logical indicating whether a random effects
#'   meta-analysis should be conducted.
#' @param hakn A logical indicating whether the method by Hartung and
#'   Knapp should be used to adjust test statistics and confidence
#'   intervals.
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2} and its
#'   square root \eqn{\tau}. Either \code{"DL"}, \code{"PM"},
#'   \code{"REML"}, \code{"ML"}, \code{"HS"}, \code{"SJ"},
#'   \code{"HE"}, or \code{"EB"}, can be abbreviated.
#' @param method.tau.ci A character string indicating which method is
#'   used to estimate the confidence interval of \eqn{\tau^2} and
#'   \eqn{\tau}. Either \code{"QP"}, \code{"BJ"}, or \code{"J"}, or
#'   \code{""}, can be abbreviated.
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
#' @param test.subgroup A logical value indicating whether to print
#'   results of test for subgroup differences.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratios and results
#'   for \code{sm="ZCOR"} are printed as correlations rather than
#'   Fisher's z transformed correlations, for example.
#' @param text.fixed A character string used in printouts and forest
#'   plot to label the pooled fixed effect estimate.
#' @param text.random A character string used in printouts and forest
#'   plot to label the pooled random effects estimate.
#' @param text.predict A character string used in printouts and forest
#'   plot to label the prediction interval.
#' @param text.w.fixed A character string used to label weights of
#'   fixed effect model.
#' @param text.w.random A character string used to label weights of
#'   random effects model.
#' @param title Title of meta-analysis / systematic review.
#' @param complab Comparison label.
#' @param outclab Outcome label.
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if \code{incr} is added to studies with zero cell
#'   frequencies).
#' @param \dots Additional arguments (to catch deprecated arguments).
#' 
#' @details
#' Cochrane Intervention reviews are based on the comparison of two
#' interventions. Each Cochrane Intervention review can have a
#' variable number of comparisons. For each comparison, a variable
#' number of outcomes can be define. For each outcome, a seperate
#' meta-analysis is conducted. Review Manager 5 (RevMan 5) is the
#' current software used for preparing and maintaining Cochrane
#' Reviews
#' (\url{https://training.cochrane.org/online-learning/core-software-cochrane-reviews/revman}).
#' 
#' This wrapper function can be used to perform meta-analysis for a
#' single outcome of a Cochrane Intervention review. Internally, R
#' functions \code{\link{metabin}}, \code{\link{metacont}}, and
#' \code{\link{metagen}} are called - depending on the definition of
#' the outcome in RevMan 5.
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
#' \emph{Review Manager (RevMan)} [Computer program]. Version 5.4.
#' The Cochrane Collaboration, 2020
#' 
#' @examples
#' # Locate export data file "Fleiss1993_CR.csv"
#' # in sub-directory of package "meta"
#' #
#' filename <- system.file("extdata", "Fleiss1993_CR.csv", package = "meta")
#' #
#' Fleiss1993_CR <- read.rm5(filename)
#' 
#' # Choose RevMan 5 settings and store old settings
#' #
#' oldset <- settings.meta("revman5")
#' 
#' # Same result as R command example(Fleiss1993bin)
#' #
#' metacr(Fleiss1993_CR)
#' 
#' # Same result as R command example(Fleiss1993cont)
#' #
#' metacr(Fleiss1993_CR, 1, 2)
#' forest(metacr(Fleiss1993_CR, 1, 2))
#' 
#' # Change summary measure to RR
#' #
#' m1 <- metacr(Fleiss1993_CR)
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
                   level = gs("level"), level.ma = gs("level.ma"),
                   fixed, random,
                   ##
                   hakn = FALSE,
                   method.tau = "DL",
                   method.tau.ci = gs("method.tau.ci"),
                   tau.common = FALSE,
                   ##
                   prediction = gs("prediction"),
                   level.predict = gs("level.predict"),
                   ##
                   swap.events, logscale,
                   ##
                   backtransf = gs("backtransf"),
                   ##
                   test.subgroup,
                   ##
                   text.fixed = gs("text.fixed"),
                   text.random = gs("text.random"),
                   text.predict = gs("text.predict"),
                   text.w.fixed = gs("text.w.fixed"),
                   text.w.random = gs("text.w.random"),
                   ##
                   title, complab, outclab,
                   ##
                   keepdata = gs("keepdata"),
                   warn = FALSE,
                   ...) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chkclass(x, "rm5")
  ##
  if (is.numeric(comp.no))
    chknumeric(comp.no, length = 1)
  else
    chkchar(comp.no, length = 1)
  ##
  if (is.numeric(outcome.no))
    chknumeric(outcome.no, length = 1)
  else
    chkchar(outcome.no, length = 1)
  ##
  chklevel(level)
  chklevel(level.ma)
  ##
  chklogical(hakn)
  method.tau <- setchar(method.tau, .settings$meth4tau)
  if (is.null(method.tau.ci))
    method.tau.ci <- if (method.tau == "DL") "J" else "QP"
  method.tau.ci <- setchar(method.tau.ci, .settings$meth4tau.ci)
  chklogical(tau.common)
  ##
  chklogical(prediction)
  chklevel(level.predict)
  ##
  if (!missing(swap.events))
    chklogical(swap.events)
  ##
  chklogical(backtransf)
  ##
  if (!is.null(text.fixed))
    chkchar(text.fixed, length = 1)
  if (!is.null(text.random))
    chkchar(text.random, length = 1)
  if (!is.null(text.predict))
    chkchar(text.predict, length = 1)
  if (!is.null(text.w.fixed))
    chkchar(text.w.fixed, length = 1)
  if (!is.null(text.w.random))
    chkchar(text.w.random, length = 1)
  ##
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
  overall <- replaceNULL(unique(x$overall[sel]), TRUE)
  if (is.na(overall))
    overall <- TRUE
  type <- unique(x$type[sel])
  if (missing(test.subgroup))
    test.subgroup <- unique(x$test.subgroup[sel])
  ##
  if (missing(method))
    method <- unique(x$method[sel])
  else
    method <- setchar(method, c("Inverse", "MH", "Peto"))
  chkchar(method, length = 1)
  ##
  if (missing(sm))
    sm <- unique(x$sm[sel])
  chkchar(sm, length = 1)
  if (sm == "PETO_OR")
    sm <- "OR"
  ##
  if (missing(fixed))
    fixed  <- unique(x$fixed[sel])
  ##
  if (missing(random))
    random <- unique(x$random[sel])
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
    if (!unique(x$logscale[sel])) {
      x$TE[sel] <- log(x$TE[sel])
      logscale <- FALSE
    }
    else
      logscale <- TRUE
  }
  else {
    if (!missing(logscale)) {
      if (!logscale)
        x$TE[sel] <- log(x$TE[sel])
    }
    else {
      if ((type == "I" & method != "Peto") & is.relative.effect(sm)) {
        warning("Assuming that values for 'TE' are on log scale. ",
                "Please use argument 'logscale = FALSE' if ",
                "values are on natural scale.",
                call. = FALSE)
        logscale <- TRUE
      }
      else
        logscale <- NA
    }
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
                 "overall", "test.subgroup",
                 "type", "method", "sm", "model", "fixed", "random",
                 "outclab", "k",
                 "event.e.pooled", "n.e.pooled",
                 "event.c.pooled", "n.c.pooled",
                 "TE.pooled", "lower.pooled", "upper.pooled",
                 "weight.pooled", "Z.pooled", "pval.TE.pooled",
                 "Q", "pval.Q", "I2", "tau2", "Q.w", "pval.Q.w", "I2.w",
                 "swap.events", "enter.n", "logscale", "label.e", "label.c",
                 "label.left", "label.right", "complab", "sel")
  ##
  varnames <- names(x)[!(names(x) %in% dropnames)]
  ##
  Q.Cochrane <- if (method == "MH" & method.tau == "DL") TRUE else FALSE
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
                      fixed = fixed, random = random,
                      hakn = hakn,
                      method.tau = method.tau, method.tau.ci = method.tau.ci,
                      tau.common = tau.common,
                      level = level, level.ma = level.ma,
                      prediction = prediction, level.predict = level.predict,
                      overall = overall,
                      subgroup = grplab,
                      subgroup.name = "grp",
                      print.subgroup.name = FALSE,
                      test.subgroup = test.subgroup,
                      backtransf = backtransf,
                      ##
                      text.fixed = text.fixed, text.random = text.random,
                      text.predict = text.predict,
                      text.w.fixed = text.w.fixed,
                      text.w.random = text.w.random,
                      ##
                      title = title,
                      complab = complab, outclab = outclab,
                      label.e = label.e, label.c = label.c,
                      label.left = label.left, label.right = label.right,
                      RR.Cochrane = TRUE,
                      Q.Cochrane = Q.Cochrane,
                      warn = warn, keepdata = keepdata)
      else
        m1 <- metabin(event.e, n.e, event.c, n.c,
                      sm = sm, method = method, studlab = studlab,
                      data = x[sel, varnames],
                      fixed = fixed, random = random,
                      hakn = hakn,
                      method.tau = method.tau, method.tau.ci = method.tau.ci,
                      tau.common = tau.common,
                      level = level, level.ma = level.ma,
                      prediction = prediction, level.predict = level.predict,
                      overall = overall,
                      subgroup = grplab,
                      subgroup.name = "grp",
                      print.subgroup.name = FALSE,
                      test.subgroup = test.subgroup,
                      backtransf = backtransf,
                      ##
                      text.fixed = text.fixed, text.random = text.random,
                      text.predict = text.predict,
                      text.w.fixed = text.w.fixed,
                      text.w.random = text.w.random,
                      ##
                      title = title,
                      complab = complab, outclab = outclab,
                      label.e = label.e, label.c = label.c,
                      label.left = label.left, label.right = label.right,
                      RR.Cochrane = TRUE,
                      Q.Cochrane = Q.Cochrane,
                      warn = warn, keepdata = keepdata)
    }
    ##
    if (type == "C")
      m1 <- metacont(n.e, mean.e, sd.e,
                     n.c, mean.c, sd.c,
                     sm = sm, studlab = studlab,
                     data = x[sel, varnames],
                     fixed = fixed, random = random,
                     hakn = hakn,
                     method.tau = method.tau, method.tau.ci = method.tau.ci,
                     tau.common = tau.common,
                     level = level, level.ma = level.ma,
                     prediction = prediction, level.predict = level.predict,
                     overall = overall,
                     subgroup = grplab,
                     subgroup.name = "grp",
                     print.subgroup.name = FALSE,
                     test.subgroup = test.subgroup,
                     ##
                     text.fixed = text.fixed, text.random = text.random,
                     text.predict = text.predict,
                     text.w.fixed = text.w.fixed,
                     text.w.random = text.w.random,
                     ##
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
                    fixed = fixed, random = random,
                    hakn = hakn,
                    method.tau = method.tau, method.tau.ci = method.tau.ci,
                    tau.common = tau.common,
                    level = level, level.ma = level.ma,
                    prediction = prediction, level.predict = level.predict,
                    overall = overall,
                    subgroup = grplab,
                    subgroup.name = "grp",
                    print.subgroup.name = FALSE,
                    test.subgroup = test.subgroup,
                    backtransf = backtransf,
                    ##
                    text.fixed = text.fixed, text.random = text.random,
                    text.predict = text.predict,
                    text.w.fixed = text.w.fixed,
                    text.w.random = text.w.random,
                    ##
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
                    fixed = fixed, random = random,
                    hakn = hakn,
                    method.tau = method.tau, method.tau.ci = method.tau.ci,
                    tau.common = tau.common,
                    level = level, level.ma = level.ma,
                    prediction = prediction, level.predict = level.predict,
                    overall = overall,
                    subgroup = grplab,
                    subgroup.name = "grp",
                    print.subgroup.name = FALSE,
                    test.subgroup = test.subgroup,
                    n.e = n.e,
                    n.c = n.c,
                    backtransf = backtransf,
                    ##
                    text.fixed = text.fixed, text.random = text.random,
                    text.predict = text.predict,
                    text.w.fixed = text.w.fixed,
                    text.w.random = text.w.random,
                    ##
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
                    fixed = fixed, random = random,
                    hakn = hakn,
                    method.tau = method.tau, method.tau.ci = method.tau.ci,
                    tau.common = tau.common,
                    level = level, level.ma = level.ma,
                    prediction = prediction, level.predict = level.predict,
                    overall = overall,
                    subgroup = grplab,
                    subgroup.name = "grp",
                    print.subgroup.name = FALSE,
                    test.subgroup = test.subgroup,
                    backtransf = backtransf,
                    ##
                    text.fixed = text.fixed, text.random = text.random,
                    text.predict = text.predict,
                    text.w.fixed = text.w.fixed,
                    text.w.random = text.w.random,
                    ##
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
                      fixed = fixed, random = random,
                      hakn = hakn,
                      method.tau = method.tau, method.tau.ci = method.tau.ci,
                      tau.common = tau.common,
                      level = level, level.ma = level.ma,
                      prediction = prediction, level.predict = level.predict,
                      overall = overall,
                      backtransf = backtransf,
                      ##
                      text.fixed = text.fixed, text.random = text.random,
                      text.predict = text.predict,
                      text.w.fixed = text.w.fixed,
                      text.w.random = text.w.random,
                      ##
                      title = title,
                      complab = complab, outclab = outclab,
                      label.e = label.e, label.c = label.c,
                      label.left = label.left, label.right = label.right,
                      RR.Cochrane = TRUE,
                      Q.Cochrane = Q.Cochrane,
                      warn = warn, keepdata = keepdata)
      else
        m1 <- metabin(event.e, n.e, event.c, n.c,
                      sm = sm, method = method, studlab = studlab,
                      data = x[sel, varnames],
                      fixed = fixed, random = random,
                      hakn = hakn,
                      method.tau = method.tau, method.tau.ci = method.tau.ci,
                      tau.common = tau.common,
                      level = level, level.ma = level.ma,
                      prediction = prediction, level.predict = level.predict,
                      overall = overall,
                      backtransf = backtransf,
                      ##
                      text.fixed = text.fixed, text.random = text.random,
                      text.predict = text.predict,
                      text.w.fixed = text.w.fixed,
                      text.w.random = text.w.random,
                      ##
                      title = title,
                      complab = complab, outclab = outclab,
                      label.e = label.e, label.c = label.c,
                      label.left = label.left, label.right = label.right,
                      RR.Cochrane = TRUE,
                      Q.Cochrane = Q.Cochrane,
                      warn = warn, keepdata = keepdata)
    }
    ##
    if (type == "C")
      m1 <- metacont(n.e, mean.e, sd.e,
                     n.c, mean.c, sd.c,
                     sm = sm, studlab = studlab,
                     data = x[sel, varnames],
                     fixed = fixed, random = random,
                     hakn = hakn,
                     method.tau = method.tau, method.tau.ci = method.tau.ci,
                     tau.common = tau.common,
                     level = level, level.ma = level.ma,
                     prediction = prediction, level.predict = level.predict,
                     overall = overall,
                     ##
                     text.fixed = text.fixed, text.random = text.random,
                     text.predict = text.predict,
                     text.w.fixed = text.w.fixed,
                     text.w.random = text.w.random,
                     ##
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
                    fixed = fixed, random = random,
                    hakn = hakn,
                    method.tau = method.tau, method.tau.ci = method.tau.ci,
                    tau.common = tau.common,
                    level = level, level.ma = level.ma,
                    prediction = prediction, level.predict = level.predict,
                    overall = overall,
                    backtransf = backtransf,
                    ##
                    text.fixed = text.fixed, text.random = text.random,
                    text.predict = text.predict,
                    text.w.fixed = text.w.fixed,
                    text.w.random = text.w.random,
                    ##
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
                    fixed = fixed, random = random,
                    hakn = hakn,
                    method.tau = method.tau, method.tau.ci = method.tau.ci,
                    tau.common = tau.common,
                    level = level, level.ma = level.ma,
                    prediction = prediction, level.predict = level.predict,
                    overall = overall,
                    n.e = n.e,
                    n.c = n.c,
                    backtransf = backtransf,
                    ##
                    text.fixed = text.fixed, text.random = text.random,
                    text.predict = text.predict,
                    text.w.fixed = text.w.fixed,
                    text.w.random = text.w.random,
                    ##
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
                    fixed = fixed, random = random,
                    hakn = hakn,
                    method.tau = method.tau, method.tau.ci = method.tau.ci,
                    tau.common = tau.common,
                    level = level, level.ma = level.ma,
                    prediction = prediction, level.predict = level.predict,
                    overall = overall,
                    backtransf = backtransf,
                    ##
                    text.fixed = text.fixed, text.random = text.random,
                    text.predict = text.predict,
                    text.w.fixed = text.w.fixed,
                    text.w.random = text.w.random,
                    ##
                    title = title,
                    complab = complab, outclab = outclab,
                    label.e = label.e, label.c = label.c,
                    label.left = label.left, label.right = label.right,
                    keepdata = keepdata)
  }
  
  if (sm == "OTHER") {
    warning('Meta-analysis not possible for sm = "OTHER".')
    return(NULL)
  }
  
  res <- m1
  ##
  attr(res, "comp.no") <- comp.no
  attr(res, "outcome.no") <- outcome.no
  attr(res, "type") <- type
  ##
  if (type == "D")
    attr(res, "swap.events") <- swap.events
    attr(res, "logscale") <- logscale
  ##
  if (!is.null(x$enter.n))
    attr(res, "enter.n") <- x$enter.n
  ##
  res
}
