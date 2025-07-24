#' Meta-analysis of outcome data from Cochrane review
#' 
#' @description
#' Wrapper function to perform meta-analysis for a single outcome of a
#' Cochrane Intervention review.
#' 
#' @param x An object of class \code{rm5} or \code{cdir} created by R
#'   function \code{read.rm5} or \code{read.cdir} .
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
#' @param common A logical indicating whether a common effect
#'   meta-analysis should be conducted.
#' @param random A logical indicating whether a random effects
#'   meta-analysis should be conducted.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2} and its
#'   square root \eqn{\tau} (see \code{\link{meta-package}}).
#' @param method.tau.ci A character string indicating which method is
#'   used to estimate the confidence interval of \eqn{\tau^2} and
#'   \eqn{\tau} (see \code{\link{meta-package}}).
#' @param level.hetstat The level used to calculate confidence intervals
#'   for heterogeneity statistics.
#' @param tau.common A logical indicating whether tau-squared should
#'   be the same across subgroups.
#' @param method.I2 A character string indicating which method is
#'   used to estimate the heterogeneity statistic I\eqn{^2}. Either
#'   \code{"Q"} or \code{"tau2"}, can be abbreviated
#'   (see \code{\link{meta-package}}).
#' @param level.ma The level used to calculate confidence intervals
#'   for meta-analysis estimates.
#' @param method.common.ci A character string indicating which method
#'   is used to calculate confidence interval and test statistic for
#'   common effect estimate (see \code{\link{meta-package}}).
#' @param method.random.ci A character string indicating which method
#'   is used to calculate confidence interval and test statistic for
#'   random effects estimate (see \code{\link{meta-package}}).
#' @param adhoc.hakn.ci A character string indicating whether an
#'   \emph{ad hoc} variance correction should be applied in the case
#'   of an arbitrarily small Hartung-Knapp variance estimate (see
#'   \code{\link{meta-package}}).
#' @param level.predict The level used to calculate prediction
#'   interval for a new study.
#' @param method.predict A character string indicating which method is
#'   used to calculate a prediction interval (see
#'   \code{\link{meta-package}}).
#' @param adhoc.hakn.pi A character string indicating whether an
#'   \emph{ad hoc} variance correction should be applied for
#'   prediction interval (see \code{\link{meta-package}}).
#' @param seed.predict A numeric value used as seed to calculate
#'   bootstrap prediction interval (see \code{\link{meta-package}}).
#' @param Q.Cochrane A logical indicating if the Mantel-Haenszel
#'   estimate is used in the calculation of the heterogeneity
#'   statistic Q which is implemented in RevMan 5.
#' @param swap.events A logical indicating whether events and
#'   non-events should be interchanged.
#' @param logscale A logical indicating whether effect estimates are
#'   entered on log-scale (ignored for \code{cdir} objects).
#' @param test.subgroup A logical value indicating whether to print
#'   results of test for subgroup differences.
#' @param prediction.subgroup A logical indicating whether prediction
#'   intervals should be printed for subgroups.
#' @param seed.predict.subgroup A numeric vector providing seeds to
#'   calculate bootstrap prediction intervals within subgroups. Must
#'   be of same length as the number of subgroups.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratios and results
#'   for \code{sm="ZCOR"} are printed as correlations rather than
#'   Fisher's z transformed correlations, for example.
#' @param rob A logical indicating whether risk of bias (RoB)
#'   assessment should be considered in meta-analysis (only for
#'   \code{read.cdir} objects).
#' @param tool Risk of bias (RoB) tool (only for \code{read.cdir}
#'   objects).
#' @param categories Possible RoB categories (only for
#'   \code{read.cdir} objects).
#' @param col Colours for RoB categories (only for \code{read.cdir}
#'   objects).
#' @param symbols Corresponding symbols for RoB categories (only for
#'   \code{read.cdir} objects).
#' @param text.common A character string used in printouts and forest
#'   plot to label the pooled common effect estimate.
#' @param text.random A character string used in printouts and forest
#'   plot to label the pooled random effects estimate.
#' @param text.predict A character string used in printouts and forest
#'   plot to label the prediction interval.
#' @param text.w.common A character string used to label weights of
#'   common effect model.
#' @param text.w.random A character string used to label weights of
#'   random effects model.
#' @param title Title of meta-analysis / systematic review.
#' @param complab Comparison label.
#' @param outclab Outcome label.
#' @param label.left Graph label on left side of null effect in forest plot.
#' @param label.right Graph label on right side of null effect in forest plot.
#' @param col.label.left The colour of the graph label on the left side of
#'   the null effect.
#' @param col.label.right The colour of the graph label on the right side of
#'   the null effect.
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if \code{incr} is added to studies with zero cell
#'   frequencies).
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param \dots Additional arguments (to catch deprecated arguments).
#' 
#' @details
#' Cochrane intervention reviews are based on the comparison of two
#' interventions. Each Cochrane intervention review can have a
#' variable number of comparisons. For each comparison, a variable
#' number of outcomes can be define. For each outcome, a separate
#' meta-analysis is conducted. Review Manager 5 (RevMan 5) was the
#' software used for preparing and maintaining Cochrane Reviews.
#' 
#' This wrapper function can be used to perform meta-analysis for a
#' single outcome of a Cochrane intervention review. Internally, R
#' functions \code{\link{metabin}}, \code{\link{metacont}}, and
#' \code{\link{metagen}} are called - depending on the definition of
#' the outcome in RevMan 5.
#'
#' Information on the risk of bias RoB) assessment can be provided
#' with arguments \code{tool}, \code{categories}, \code{col} and
#' \code{symbols}. This is not useful if an overall RoB assessment has
#' been done. In this case use \code{\link{rob}} to add the full
#' flexible RoB information to a \code{\link{metacr}} object.
#' 
#' Note, it is recommended to choose the RevMan 5 settings before
#' executing \code{metacr}, i.e., \code{settings.meta("revman5")}.
#' 
#' @return
#' An object of class \code{"meta"} and - depending on outcome type
#' utilised in Cochrane intervention review for selected outcome -
#' \code{"metabin"}, \code{"metacont"}, or \code{"metagen"} with
#' corresponding generic functions (see \code{\link{meta-object}}).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{meta-package}}, \code{\link{rob}},
#'   \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{metagen}}, \code{\link{read.cdir}},
#'   \code{\link{read.rm5}}, \code{\link{settings.meta}}
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
#' oldset <- settings.meta("revman5", quietly = FALSE)
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
                   method, sm, level = gs("level"),
                   ##
                   common, random,
                   prediction = gs("prediction") | !missing(method.predict),
                   ##
                   method.tau = "DL",
                   method.tau.ci = gs("method.tau.ci"),
                   level.hetstat = gs("level.hetstat"),
                   tau.common = FALSE,
                   #
                   method.I2 = gs("method.I2"),
                   #
                   level.ma = gs("level.ma"),
                   method.common.ci = "classic",
                   method.random.ci = "classic",
                   adhoc.hakn.ci = gs("adhoc.hakn.ci"),
                   ##
                   level.predict = gs("level.predict"),
                   method.predict = gs("method.predict"),
                   adhoc.hakn.pi = gs("adhoc.hakn.pi"),
                   seed.predict = NULL,
                   ##
                   Q.Cochrane, swap.events, logscale,
                   ##
                   backtransf = gs("backtransf"),
                   ##
                   test.subgroup,
                   prediction.subgroup = gs("prediction.subgroup"),
                   seed.predict.subgroup = NULL,
                   ##
                   rob = NULL,
                   tool = NULL,
                   categories = NULL,
                   col = NULL,
                   symbols = NULL,
                   ##
                   text.common = gs("text.common"),
                   text.random = gs("text.random"),
                   text.predict = gs("text.predict"),
                   text.w.common = gs("text.w.common"),
                   text.w.random = gs("text.w.random"),
                   ##
                   title, complab, outclab,
                   #
                   label.left, label.right,
                   col.label.left = gs("col.label.left"),
                   col.label.right = gs("col.label.right"),
                   #
                   keepdata = gs("keepdata"),
                   warn = FALSE,
                   warn.deprecated = gs("warn.deprecated"),
                   ...) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chkclass(x, c("rm5", "cdir"))
  ##
  if (inherits(x, "cdir")) {
    x.list <- x
    if (!is.null(x$rob)) {
      robdata <- x$rob
      ##
      if (is.null(rob))
        rob <- attr(robdata, "rob")
      chklogical(rob)
      ##     
      domains <- attr(robdata, "domains")
      n.domains <- length(domains)
      ##
      if (is.null(tool))
        tool <- attr(robdata, "tool")
      if (is.null(categories))
        categories <- attr(robdata, "categories")
      if (is.null(col))
        col <- attr(robdata, "col")
      if (is.null(symbols))
        symbols <- attr(robdata, "symbols")
      ##
      robdata <- robdata[, !grepl(".details", names(robdata))]
    }
    else
      rob <- FALSE
    ##
    x <- x$data
    ##
    if (rob)
      x <- merge(x, robdata, all.x = TRUE)
  }
  else
    rob <- FALSE
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
  chklogical(warn.deprecated)
  args <- list(...)
  ##
  missing.common <- missing(common)
  missing.fixed <- is.na(argid(names(args), "fixed"))
  if (!missing.fixed) {
    common <-
      deprecated(common, missing.common, args, "fixed",
                 warn.deprecated)
    missing.common <- FALSE
  }
  ##
  chklevel(level)
  chklevel(level.ma)
  #
  method.common.ci <- setchar(method.common.ci, gs("meth4common.ci"))
  #
  method.tau <- setchar(method.tau, gs("meth4tau"))
  if (is.null(method.tau.ci))
    method.tau.ci <- if (method.tau == "DL") "J" else "QP"
  method.tau.ci <- setchar(method.tau.ci, gs("meth4tau.ci"))
  chklogical(tau.common)
  #
  method.I2 <- setchar(method.I2, gs("meth4i2"))
  #
  chklogical(prediction)
  chklevel(level.predict)
  chklogical(prediction.subgroup)
  if (!is.null(seed.predict))
    chknumeric(seed.predict)
  if (!is.null(seed.predict.subgroup))
    chknumeric(seed.predict.subgroup)
  ##
  missing.Q.Cochrane <- missing(Q.Cochrane)
  if (!missing.Q.Cochrane)
    chklogical(Q.Cochrane)
  ##
  if (!missing(swap.events))
    chklogical(swap.events)
  ##
  chklogical(backtransf)
  ##
  text.common <-
    deprecated(text.common, missing(text.common), args, "text.fixed",
               warn.deprecated)
  if (!is.null(text.common))
    chkchar(text.common, length = 1)
  ##
  if (!is.null(text.random))
    chkchar(text.random, length = 1)
  if (!is.null(text.predict))
    chkchar(text.predict, length = 1)
  ##
  text.w.common <-
    deprecated(text.w.common, missing(text.w.common), args, "text.w.fixed",
               warn.deprecated)
  if (!is.null(text.w.common))
    chkchar(text.w.common, length = 1)
  ##
  if (!is.null(text.w.random))
    chkchar(text.w.random, length = 1)
  ##
  chklogical(keepdata)
  chklogical(warn)
  ##
  method.random.ci <-
    deprecated(method.random.ci, missing(method.random.ci),
               args, "hakn", warn.deprecated)
  if (is.logical(method.random.ci))
    if (method.random.ci)
      method.random.ci <- "HK"
    else
      method.random.ci <- "classic"
  method.random.ci <- setchar(method.random.ci, gs("meth4random.ci"))
  
  
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
  #
  if (missing(label.left))
    label.left  <- unique(x$label.left[sel])
  else
    chkchar(label.left, length = 1)
  #
  if (missing(label.right))
    label.right  <- unique(x$label.right[sel])
  else
    chkchar(label.right, length = 1)
  #
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
  if (missing.common)
    common <- unique(x$common[sel])
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

  if (inherits(x, "rm5")) {
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
        if ((type == "I" & method != "Peto") & is_relative_effect(sm)) {
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
  }
  else
    logscale <- NULL
  
  
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
                 "type", "method", "sm", "model", "common", "random",
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
  if (missing.Q.Cochrane)
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
                      studlab = studlab,
                      data = x[sel, varnames],
                      ##
                      method = method, sm = sm, level = level,
                      ##
                      common = common, random = random,
                      overall = overall,
                      prediction = prediction,
                      ##
                      method.tau = method.tau, method.tau.ci = method.tau.ci,
                      level.hetstat = level.hetstat,
                      tau.common = tau.common,
                      #
                      method.I2 = method.I2,
                      #
                      level.ma = level.ma,
                      method.common.ci = method.common.ci,
                      method.random.ci = method.random.ci,
                      adhoc.hakn.ci = adhoc.hakn.ci,
                      ##
                      level.predict = level.predict,
                      method.predict = method.predict,
                      adhoc.hakn.pi = adhoc.hakn.pi,
                      seed.predict = seed.predict,
                      ##
                      subgroup = grplab,
                      subgroup.name = "grp",
                      print.subgroup.name = FALSE,
                      test.subgroup = test.subgroup,
                      prediction.subgroup = prediction.subgroup,
                      seed.predict.subgroup = seed.predict.subgroup,
                      ##
                      backtransf = backtransf,
                      ##
                      text.common = text.common, text.random = text.random,
                      text.predict = text.predict,
                      text.w.common = text.w.common,
                      text.w.random = text.w.random,
                      ##
                      title = title,
                      complab = complab, outclab = outclab,
                      #
                      label.e = label.e, label.c = label.c,
                      label.left = label.left, label.right = label.right,
                      col.label.left = col.label.left,
                      col.label.right = col.label.right,
                      #
                      RR.Cochrane = TRUE,
                      Q.Cochrane = Q.Cochrane,
                      ##
                      warn = warn, keepdata = keepdata)
      else
        m1 <- metabin(event.e, n.e, event.c, n.c, studlab = studlab,
                      data = x[sel, varnames],
                      ##
                      method = method, sm = sm, level = level,
                      ##
                      common = common, random = random,
                      overall = overall,
                      prediction = prediction,
                      ##
                      method.tau = method.tau, method.tau.ci = method.tau.ci,
                      level.hetstat = level.hetstat,
                      tau.common = tau.common,
                      #
                      method.I2 = method.I2,
                      #
                      level.ma = level.ma,
                      method.common.ci = method.common.ci,
                      method.random.ci = method.random.ci,
                      adhoc.hakn.ci = adhoc.hakn.ci,
                      ##
                      level.predict = level.predict,
                      method.predict = method.predict,
                      adhoc.hakn.pi = adhoc.hakn.pi,
                      seed.predict = seed.predict,
                      ##
                      subgroup = grplab,
                      subgroup.name = "grp",
                      print.subgroup.name = FALSE,
                      test.subgroup = test.subgroup,
                      prediction.subgroup = prediction.subgroup,
                      seed.predict.subgroup = seed.predict.subgroup,
                      ##
                      backtransf = backtransf,
                      ##
                      text.common = text.common, text.random = text.random,
                      text.predict = text.predict,
                      text.w.common = text.w.common,
                      text.w.random = text.w.random,
                      ##
                      title = title,
                      complab = complab, outclab = outclab,
                      #
                      label.e = label.e, label.c = label.c,
                      label.left = label.left, label.right = label.right,
                      col.label.left = col.label.left,
                      col.label.right = col.label.right,
                      #
                      RR.Cochrane = TRUE,
                      Q.Cochrane = Q.Cochrane,
                      ##
                      warn = warn, keepdata = keepdata)
    }
    ##
    if (type == "C")
      m1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, studlab = studlab,
                     data = x[sel, varnames],
                     ##
                     sm = sm, level = level,
                     ##
                     common = common, random = random,
                     overall = overall,
                     prediction = prediction,
                     ##
                     method.tau = method.tau, method.tau.ci = method.tau.ci,
                     level.hetstat = level.hetstat,
                     tau.common = tau.common,
                     #
                     method.I2 = method.I2,
                     #
                     level.ma = level.ma,
                     method.common.ci = method.common.ci,
                     method.random.ci = method.random.ci,
                     adhoc.hakn.ci = adhoc.hakn.ci,
                     ##
                     level.predict = level.predict,
                     method.predict = method.predict,
                     adhoc.hakn.pi = adhoc.hakn.pi,
                     seed.predict = seed.predict,
                     ##
                     subgroup = grplab,
                     subgroup.name = "grp",
                     print.subgroup.name = FALSE,
                     test.subgroup = test.subgroup,
                     prediction.subgroup = prediction.subgroup,
                     seed.predict.subgroup = seed.predict.subgroup,
                     ##
                     text.common = text.common, text.random = text.random,
                     text.predict = text.predict,
                     text.w.common = text.w.common,
                     text.w.random = text.w.random,
                     ##
                     title = title,
                     complab = complab, outclab = outclab,
                     #
                     label.e = label.e, label.c = label.c,
                     label.left = label.left, label.right = label.right,
                     col.label.left = col.label.left,
                     col.label.right = col.label.right,
                     #
                     keepdata = keepdata)
    ##
    if (type == "P")
      m1 <- metagen(O.E / V, sqrt(1 / V), studlab = studlab,
                    data = x[sel, varnames],
                    ##
                    sm = sm, level = level,
                    ##
                    common = common, random = random,
                    overall = overall,
                    prediction = prediction,
                    ##
                    method.tau = method.tau, method.tau.ci = method.tau.ci,
                    level.hetstat = level.hetstat,
                    tau.common = tau.common,
                    #
                    method.I2 = method.I2,
                    #
                    level.ma = level.ma,
                    method.common.ci = method.common.ci,
                    method.random.ci = method.random.ci,
                    adhoc.hakn.ci = adhoc.hakn.ci,
                    ##
                    level.predict = level.predict,
                    method.predict = method.predict,
                    adhoc.hakn.pi = adhoc.hakn.pi,
                    seed.predict = seed.predict,
                    ##
                    subgroup = grplab,
                    subgroup.name = "grp",
                    print.subgroup.name = FALSE,
                    test.subgroup = test.subgroup,
                    prediction.subgroup = prediction.subgroup,
                    seed.predict.subgroup = seed.predict.subgroup,
                    ##
                    backtransf = backtransf,
                    ##
                    text.common = text.common, text.random = text.random,
                    text.predict = text.predict,
                    text.w.common = text.w.common,
                    text.w.random = text.w.random,
                    ##
                    title = title,
                    complab = complab, outclab = outclab,
                    #
                    label.e = label.e, label.c = label.c,
                    label.left = label.left, label.right = label.right,
                    col.label.left = col.label.left,
                    col.label.right = col.label.right,
                    #
                    keepdata = keepdata)
    ##
    if (type == "I" & method != "Peto")
      m1 <- metagen(TE, seTE, studlab = studlab,
                    data = x[sel, varnames],
                    ##
                    sm = sm, level = level,
                    ##
                    common = common, random = random,
                    overall = overall,
                    prediction = prediction,
                    ##
                    method.tau = method.tau, method.tau.ci = method.tau.ci,
                    level.hetstat = level.hetstat,
                    tau.common = tau.common,
                    #
                    method.I2 = method.I2,
                    #
                    level.ma = level.ma,
                    method.common.ci = method.common.ci,
                    method.random.ci = method.random.ci,
                    adhoc.hakn.ci = adhoc.hakn.ci,
                    ##
                    level.predict = level.predict,
                    method.predict = method.predict,
                    adhoc.hakn.pi = adhoc.hakn.pi,
                    seed.predict = seed.predict,
                    ##
                    subgroup = grplab,
                    subgroup.name = "grp",
                    print.subgroup.name = FALSE,
                    test.subgroup = test.subgroup,
                    prediction.subgroup = prediction.subgroup,
                    seed.predict.subgroup = seed.predict.subgroup,
                    ##
                    n.e = n.e,
                    n.c = n.c,
                    ##
                    backtransf = backtransf,
                    ##
                    text.common = text.common, text.random = text.random,
                    text.predict = text.predict,
                    text.w.common = text.w.common,
                    text.w.random = text.w.random,
                    ##
                    title = title,
                    complab = complab, outclab = outclab,
                    #
                    label.e = label.e, label.c = label.c,
                    label.left = label.left, label.right = label.right,
                    col.label.left = col.label.left,
                    col.label.right = col.label.right,
                    #
                    keepdata = keepdata)
    ##
    if (type == "I" & method == "Peto")
      m1 <- metagen(O.E / V, sqrt(1 / V), studlab = studlab,
                    data = x[sel, varnames],
                    ##
                    sm = sm, level = level,
                    ##
                    common = common, random = random,
                    overall = overall,
                    prediction = prediction,
                    ##
                    method.tau = method.tau, method.tau.ci = method.tau.ci,
                    level.hetstat = level.hetstat,
                    tau.common = tau.common,
                    #
                    method.I2 = method.I2,
                    #
                    level.ma = level.ma,
                    method.common.ci = method.common.ci,
                    method.random.ci = method.random.ci,
                    adhoc.hakn.ci = adhoc.hakn.ci,
                    ##
                    level.predict = level.predict,
                    method.predict = method.predict,
                    adhoc.hakn.pi = adhoc.hakn.pi,
                    seed.predict = seed.predict,
                    ##
                    subgroup = grplab,
                    subgroup.name = "grp",
                    print.subgroup.name = FALSE,
                    test.subgroup = test.subgroup,
                    prediction.subgroup = prediction.subgroup,
                    seed.predict.subgroup = seed.predict.subgroup,
                    ##
                    backtransf = backtransf,
                    ##
                    text.common = text.common, text.random = text.random,
                    text.predict = text.predict,
                    text.w.common = text.w.common,
                    text.w.random = text.w.random,
                    ##
                    title = title,
                    complab = complab, outclab = outclab,
                    #
                    label.e = label.e, label.c = label.c,
                    label.left = label.left, label.right = label.right,
                    col.label.left = col.label.left,
                    col.label.right = col.label.right,
                    #
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
                      studlab = studlab,
                      data = x[sel, varnames],
                      ##
                      method = method, sm = sm, level = level,
                      ##
                      common = common, random = random,
                      overall = overall,
                      prediction = prediction,
                      ##
                      method.tau = method.tau, method.tau.ci = method.tau.ci,
                      level.hetstat = level.hetstat,
                      tau.common = tau.common,
                      #
                      method.I2 = method.I2,
                      #
                      level.ma = level.ma,
                      method.common.ci = method.common.ci,
                      method.random.ci = method.random.ci,
                      adhoc.hakn.ci = adhoc.hakn.ci,
                      ##
                      level.predict = level.predict,
                      method.predict = method.predict,
                      adhoc.hakn.pi = adhoc.hakn.pi,
                      seed.predict = seed.predict,
                      ##
                      backtransf = backtransf,
                      ##
                      text.common = text.common, text.random = text.random,
                      text.predict = text.predict,
                      text.w.common = text.w.common,
                      text.w.random = text.w.random,
                      ##
                      title = title,
                      complab = complab, outclab = outclab,
                      #
                      label.e = label.e, label.c = label.c,
                      label.left = label.left, label.right = label.right,
                      col.label.left = col.label.left,
                      col.label.right = col.label.right,
                      #
                      RR.Cochrane = TRUE,
                      Q.Cochrane = Q.Cochrane,
                      ##
                      warn = warn, keepdata = keepdata)
      else
        m1 <- metabin(event.e, n.e, event.c, n.c, studlab = studlab,
                      data = x[sel, varnames],
                      ##
                      method = method, sm = sm, level = level,
                      ##
                      common = common, random = random,
                      overall = overall,
                      prediction = prediction,
                      ##
                      method.tau = method.tau, method.tau.ci = method.tau.ci,
                      level.hetstat = level.hetstat,
                      tau.common = tau.common,
                      #
                      method.I2 = method.I2,
                      #
                      level.ma = level.ma,
                      method.common.ci = method.common.ci,
                      method.random.ci = method.random.ci,
                      adhoc.hakn.ci = adhoc.hakn.ci,
                      ##
                      level.predict = level.predict,
                      method.predict = method.predict,
                      adhoc.hakn.pi = adhoc.hakn.pi,
                      seed.predict = seed.predict,
                      ##
                      backtransf = backtransf,
                      ##
                      text.common = text.common, text.random = text.random,
                      text.predict = text.predict,
                      text.w.common = text.w.common,
                      text.w.random = text.w.random,
                      ##
                      title = title,
                      complab = complab, outclab = outclab,
                      #
                      label.e = label.e, label.c = label.c,
                      label.left = label.left, label.right = label.right,
                      col.label.left = col.label.left,
                      col.label.right = col.label.right,
                      #
                      RR.Cochrane = TRUE,
                      Q.Cochrane = Q.Cochrane,
                      warn = warn, keepdata = keepdata)
    }
    ##
    if (type == "C")
      m1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, studlab = studlab,
                     data = x[sel, varnames],
                     ##
                     sm = sm, level = level,
                     ##
                     common = common, random = random,
                     overall = overall,
                     prediction = prediction,
                     ##
                     method.tau = method.tau, method.tau.ci = method.tau.ci,
                     level.hetstat = level.hetstat,
                     tau.common = tau.common,
                     #
                     method.I2 = method.I2,
                     #
                     level.ma = level.ma,
                     method.common.ci = method.common.ci,
                     method.random.ci = method.random.ci,
                     adhoc.hakn.ci = adhoc.hakn.ci,
                     ##
                     level.predict = level.predict,
                     method.predict = method.predict,
                     adhoc.hakn.pi = adhoc.hakn.pi,
                     seed.predict = seed.predict,
                     ##
                     text.common = text.common, text.random = text.random,
                     text.predict = text.predict,
                     text.w.common = text.w.common,
                     text.w.random = text.w.random,
                     ##
                     title = title,
                     complab = complab, outclab = outclab,
                     #
                     label.e = label.e, label.c = label.c,
                     label.left = label.left, label.right = label.right,
                     col.label.left = col.label.left,
                     col.label.right = col.label.right,
                     #
                     keepdata = keepdata)
    ##
    if (type == "P")
      m1 <- metagen(O.E / V, sqrt(1 / V), studlab = studlab,
                    data = x[sel, varnames],
                    ##
                    sm = sm, level = level,
                    ##
                    common = common, random = random,
                    overall = overall,
                    prediction = prediction,
                    ##
                    method.tau = method.tau, method.tau.ci = method.tau.ci,
                    level.hetstat = level.hetstat,
                    tau.common = tau.common,
                    #
                    method.I2 = method.I2,
                    #
                    level.ma = level.ma,
                    method.common.ci = method.common.ci,
                    method.random.ci = method.random.ci,
                    adhoc.hakn.ci = adhoc.hakn.ci,
                    ##
                    level.predict = level.predict,
                    method.predict = method.predict,
                    adhoc.hakn.pi = adhoc.hakn.pi,
                    seed.predict = seed.predict,
                    ##
                    backtransf = backtransf,
                    ##
                    text.common = text.common, text.random = text.random,
                    text.predict = text.predict,
                    text.w.common = text.w.common,
                    text.w.random = text.w.random,
                    ##
                    title = title,
                    complab = complab, outclab = outclab,
                    #
                    label.e = label.e, label.c = label.c,
                    label.left = label.left, label.right = label.right,
                    col.label.left = col.label.left,
                    col.label.right = col.label.right,
                    #
                    keepdata = keepdata)
    ##
    if (type == "I" & method != "Peto")
      m1 <- metagen(TE, seTE, studlab = studlab,
                    data = x[sel, varnames],
                    ##
                    sm = sm, level = level,
                    ##
                    common = common, random = random,
                    overall = overall,
                    prediction = prediction,
                    ##
                    method.tau = method.tau, method.tau.ci = method.tau.ci,
                    level.hetstat = level.hetstat,
                    tau.common = tau.common,
                    #
                    method.I2 = method.I2,
                    #
                    level.ma = level.ma,
                    method.common.ci = method.common.ci,
                    method.random.ci = method.random.ci,
                    adhoc.hakn.ci = adhoc.hakn.ci,
                    ##
                    level.predict = level.predict,
                    method.predict = method.predict,
                    adhoc.hakn.pi = adhoc.hakn.pi,
                    seed.predict = seed.predict,
                    ##
                    n.e = n.e,
                    n.c = n.c,
                    ##
                    backtransf = backtransf,
                    ##
                    text.common = text.common, text.random = text.random,
                    text.predict = text.predict,
                    text.w.common = text.w.common,
                    text.w.random = text.w.random,
                    ##
                    title = title,
                    complab = complab, outclab = outclab,
                    #
                    label.e = label.e, label.c = label.c,
                    label.left = label.left, label.right = label.right,
                    col.label.left = col.label.left,
                    col.label.right = col.label.right,
                    #
                    keepdata = keepdata)
    ##
    if (type == "I" & method == "Peto")
      m1 <- metagen(O.E / V, sqrt(1 / V), studlab = studlab,
                    data = x[sel, varnames],
                    ##
                    sm = sm, level = level,
                    ##
                    common = common, random = random,
                    overall = overall,
                    prediction = prediction,
                    ##
                    method.tau = method.tau, method.tau.ci = method.tau.ci,
                    level.hetstat = level.hetstat,
                    tau.common = tau.common,
                    #
                    method.I2 = method.I2,
                    #
                    level.ma = level.ma,
                    method.common.ci = method.common.ci,
                    method.random.ci = method.random.ci,
                    adhoc.hakn.ci = adhoc.hakn.ci,
                    ##
                    level.predict = level.predict,
                    method.predict = method.predict,
                    adhoc.hakn.pi = adhoc.hakn.pi,
                    seed.predict = seed.predict,
                    ##
                    backtransf = backtransf,
                    ##
                    text.common = text.common, text.random = text.random,
                    text.predict = text.predict,
                    text.w.common = text.w.common,
                    text.w.random = text.w.random,
                    ##
                    title = title,
                    complab = complab, outclab = outclab,
                    #
                    label.e = label.e, label.c = label.c,
                    label.left = label.left, label.right = label.right,
                    col.label.left = col.label.left,
                    col.label.right = col.label.right,
                    #
                    keepdata = keepdata)
  }
  
  if (sm == "OTHER") {
    warning('Meta-analysis not possible for sm = "OTHER".')
    return(NULL)
  }


  ##
  ##
  ## (5) Add risk of bias assessment
  ##
  ##
  if (rob) {
    rd <- m1$data
    m1 <- rob(
      studlab = rd$studlab,
      item1 = if (n.domains >= 1) rd$D1 else NULL,
      item2 = if (n.domains >= 2) rd$D2 else NULL,
      item3 = if (n.domains >= 3) rd$D3 else NULL,
      item4 = if (n.domains >= 4) rd$D4 else NULL,
      item5 = if (n.domains >= 5) rd$D5 else NULL,
      item6 = if (n.domains >= 6) rd$D6 else NULL,
      item7 = if (n.domains >= 7) rd$D7 else NULL,
      item8 = if (n.domains >= 8) rd$D8 else NULL,
      item9 = if (n.domains >= 9) rd$D9 else NULL,
      item10 = if (n.domains >= 10) rd$D10 else NULL,
      weight = if (m1$random) m1$w.random else m1$w.common,
      ##
      tool = tool,
      categories = categories,
      col = col,
      symbols = symbols,
      domains = domains,
      ##
      data = m1, warn = FALSE)
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
