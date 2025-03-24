#' Forest plot to display the result of a cumulative meta-analysis
#' 
#' @description
#' Draws a forest plot in the active graphics window (using grid
#' graphics system).
#' 
#' @aliases forest.metacum
#' 
#' @param x An object of class \code{\link{metacum}}.
#' @param leftcols A character vector specifying (additional) columns
#'   to be plotted on the left side of the forest plot or a logical
#'   value.
#' @param leftlabs A character vector specifying labels for
#'   (additional) columns on left side of the forest plot.
#' @param rightcols A character vector specifying (additional) columns
#'   to be plotted on the right side of the forest plot or a logical
#'   value.
#' @param rightlabs A character vector specifying labels for
#'   (additional) columns on right side of the forest plot.
#' @param prediction A logical indicating whether prediction
#'   intervals should be printed.
#' @param just.addcols Justification of text for additional columns
#'   (possible values: "left", "right", "center").
#' @param smlab A label for the summary measure (printed at top of
#'   figure).
#' @param type A character string or vector specifying how to
#'   plot treatment effects and confidence intervals for cumulative
#'   meta-analysis results.
#' @param lab.NA A character string to label missing values.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in forest plots. If \code{backtransf = TRUE},
#'   results for \code{sm = "OR"} are presented as odds ratios rather
#'   than log odds ratios, for example.
#' @param big.mark A character used as thousands separator.
#' @param digits Minimal number of significant digits for treatment
#'   effects, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance.
#' @param digits.tau Minimal number of significant digits for square
#'   root of between-study variance.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   statistic.
#' @param col The colour for cumulative meta-analysis results (only considered
#'   if \code{type = "square"}).
#' @param col.bg The background colour for squares, diamonds and prediction
#'   intervals of cumulative meta-analysis results.
#' @param col.border The colour for the outer lines of squares, diamonds and
#'   prediction intervals of cumulative meta-analysis results.
#' @param addrows.below.overall A numeric value indicating how many
#'   empty rows are printed between meta-analysis results and
#'   meta-analysis details.
#' @param details A logical specifying whether details on statistical
#'   methods should be printed.
#' @param \dots Additional graphical arguments (passed on to
#'   \code{\link{forest.meta}}).
#' 
#' 
#' 
#' @details
#' A forest plot, also called confidence interval plot, is drawn in
#' the active graphics window. The forest functions in R package
#' \bold{meta} are based on the grid graphics system. In order to
#' print the forest plot, resize the graphics window and either use
#' \code{\link{dev.copy2eps}} or \code{\link{dev.copy2pdf}}. Another
#' possibility is to create a file using \code{\link{pdf}},
#' \code{\link{png}}, or \code{\link{svg}} and to specify the width
#' and height of the graphic (see \code{\link{forest.meta}} examples).
#' 
#' The arguments \code{leftcols} and \code{rightcols} can be used to
#' specify columns which are plotted on the left and right side of the
#' forest plot, respectively.
#' 
#' The arguments \code{leftlabs} and \code{rightlabs} can be used to
#' specify column headings which are plotted on left and right side of
#' the forest plot, respectively. For certain columns predefined
#' labels exist. For other columns, the column name will be used as a
#' label. It is possible to only provide labels for new columns (see
#' \code{\link{forest.meta}} examples). Otherwise the length of
#' \code{leftlabs} and \code{rightlabs} must be the same as the number
#' of printed columns, respectively. The value \code{NA} can be used
#' to specify columns which should use default labels.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{forest.meta}}, \code{\link{metacum}},
#'   \code{\link{settings.meta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' data(Fleiss1993bin)
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac,
#'   data = Fleiss1993bin, studlab = study, sm = "RR", method = "I")
#' m1
#' metacum(m1)
#' metacum(m1, pooled = "random")
#' 
#' forest(metacum(m1))
#' forest(metacum(m1, pooled = "random"))
#' forest(metacum(m1, pooled = "random", prediction = TRUE))
#'
#' @method forest metacum
#' @export


forest.metacum <- function(x,
                           leftcols = "studlab",
                           leftlabs = rep(NA, length(leftcols)),
                           rightcols =
                             c("effect", "ci", "pval", "tau2", "tau", "I2"),
                           rightlabs = rep(NA, length(rightcols)),
                           prediction = x$prediction,
                           just.addcols = "right",
                           smlab = "Cumulative Meta-Analysis",
                           type = "square",
                           lab.NA = ".",
                           #
                           backtransf = x$backtransf,
                           #
                           big.mark = gs("big.mark"),
                           digits = gs("digits.forest"),
                           digits.tau2 = gs("digits.tau2"),
                           digits.tau = gs("digits.tau"),
                           digits.I2 = gs("digits.I2"),
                           #
                           col = gs("col.study"),
                           col.bg = gs("col.square"),
                           col.border =
                             if (type == "square") gs("col.square.lines")
                             else "black",
                           #
                           addrows.below.overall = 1L * details,
                           details = gs("forest.details"),
                           ...) {
  
  #
  #
  # (1) Check and set arguments
  #
  #
  chkclass(x, "metacum")
  x <- updateversion(x)
  #
  type <-
    setchar(type, c("square", "diamond", "predict", "circle", "squarediamond"))
  #
  just.addcols <- setchar(just.addcols, c("left", "center", "right"))
  #
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  common <- x$pooled == "common"
  #
  if (!missing(prediction))
    prediction <- catch("prediction", mc, x, sfsp)
  #
  k.all <- length(x$TE)
  #
  if (length(prediction) > 1 && length(prediction) != k.all)
    stop("Argument 'prediction' must be of length 1 or number of studies.",
         call. = FALSE)
  #
  if (!is.logical(prediction))
    stop("Argument 'prediction' must be of type logical.",
         call. = FALSE)
  #
  if (length(prediction) == 1)
    prediction <- rep(prediction, k.all)
  #
  prediction <- prediction & !common
  #
  chklogical(backtransf)
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  chknumeric(addrows.below.overall, length = 1)
  #
  chklogical(details)
  #
  print.tau2 <- any(c("tau2" %in% leftcols, "tau2" %in% rightcols))
  print.tau <- any(c("tau" %in% leftcols, "tau" %in% rightcols))
  print.I2 <- any(c("I2" %in% leftcols, "I2" %in% rightcols))
  #
  x.tmp <- x
  x.tmp$prediction <- any(prediction)
  class(x.tmp) <- c(class(x.tmp), "metacum")
  #
  text.details <-
    catmeth(x.tmp,
            x$pooled == "common", x$pooled == "random", any(prediction),
            TRUE, FALSE,
            #
            func.transf = x$func.transf,
            backtransf = backtransf, func.backtransf = x$func.backtransf,
            #
            big.mark = big.mark, digits = digits,
            digits.tau = digits.tau,
            text.tau = gs("text.tau"), text.tau2 = gs("text.tau2"),
            #
            print.tau2 = print.tau2 | x$pooled == "random",
            print.tau2.ci = FALSE,
            print.tau = print.tau | x$pooled == "random",
            print.tau.ci = FALSE,
            #
            print.I2 = print.I2, text.I2 = gs("text.I2"),
            #
            print.df = TRUE, prediction.subgroup = FALSE,
            #
            forest = TRUE)
  #
  if (any(prediction)) {
    TE <- as.vector(
      matrix(c(x$TE, rep(NA, k.all)),
             ncol = k.all, byrow = TRUE))
    #
    seTE <- as.vector(
      matrix(c(x$seTE, rep(NA, k.all)),
             ncol = k.all, byrow = TRUE))
    #
    lower <- as.vector(
      matrix(c(x$lower, x$lower.predict),
               ncol = k.all, byrow = TRUE))
    #
    upper <- as.vector(
      matrix(c(x$upper,
               x$upper.predict),
               ncol = k.all, byrow = TRUE))
    #
    studlab <- as.vector(
      matrix(c(x$studlab,
               rep("", k.all)),
               ncol = k.all, byrow = TRUE))
    #
    sel.pred <- as.vector(
      matrix(c(rep(TRUE, k.all), prediction),
             ncol = k.all, byrow = TRUE))
    #
    m <- metagen(TE[sel.pred], seTE[sel.pred], studlab = studlab[sel.pred],
                 common = common, random = !common,
                 prediction = any(prediction),
                 sm = x$sm,
                 backtransf = backtransf,
                 func.backtransf = x$func.backtransf)
    #
    m$lower <- lower[sel.pred]
    m$upper <- upper[sel.pred]
    #
    m$pval <- as.vector(
      matrix(c(x$pval, rep(NA, k.all)),
             ncol = k.all, byrow = TRUE))[sel.pred]
    #
    m$I2 <- as.vector(
      matrix(c(x$I2, rep(NA, k.all)),
             ncol = k.all, byrow = TRUE))[sel.pred]
    m$tau2 <- as.vector(
      matrix(c(x$tau2, rep(NA, k.all)),
             ncol = k.all, byrow = TRUE))[sel.pred]
    m$tau <- as.vector(
      matrix(c(x$tau, rep(NA, k.all)),
             ncol = k.all, byrow = TRUE))[sel.pred]
    #
    type.study <- rep(c(type, "predict"), k.all)[sel.pred]
  }
  #
  else {
    TE <- x$TE
    seTE <- x$seTE
    #
    lower <- x$lower
    upper <- x$upper
    #
    studlab <- x$studlab
    #
    m <- metagen(TE, seTE, studlab = studlab,
                 common = common, random = !common,
                 prediction = any(prediction),
                 sm = x$sm,
                 backtransf = backtransf,
                 func.backtransf = x$func.backtransf)
    #
    m$lower <- lower
    m$upper <- upper
    #
    m$I2 <- x$I2
    m$tau2 <- x$tau2
    m$tau <- x$tau
    #
    type.study <- type
  }
  #
  m$TE.common <- x$TE.pooled
  m$lower.common <- x$lower.pooled
  m$upper.common <- x$upper.pooled
  m$statistic.common <- x$statistic.pooled
  m$pval.common <- x$pval.pooled
  #
  m$TE.random <- x$TE.pooled
  m$lower.random <- x$lower.pooled
  m$upper.random <- x$upper.pooled
  m$statistic.random <- x$statistic.pooled
  m$pval.random <- x$pval.pooled
  #
  m$df.random <- x$df.random.pooled
  #
  if (any(prediction)) {
    m$lower.predict <- x$lower.predict.pooled
    m$upper.predict <- x$upper.predict.pooled
    m$method.predict <- x$method.predict
  }
  #
  m$level <- x$level.ma <- x$level.ma
  m$level.predict <- x$level.predict
  #
  m$method <- x$method
  m$method.random <- x$method.random
  #
  m$method.random.ci <- x$method.random.ci
  m$adhoc.hakn.ci <- x$adhoc.hakn.ci
  #
  m$method.tau <- x$method.tau
  m$method.tau.ci <- x$method.tau.ci
  #
  m$tau.preset <- x$tau.preset
  m$TE.tau <- x$TE.tau
  #
  m$method.I2 <- x$method.I2
  #
  m$k <- x$k.pooled
  m$k.study <- x$k.study.pooled
  m$k.all <- x$k.all.pooled
  m$k <- x$k.TE.pooled
  #
  if (any(rightcols %in% c("ci", "effect.ci")) |
      any(leftcols %in% c("ci", "effect.ci"))) {
    level.ma <- x$level.ma
    level.predict <- x$level.predict
    #
    if (any(prediction)) {
      if (level.ma == level.predict)
        ci.lab <- paste0(100 * level.ma, "%-CI/PI")
      else
        ci.lab <-
          paste0(100 * level.ma, "%-CI / ", 100 * level.predict, "%-PI")
    }
    else
      ci.lab <- paste0(100 * level.ma, "%-CI")
    #
    sel.left <- leftcols == "ci"
    #
    if (any(sel.left) && is.na(leftlabs[sel.left]))
      leftlabs[sel.left] <- ci.lab
    #
    sel.right <- rightcols == "ci"
    #
    if (any(sel.right) && is.na(rightlabs[sel.right]))
      rightlabs[sel.right] <- ci.lab
    #
    sel.left <- leftcols == "effect.ci"
    #
    if (any(sel.left) && is.na(leftlabs[sel.left]))
      leftlabs[sel.left] <-
        paste(smlab(x$sm, backtransf, x$pscale, x$irscale), ci.lab)
    #
    sel.right <- rightcols == "effect.ci"
    #
    if (any(sel.right) && is.na(rightlabs[sel.right]))
      rightlabs[sel.right] <-
        paste(smlab(x$sm, backtransf, x$pscale, x$irscale), ci.lab)
  }
  #
  m$.text.details.methods = text.details
  #
  data.p <- data.frame(pval = x$pval.pooled,
                       tau2 = x$tau2.pooled,
                       tau = x$tau.pooled,
                       I2 = x$I2.pooled)
  #
  dots_list <- drop_from_dots(list(...),
                              c("col.study", "col.square", "col.square.lines",
                                "overall.hetstat", "overall.hetstat",
                                "data.pooled"),
                              c("col", "col.bg", "col.border",
                                "", "",
                                ""))
  #
  args_list <-
    list(x = m,
         leftcols = leftcols,
         rightcols = rightcols, rightlabs = rightlabs,
         overall.hetstat = FALSE,
         type.study = type.study,
         weight.study = "same",
         lab.NA = lab.NA, smlab = smlab,
         data.pooled = data.p,
         just.addcols = just.addcols,
         #
         backtransf = backtransf,
         #
         big.mark = big.mark,
         #
         digits = digits,
         digits.tau2 = digits.tau2,
         digits.tau = digits.tau,
         digits.I2 = digits.I2,
         #
         addrows.below.overall = addrows.below.overall,
         #
         col.study = col,
         col.square = col.bg, col.square.lines = col.border,
         #
         details = details)
  #
  res <- do.call("forest.meta", c(args_list, dots_list))
  return(invisible(res))
  res <- forest(m,
                leftcols = leftcols,
                rightcols = rightcols, rightlabs = rightlabs,
                overall.hetstat = FALSE,
                type.study = type.study,
                weight.study = "same",
                lab.NA = lab.NA, smlab = smlab,
                data.pooled = data.p,
                just.addcols = just.addcols,
                #
                backtransf = backtransf,
                #
                big.mark = big.mark,
                #
                digits = digits,
                digits.tau2 = digits.tau2,
                digits.tau = digits.tau,
                digits.I2 = digits.I2,
                #
                addrows.below.overall = addrows.below.overall,
                #
                col.study = col,
                col.square = col.bg, col.square.lines = col.border,
                #
                details = details,
                ...)
  #
  invisible(res)
}
