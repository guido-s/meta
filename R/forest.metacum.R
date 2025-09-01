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
#' @param overall A logical indicating whether overall results should be
#'   shown.
#' @param just.addcols Justification of text for additional columns
#'   (possible values: "left", "right", "center").
#' @param smlab A label for the summary measure (printed at top of
#'   figure).
#' @param type A character string or vector specifying how to
#'   plot treatment effects and confidence intervals for cumulative
#'   meta-analysis results.
#' @param layout A character string specifying the layout of the
#'   forest plot (see \code{\link{forest.meta}}).
#' @param lab.NA A character string to label missing values.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in forest plots. If \code{backtransf = TRUE},
#'   results for \code{sm = "OR"} are presented as odds ratios rather
#'   than log odds ratios, for example.
#' @param big.mark A character used as thousands separator.
#' @param digits Minimal number of significant digits for treatment
#'   effects, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for
#'   p-values.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance.
#' @param digits.tau Minimal number of significant digits for square
#'   root of between-study variance.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   statistic.
#' @param digits.cid Minimal number of significant digits for
#'   CID / decision thresholds, see \code{print.default}.
#' @param digits.percent Minimal number of significant digits for
#'   probabilities, printed as percentages, see \code{print.default}.
#' @param col The colour for cumulative meta-analysis results (only considered
#'   if \code{type = "square"}).
#' @param col.bg The background colour for squares and diamonds of
#'   cumulative meta-analysis results.
#' @param col.border The colour for the outer lines of squares and diamonds of
#'   cumulative meta-analysis results.
#' @param col.bg.predict The background colour for prediction intervals of
#'   cumulative meta-analysis results.
#' @param col.border.predict The colour for the outer lines of prediction
#'   intervals of cumulative meta-analysis results.
#' @param addrows.below.overall A numeric value indicating how many
#'   empty rows are printed between meta-analysis results and
#'   meta-analysis details.
#' @param details A logical specifying whether details on statistical
#'   methods should be printed.
#' @param \dots Additional graphical arguments (passed on to
#'   \code{\link{forest.meta}}).
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
                           #
                           leftcols = NULL, leftlabs = NULL,
                           rightcols = NULL, rightlabs = NULL,
                           #
                           prediction = x$prediction,
                           overall = x$overall,
                           just.addcols = "right",
                           smlab = "Cumulative Meta-Analysis",
                           type = "square",
                           layout = gs("layout"),
                           lab.NA = ".",
                           #
                           backtransf = x$backtransf,
                           #
                           big.mark = gs("big.mark"),
                           digits = gs("digits.forest"),
                           digits.pval = gs("digits.pval"),
                           digits.tau2 = gs("digits.tau2"),
                           digits.tau = gs("digits.tau"),
                           digits.I2 = gs("digits.I2"),
                           digits.cid = gs("digits.cid"),
                           digits.percent = 1,
                           #
                           col = gs("col.study"),
                           col.bg = 
                             ifelse(type == "diamond",
                                  gs("col.diamond"), gs("col.square")),
                           col.border =
                             ifelse(type == "diamond",
                                    gs("col.diamond.lines"),
                                    gs("col.square.lines")),
                           col.bg.predict = gs("col.predict"),
                           col.border.predict = gs("col.predict.lines"),
                           #
                           addrows.below.overall = 1L * details,
                           details = gs("forest.details"),
                           ...) {
  
  #
  #
  # (1) Check and set arguments
  #
  #
  
  chkclass(x, c("metacum", "metainf"))
  x <- updateversion(x)
  #
  type <- setchar(type, c("square", "diamond", "circle", "squarediamond"))
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
  chklogical(overall)
  #
  layout <- setchar(layout, c("meta", "BMJ", "RevMan5", "JAMA"))
  #
  chklogical(backtransf)
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.pval, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  chknumeric(digits.cid, min = 0, length = 1)
  chknumeric(digits.percent, min = 0, length = 1)
  chknumeric(addrows.below.overall, length = 1)
  #
  chklogical(details)
  #
  missing.leftcols <- missing(leftcols)
  missing.rightcols <- missing(rightcols)
  #
  missing.col.bg <- missing(col.bg)
  missing.col.border <- missing(col.border)
  
  
  avail.prop.cid.below.null <-
    !is.null(x$prop.cid.below.null) && !(all(is.na(x$prop.cid.below.null)))
  avail.prop.cid.above.null <-
    !is.null(x$prop.cid.above.null) && !(all(is.na(x$prop.cid.above.null)))
  #
  avail.prop.cid <- avail.prop.cid.below.null | avail.prop.cid.above.null
  #
  pvalNA <- all(is.na(x$pval))
  #
  if (is.null(leftcols))
    leftcols <- "studlab"
  #
  if (is.null(leftlabs))
    leftlabs <- rep(NA, length(leftcols))
  #
  if (is.null(rightcols)) {
    rightcols <- c("effect", "ci", if (!pvalNA) "pval", "tau2", "tau", "I2")
    #
    if (avail.prop.cid.below.null)
      rightcols <- c(rightcols, "prop.cid.below.null")
    #
    if (avail.prop.cid.above.null)
      rightcols <- c(rightcols, "prop.cid.above.null")
  }
  #
  if (is.null(rightlabs))
    rightlabs <- rep(NA, length(rightcols))
  #
  print.tau2 <- any(c("tau2" %in% leftcols, "tau2" %in% rightcols))
  print.tau <- any(c("tau" %in% leftcols, "tau" %in% rightcols))
  print.I2 <- any(c("I2" %in% leftcols, "I2" %in% rightcols))
  #
  print.cid.below.null <- any(c("prop.cid.below.null" %in% leftcols,
                                "prop.cid.below.null" %in% rightcols))
  #
  print.cid.above.null <- any(c("prop.cid.above.null" %in% leftcols,
                                "prop.cid.above.null" %in% rightcols))
  #
  pval <- formatPT(x$pval, digits = digits.pval, lab.NA = lab.NA)
  tau2 <- formatPT(x$tau2, digits = digits.tau2, lab.NA = lab.NA)
  tau <- formatPT(x$tau, digits = digits.tau2, lab.NA = lab.NA)
  I2 <- ifelse(is.na(x$I2), lab.NA,
               paste0(formatPT(100 * x$I2, digits = digits.I2,
                               lab.NA = lab.NA), "%"))
  #
  if (avail.prop.cid.below.null) {
    x$prop.cid.below.null <-
      ifelse(is.na(x$prop.cid.below.null), lab.NA,
             paste0(formatPT(100 * x$prop.cid.below.null,
                             digits = digits.percent), "%"))
  }
  #
  if (avail.prop.cid.above.null) {
    x$prop.cid.above.null <-
      ifelse(is.na(x$prop.cid.above.null), lab.NA,
             paste0(formatPT(100 * x$prop.cid.above.null,
                             digits = digits.percent), "%"))
  }
  #
  x.tmp <- x
  x.tmp$prediction <- any(prediction)
  class(x.tmp) <-
    c(class(x.tmp), if (inherits(x, "metainf")) "metainf" else "metacum")
  #
  text.details <-
    catmeth(x.tmp,
            x$common, x$random, any(prediction), overall, TRUE,
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
  if (avail.prop.cid)
    svd <- x$small.values == "desirable"
  #
  if (avail.prop.cid.below.null) {
    text.details <-
      paste0(text.details,
             paste0("\n- Lower decision threshold (",
                    if (svd) "beneficial " else "harmful ",
                    "effects): ",
                    formatN(x$cid.below.null, digits = digits.cid,
                            big.mark = big.mark)))
  }
  #
  if (avail.prop.cid.above.null) {
    text.details <-
      paste0(text.details,
             paste0("\n- Upper decision threshold (",
                    if (svd) "harmful " else "beneficial ",
                    "effects): ",
                    formatN(x$cid.above.null, digits = digits.cid,
                            big.mark = big.mark)))
  }
  
  
  #
  # Print prediction intervals in separate rows
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
    if (print.cid.below.null)
      prop.cid.below.null <- as.vector(
        matrix(c(x$prop.cid.below.null,
                 rep("", k.all)),
               ncol = k.all, byrow = TRUE))
    #
    if (print.cid.above.null)
      prop.cid.above.null <- as.vector(
        matrix(c(x$prop.cid.above.null,
                 rep("", k.all)),
               ncol = k.all, byrow = TRUE))
    #
    if (!is.null(col.bg))
      col.bg <- as.vector(
        matrix(c(rep(col.bg, k.all),
                 rep(col.bg.predict, k.all)),
               ncol = k.all, byrow = TRUE))
    #
    if (!is.null(col.border))
      col.border <- as.vector(
        matrix(c(rep(col.border, k.all),
                 rep(col.border.predict, k.all)),
               ncol = k.all, byrow = TRUE))
    #
    sel.pred <- as.vector(
      matrix(c(rep(TRUE, k.all), prediction),
             ncol = k.all, byrow = TRUE))
    #
    m <- metagen(TE[sel.pred], seTE[sel.pred], studlab = studlab[sel.pred],
                 common = common, random = !common,
                 prediction = any(prediction),
                 overall = overall,
                 sm = x$sm,
                 backtransf = backtransf,
                 func.backtransf = x$func.backtransf,
                 #
                 null.effect = x$null.effect,
                 #
                 label.left = x$label.left,
                 label.right = x$label.right,
                 #
                 method.tau = "DL")
    #
    m$lower <- lower[sel.pred]
    m$upper <- upper[sel.pred]
    #
    m$pval <- as.vector(
      matrix(c(pval, rep("", k.all)),
             ncol = k.all, byrow = TRUE))[sel.pred]
    #
    m$I2 <- as.vector(
      matrix(c(I2, rep("", k.all)),
             ncol = k.all, byrow = TRUE))[sel.pred]
    m$tau2 <- as.vector(
      matrix(c(tau2, rep("", k.all)),
             ncol = k.all, byrow = TRUE))[sel.pred]
    m$tau <- as.vector(
      matrix(c(tau, rep("", k.all)),
             ncol = k.all, byrow = TRUE))[sel.pred]
    #
    if (print.cid.below.null)
      m$prop.cid.below.null <- prop.cid.below.null[sel.pred]
    #
    if (print.cid.above.null)
      m$prop.cid.above.null <- prop.cid.above.null[sel.pred]
    #
    type.study <- rep(c(type, "predict"), k.all)[sel.pred]
    #
    if (!is.null(col.bg))
      col.bg <- col.bg[sel.pred]
    if (!is.null(col.border))
      col.border <- col.border[sel.pred]
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
                 overall = overall,
                 sm = x$sm,
                 backtransf = backtransf,
                 func.backtransf = x$func.backtransf,
                 #
                 null.effect = x$null.effect,
                 #
                 label.left = x$label.left,
                 label.right = x$label.right,
                 #
                 method.tau = "DL")
    #
    m$lower <- lower
    m$upper <- upper
    #
    m$pval <- pval
    m$I2 <- I2
    m$tau2 <- tau2
    m$tau <- tau
    #
    if (print.cid.below.null)
      m$prop.cid.below.null <- x$prop.cid.below.null
    #
    if (print.cid.below.null)
      m$prop.cid.above.null <- x$prop.cid.above.null
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
  # Set column labels for decision threshold probabilites
  #
  if (print.cid.below.null) {
    sel.left <- leftcols == "prop.cid.below.null"
    #
    if (any(sel.left) && is.na(leftlabs[sel.left]))
      leftlabs[sel.left] <-
        paste0("P(",
               if (x$small.values == "desirable") "benefit" else "harm",
               ")")
    #
    sel.right <- rightcols == "prop.cid.below.null"
    #
    if (any(sel.right) && is.na(rightlabs[sel.right]))
      rightlabs[sel.right] <-
      paste0("P(",
             if (x$small.values == "desirable") "benefit" else "harm",
             ")")
  }
  #
  if (print.cid.above.null) {
    sel.left <- leftcols == "prop.cid.above.null"
    #
    if (any(sel.left) && is.na(leftlabs[sel.left]))
      leftlabs[sel.left] <-
        paste0("P(",
               if (x$small.values == "desirable") "harm" else "benefit",
               ")")
    #
    sel.right <- rightcols == "prop.cid.above.null"
    #
    if (any(sel.right) && is.na(rightlabs[sel.right]))
      rightlabs[sel.right] <-
      paste0("P(",
             if (x$small.values == "desirable") "harm" else "benefit",
             ")")
  }
  #
  m$.text.details.methods <- text.details
  #
  # Move columns to left side of forest plot for JAMA and RevMan5 layouts
  #
  if (missing.leftcols & missing.rightcols & layout %in% c("JAMA", "RevMan5")) {
    leftcols <- c(leftcols, rightcols[-(1:2)], rightcols[1:2])
    rightcols <- NULL
    leftlabs <- c(leftlabs, rightlabs[-(1:2)], rightlabs[1:2])
    rightlabs <- NULL
  }
  #
  # Use default colours for JAMA and RevMan5 layouts
  #
  if (missing.col.bg & layout %in% c("JAMA", "RevMan5"))
    col.bg <- NULL
  #
  if (missing.col.border & layout %in% c("JAMA", "RevMan5"))
    col.border <- NULL
  
  
  data.p <-
    data.frame(pval = formatPT(x$pval.pooled, digits = digits.pval,
                               lab.NA = lab.NA),
               tau2 = formatPT(x$tau2.pooled, digits = digits.tau2,
                               lab.NA = lab.NA),
               tau = formatPT(x$tau.pooled, digits = digits.tau,
                              lab.NA = lab.NA),
               I2 = ifelse(is.na(x$I2.pooled), lab.NA,
                           paste0(formatPT(100 * x$I2.pooled,
                                           digits = digits.I2,
                                           lab.NA = lab.NA), "%")))
  #
  if (avail.prop.cid.below.null) {
    data.p$prop.cid.below.null <-
      ifelse(is.na(x$prop.cid.below.null.pooled), lab.NA,
             paste0(formatPT(100 * x$prop.cid.below.null.pooled,
                             digits = digits.percent), "%"))
  }
  #
  if (avail.prop.cid.above.null) {
    data.p$prop.cid.above.null <-
      ifelse(is.na(x$prop.cid.above.null.pooled), lab.NA,
             paste0(formatPT(100 * x$prop.cid.above.null.pooled,
                             digits = digits.percent), "%"))
  }
  
  
  class(m) <- c(class(m), class(x))
  m$classes <- x$classes
  #
  dots_list <- drop_from_dots(list(...),
                              c("col.study", "col.square", "col.square.lines",
                                "overall.hetstat", "overall.hetstat",
                                "data.pooled"),
                              c("col", "col.bg", "col.border",
                                "col.bg.predict", "col.border.predict",
                                ""))
  #
  args_list <-
    list(x = m,
         leftcols = leftcols, leftlabs = leftlabs,
         rightcols = rightcols, rightlabs = rightlabs,
         overall.hetstat = FALSE,
         type.study = type.study,
         weight.study = "same",
         lab.NA = lab.NA, smlab = smlab,
         data.pooled = data.p,
         just.addcols = just.addcols,
         #
         layout = layout,
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
  #
  invisible(res)
}
