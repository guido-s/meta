#' Forest plot to display the result of a cumulative meta-analysis
#' 
#' @description
#' Draws a forest plot showing the results as studies are cumulatively added,
#' using the grid graphics system.
#' 
#' @aliases forest.metacum
#' 
#' @param x An object of class \code{\link{metacum}}.
#' @param layout A character string specifying the layout of the
#'   forest plot (see \code{\link{forest.meta}}).
#' @param prediction A logical indicating whether prediction
#'   intervals should be printed.
#' @param overall A logical indicating whether overall results should be
#'   shown.
#' @param leftcols A character vector specifying columns to be printed on the
#'   left side of the forest plot or a logical value; see
#'   \code{\link{forest.meta}}.
#' @param rightcols A character vector specifying columns to be printed on the
#'   right side of the forest plot or a logical value; see
#'   \code{\link{forest.meta}}.
#' @param leftlabs A character vector specifying labels for columns on the left
#'   side of the forest plot (see Details).
#' @param rightlabs A character vector specifying labels for columns on the
#'   right side of the forest plot (see Details).
#' @param lab.NA A character string to label missing values.
#' @param col A named list specifying colours used in the forest
#'   plot; see \code{\link{forest-formatting}}.
#' @param print A named list specifying which information to print in the forest
#'   plot; see \code{\link{forest-formatting}}.
#' @param type A character string or vector specifying how to
#'   plot treatment effects and confidence intervals for cumulative
#'   meta-analysis results.
#' @param digits A named list specifying digits used for displayed results; a
#'   scalar value is also accepted. See \code{\link{forest-formatting}}.
#' @param details A logical specifying whether details on statistical
#'   methods should be printed.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in forest plots. If \code{backtransf = TRUE},
#'   results for \code{sm = "OR"} are presented as odds ratios rather
#'   than log odds ratios, for example.
#' @param big.mark A character used as thousands separator.
#' @param main A title for the forest plot.
#' @param filename File name to store forest plot.
#' @param path Path of the directory to save plot to: `path` and `filename` are
#'   combined to create the fully qualified file name. Defaults to the working
#'   directory.
#' @param device Device to use. Can either be a device function (e.g.,
#'   \code{\link{png}}), or one of \code{"eps"}, \code{"ps"},
#'   \code{"tex"} (pictex), \code{"pdf"}, \code{"jpeg"}, \code{"tiff"},
#'   \code{"png"}, \code{"bmp"}, \code{"svg"} or \code{"wmf"} (Windows only).
#'   By default, the device is guessed based on the `filename` extension.
#' @param device.args List with additional graphical parameters passed on
#'   to graphics function.
#' @param units Units of the `width` and `height` of the file containing the
#'   forest plot. One of \code{"in"} (for inches, default), \code{"cm"},
#'   \code{"mm"}, or \code{"px"} (for pixels), can be abbreviated.
#' @param height Height of the graphics file.
#' @param width Width of the graphics file.
#' @param dpi Plot resolution.
#' @param dev.off A logical to specify whether current graphics device
#'   should be shut down.
#' @param new A logical value indicating whether a new figure should
#'   be printed in an existing graphics window.
#' @param \dots Additional graphical arguments (passed on to
#'   \code{\link{forest.meta}}).
#' 
#' @details
#' Each row of the forest plot represents a meta-analysis based on the first
#' \code{k} studies in the cumulative sequence. Results are shown in the order
#' in which studies were added to the cumulative meta-analysis. Common-effect
#' or random-effects estimates, prediction intervals, heterogeneity statistics,
#' and study-level information can be selected with the corresponding
#' arguments.
#'
#' A forest plot, also called a confidence interval plot, is drawn in the
#' active graphics window using the grid graphics system. To save the plot,
#' specify a file name with argument \code{filename}; the graphics device is
#' selected from the file extension by default (see \code{\link{forest.meta}}
#' for examples).
#' 
#' The named list \code{print} determines which result columns are shown in
#' the forest plot. The argument \code{layout} determines the default
#' placement of these columns on the left or right side. The columns and their
#' placement are determined internally and cannot be changed by the user.
#'
#' Arguments \code{leftcols} and \code{rightcols} can be used to specify
#' columns printed on the left or right side, respectively. Arguments
#' \code{leftlabs} and \code{rightlabs} can be used to provide headings for
#' these columns. Use \code{NA} for a column heading to retain its default
#' label.
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
#' ma <- metabin(d.asp, n.asp, d.plac, n.plac,
#'   data = Fleiss1993bin, studlab = study, sm = "RR", method = "I")
#' ma
#' metacum(ma)
#' metacum(ma, pooled = "random")
#' 
#' forest(metacum(ma))
#' forest(metacum(ma, pooled = "random"))
#' forest(metacum(ma, pooled = "random", prediction = TRUE))
#'
#' @method forest metacum
#' @export

forest.metacum <- function(x,
                           #
                           layout = gs("layout"),
                           #
                           prediction = x$prediction,
                           overall = x$overall,
                           #
                           leftcols = NULL, rightcols = NULL,
                           leftlabs = NULL, rightlabs = NULL,
                           #
                           lab.NA = gs("lab.NA"),
                           #
                           col = NULL,
                           print = NULL,
                           type = "square",
                           digits = NULL,
                           #
                           details = gs("forest.details"),
                           #
                           backtransf = x$backtransf,
                           big.mark = gs("big.mark"),
                           #
                           main = "Cumulative Meta-Analysis",
                           #
                           # Output
                           #
                           filename = NULL,
                           path = NULL,
                           device = NULL, device.args = NULL,
                           units = "in",
                           height = NULL, width = NULL,
                           dpi = 300,
                           dev.off = NULL,
                           new = TRUE,
                           #
                           ...) {
  
  #
  #
  # (1) Check and set arguments
  #
  #
  
  chkclass(x, c("metacum", "metainf"))
  x <- updateversion(x)
  #
  layout <- setchar(layout, c(gs("layout4forest"), "subgroup"), length = 1)
  #
  common <- x$pooled == "common"
  #
  k.all <- length(x$TE)
  #
  # Drop internally used arguments from '...'
  #
  drop_args <- c(
    "col.study", "col.square", "col.square.lines", "type.study",
    "weight.study", "common", "random", "overall.hetstat",
    "test.overall", "test.overall.common", "test.overall.random",
    "data.pooled")
  #
  replace_args <- list(col.study = "col", col.square = "col.bg",
                       col.square.lines = "col.border")
  #
  args <- drop_from_dots(list(...), drop_args, replace_args)
  #
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  if (!missing(prediction))
    prediction <- catch("prediction", mc, x, sfsp)
  #
  chklogical(prediction, length = c(1, k.all))
  #
  if (length(prediction) == 1)
    prediction <- rep(prediction, k.all)
  #
  prediction <- prediction & !common
  #
  chklogical(overall)
  #
  just.addcols <- "right"
  #
  type <- setchar(type, c("square", "diamond", "circle", "squarediamond"))
  #
  # Define colours
  #
  col <- as_named_list(col)
  col.list <- merge_legacy(col, args, col_defaults(x, "metacum"), "col")
  #
  col <- col.list$res
  args <- col.list$args
  #
  cols <- col.list$res
  #
  col <- cols$col
  col.bg <- cols$col.bg
  col.border <- cols$col.border
  col.bg.predict <- cols$bg.predict
  col.border.predict <- cols$col.border.predict
  #
  # Colour schemes for layouts:
  # - colors[1] - vertical line for common effect or random effects model
  # - colors[2] - diamond for meta-analysis results
  # - colors[3] - outer lines of diamonds
  # - colors[4] - square  colour
  # - colors[5] - outer lines of squares
  # - colors[6] - line color (axes, reference line, header lines)
  # - colors[7] - prediction interval
  # - colors[8] - outer lines of prediction intervals
  # - colors[9] - subgroups
  # - colors[10] - color within squares or circles
  # - colors[11] - color of circles and outer lines of circles
  #
  colors <- layout_colors(layout)
  #
  col <- replaceNULL(col, colors[1])
  #
  if (missing(col.bg) || is.null(col.bg))
    col.bg <- ifelse(type == "diamond", colors[2],
                     ifelse(type == "circle", colors[11], colors[4]))
  #
  if (missing(col.border) || is.null(col.border))
    col.border <- ifelse(type == "diamond", colors[3],
                         ifelse(type == "circle", colors[11], colors[5]))
  #
  if (missing(col.bg.predict) || is.null(col.bg.predict))
    col.bg.predict <- colors[7]
  #
  if (missing(col.border.predict) || is.null(col.border.predict))
    col.border.predict = colors[8]
  #
  # Print information
  #
  print.list <-
    merge_legacy(print, args, print_defaults(x, "metacum"), "print")
  #
  print <- print.list$res
  args <- print.list$args
  #
  print.pval <- print$pval
  print.stat <- print$stat
  print.I2 <- print$I2
  print.I2.ci <- print$I2.ci
  print.tau2 <- print$tau2
  print.tau2.ci <- print$tau2.ci
  print.tau <- print$tau
  print.tau.ci <- print$tau.ci
  print.Q <- print$Q
  print.pval.Q <- print$pval.Q
  print.cidprop <- print$cidprop
  #
  chklogical(print.pval)
  chklogical(print.stat)
  chklogical(print.I2)
  chklogical(print.I2.ci)
  chklogical(print.tau2)
  chklogical(print.tau2.ci)
  chklogical(print.tau)
  chklogical(print.tau.ci)
  chklogical(print.Q)
  chklogical(print.pval.Q)
  chklogical(print.cidprop)
  #
  print.tau2.ci <- print.tau2.ci && print.tau2
  print.tau.ci <- print.tau.ci && print.tau
  print.I2.ci <- print.I2.ci && print.I2
  #
  print.cidprop.below.null <- print.cidprop &&
    !is.null(x$prop.cid.below.null) && !(all(is.na(x$prop.cid.below.null)))
  print.cidprop.above.null <- print.cidprop &&
    !is.null(x$prop.cid.above.null) && !(all(is.na(x$prop.cid.above.null)))
  #
  print.cidprop <- print.cidprop.below.null || print.cidprop.above.null
  #
  chklogical(backtransf)
  #
  # Digits
  #
  digits <- as_named_list(digits)
  digits.list <-
    merge_legacy(digits, args, digits_defaults("metacum"), "digits")
  #
  digitss <- digits.list$res
  args <- digits.list$args
  #
  digits <- digitss$digits
  digits.pval <- digitss$pval
  digits.stat <- digitss$stat
  digits.tau2 <- digitss$tau2
  digits.tau <- digitss$tau
  digits.I2 <- digitss$I2
  digits.cid <- digitss$cid
  digits.percent <- digitss$percent
  #
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.pval, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  chknumeric(digits.cid, min = 0, length = 1)
  chknumeric(digits.percent, min = 0, length = 1)
  #
  details <- replaceNULL(details, layout == "meta")
  chklogical(details)
  #
  print.pval <- print.pval && ! all(is.na(x$pval))
  print.stat <- print.stat && !all(is.na(x$statistic))
  
  
  #
  #
  # (2) Columns and labels shown in forest plot
  #
  #
  
  # Columns and labels on the left side
  #
  if (is.null(leftcols)) {
    leftcols <- "studlab"
    #
    if (layout != "meta8") {
      if (print.tau2)
        leftcols <- c(leftcols, "tau2")
      #
      if (print.tau2.ci)
        leftcols <- c(leftcols, "tau2.ci")
      #
      if (print.tau)
        leftcols <- c(leftcols, "tau")
      #
      if (print.tau.ci)
        leftcols <- c(leftcols, "tau.ci")
      #
      if (print.I2)
        leftcols <- c(leftcols, "I2")
      #
      if (print.I2.ci)
        leftcols <- c(leftcols, "I2.ci")
      #
      if (print.Q)
        leftcols <- c(leftcols, "Q")
      #
      if (print.pval.Q)
        leftcols <- c(leftcols, "pval.Q")
      #
      if (print.cidprop.below.null)
        leftcols <- c(leftcols, "prop.cid.below.null")
      #
      if (print.cidprop.above.null)
        leftcols <- c(leftcols, "prop.cid.above.null")
      #
      if (layout == "RevMan5" || layout == "JAMA")
        leftcols <- c(leftcols, "effect.ci")
      else if (layout != "BMJ")
        leftcols <- c(leftcols, "effect", "ci")
      #
      if (print.stat && !(layout == "BMJ" || layout == "meta8"))
        leftcols <- c(leftcols, "statistic")
      #
      if (print.pval && !(layout == "BMJ" || layout == "meta8"))
        leftcols <- c(leftcols, "pval")
    }
  }
  #
  # Argument leftlabs
  #
  if (is.null(leftlabs))
    leftlabs <- rep_len(NA, length(leftcols))
  else {
    if (length(leftlabs) != length(leftcols))
      stop("Argument 'leftlabs' must provide ", length(leftcols),
           " labels (leftcols: ",
           paste(leftcols, collapse = ", "), ")",
           call. = FALSE)
  }
  #
  if (print.tau2.ci && "tau2.ci" %in% leftcols &&
      is.na(leftlabs[leftcols == "tau2.ci"]))
    leftlabs[leftcols == "tau2.ci"] <- percent_label(x$level.hetstat)
  #
  if (print.tau.ci && "tau.ci" %in% leftcols &&
      is.na(leftlabs[leftcols == "tau.ci"]))
    leftlabs[leftcols == "tau.ci"] <- percent_label(x$level.hetstat)
  #
  if (print.I2.ci && "I2.ci" %in% leftcols &&
      is.na(leftlabs[leftcols == "I2.ci"]))
    leftlabs[leftcols == "I2.ci"] <- percent_label(x$level.hetstat)
  
  
  # Columns and labels on the right side
  #
  if (is.null(leftcols)) {
    if (layout == "BMJ") {
      rightcols <- c("effect.ci")
      #
      if (print.stat)
        rightcols <- c(rightcols, "statistic")
      #
      if (print.pval)
        rightcols <- c(rightcols, "pval")
    }
    else if (layout == "meta8") {
      rightcols <- c("effect", "ci")
      #
      if (print.stat)
        rightcols <- c(rightcols, "statistic")
      #
      if (print.pval)
        rightcols <- c(rightcols, "pval")
      #
      if (print.tau2)
        rightcols <- c(rightcols, "tau2")
      #
      if (print.tau2.ci)
        leftcols <- c(leftcols, "tau2.ci")
      #
      if (print.tau)
        rightcols <- c(rightcols, "tau")
      #
      if (print.tau.ci)
        leftcols <- c(leftcols, "tau.ci")
      #
      if (print.I2)
        rightcols <- c(rightcols, "I2")
      #
      if (print.I2.ci)
        leftcols <- c(leftcols, "I2.ci")
      #
      if (print.cidprop.below.null)
        rightcols <- c(rightcols, "prop.cid.below.null")
      #
      if (print.cidprop.above.null)
        rightcols <- c(rightcols, "prop.cid.above.null")
    }
    else
      rightcols <- FALSE
  }
  #
  # Argument rightlabs
  #
  if (is.null(rightlabs))
    rightlabs <- rep_len(NA, length(rightcols))
  else {
    if (length(rightlabs) != length(rightcols))
      stop("Argument 'rightlabs' must provide ", length(rightcols),
           " labels (rightcols: ",
           paste(rightcols, collapse = ", "), ")",
           call. = FALSE)
  }
  #
  if (print.tau2.ci && "tau2.ci" %in% rightcols &&
      is.na(rightlabs[rightcols == "tau2.ci"]))
    rightlabs[rightcols == "tau2.ci"] <- percent_label(x$level.hetstat)
  #
  if (print.tau.ci && "tau.ci" %in% 
      rightcols && is.na(rightlabs[rightcols == "tau.ci"]))
    rightlabs[rightcols == "tau.ci"] <- percent_label(x$level.hetstat)
  #
  if (print.I2.ci && "I2.ci" %in% rightcols &&
      is.na(rightlabs[rightcols == "I2.ci"]))
    rightlabs[rightcols == "I2.ci"] <- percent_label(x$level.hetstat)
  
    
  pval <- formatPT(x$pval, digits = digits.pval, lab.NA = lab.NA)
  statistic <- formatN(x$statistic, digits = digits.stat, lab.NA = lab.NA,
                       big.mark = big.mark)
  tau2 <- formatPT(x$tau2, digits = digits.tau2, lab.NA = lab.NA)
  tau <- formatPT(x$tau, digits = digits.tau2, lab.NA = lab.NA)
  #
  I2 <- percent_label(x$I2, "", digits.I2, formatPT, lab.NA = lab.NA)
  #
  if (print.tau2.ci) {
    tau2.ci <-
      formatCI(x$lower.tau2, x$upper.tau2,
               digits = digits.tau2, lab.NA = lab.NA, na.lab = TRUE)
  }
  #
  if (print.tau.ci) {
    tau.ci <-
      formatCI(x$lower.tau, x$upper.tau,
               digits = digits.tau, lab.NA = lab.NA, na.lab = TRUE)
  }
  #
  if (print.I2.ci) {
    I2.ci <-
      formatCI(x$lower.I2, x$upper.I2, scale = 100, suffix = "%",
               digits = digits.I2, lab.NA = lab.NA, na.lab = TRUE)
  }
  #
  if (print.cidprop.below.null) {
    x$prop.cid.below.null <-
      percent_label(x$prop.cid.below.null, "", digits.percent, formatPT,
                    lab.NA = lab.NA)
  }
  #
  if (print.cidprop.above.null) {
    x$prop.cid.above.null <-
      percent_label(x$prop.cid.above.null, "", digits.percent, formatPT,
                    lab.NA = lab.NA)
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
            print.tau2.ci = print.tau2.ci,
            print.tau = print.tau | x$pooled == "random",
            print.tau.ci = print.tau.ci,
            #
            print.I2 = print.I2, text.I2 = gs("text.I2"),
            #
            print.df = TRUE, prediction.subgroup = FALSE,
            #
            forest = TRUE)
  #
  if (print.cidprop)
    svd <- x$small.values == "desirable"
  #
  if (print.cidprop.below.null) {
    text.details <-
      paste0(text.details,
             paste0("\n- Lower decision threshold (",
                    if (svd) "beneficial " else "harmful ",
                    "effects): ",
                    formatN(x$cid.below.null, digits = digits.cid,
                            big.mark = big.mark)))
  }
  #
  if (print.cidprop.above.null) {
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
    if (print.cidprop)
      prop.cid.below.null <- as.vector(
        matrix(c(x$prop.cid.below.null,
                 rep("", k.all)),
               ncol = k.all, byrow = TRUE))
    #
    if (print.cidprop)
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
    m <- metagen(TE[sel.pred], seTE[sel.pred],
                 lower = lower[sel.pred],
                 upper = upper[sel.pred],
                 studlab = studlab[sel.pred],
                 #
                 common = common, random = !common,
                 prediction = any(prediction),
                 overall = overall,
                 #
                 sm = x$sm,
                 transf = TRUE,
                 backtransf = backtransf,
                 func.backtransf = x$func.backtransf,
                 #
                 pscale = x$pscale, irscale = x$irscale, irunit = x$irunit,
                 #
                 level = x$level,
                 level.ma = x$level.ma,
                 level.predict = x$level.predict,
                 level.hetstat = x$level.hetstat,
                 #
                 method.I2 = x$method.I2,
                 #
                 null.effect = x$null.effect,
                 #
                 label.left = x$label.left,
                 label.right = x$label.right,
                 #
                 method.tau = "DL")
    #
    if (x$sm == "PFT") {
      m$n.harmonic.mean <- as.vector(
        matrix(c(x$n.harmonic.mean, x$n.harmonic.mean),
               ncol = k.all, byrow = TRUE))[sel.pred]
      m$n.harmonic.mean.pooled <- x$n.harmonic.mean.pooled
    }
    else if (x$sm == "IRFT") {
      m$t.harmonic.mean <- as.vector(
        matrix(c(x$t.harmonic.mean, x$t.harmonic.mean),
               ncol = k.all, byrow = TRUE))[sel.pred]
      m$t.harmonic.mean.pooled <- x$t.harmonic.mean.pooled
    }
    #
    m$pval <- as.vector(
      matrix(c(pval, rep("", k.all)),
             ncol = k.all, byrow = TRUE))[sel.pred]
    #
    m$statistic <- as.vector(
      matrix(c(statistic, rep("", k.all)),
             ncol = k.all, byrow = TRUE))[sel.pred]
    #
    m$I2 <- as.vector(matrix(c(I2, rep("", k.all)),
                             ncol = k.all, byrow = TRUE))[sel.pred]
    #
    m$tau2 <- as.vector(matrix(c(tau2, rep("", k.all)),
                               ncol = k.all, byrow = TRUE))[sel.pred]
    #
    m$tau <- as.vector(matrix(c(tau, rep("", k.all)),
                              ncol = k.all, byrow = TRUE))[sel.pred]
    #
    if (print.cidprop)
      m$prop.cid.below.null <- prop.cid.below.null[sel.pred]
    #
    if (print.cidprop)
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
    m <- metagen(TE, seTE,
                 lower = lower, upper = upper,
                 studlab = studlab,
                 #
                 common = common, random = !common,
                 prediction = any(prediction),
                 overall = overall,
                 #
                 sm = x$sm,
                 #
                 transf = TRUE,
                 backtransf = backtransf,
                 func.backtransf = x$func.backtransf,
                 #
                 null.effect = x$null.effect,
                 #
                 pscale = x$pscale, irscale = x$irscale, irunit = x$irunit,
                 #
                 level = x$level,
                 level.ma = x$level.ma,
                 level.predict = x$level.predict,
                 level.hetstat = x$level.hetstat,
                 #
                 method.I2 = x$method.I2,
                 #
                 label.left = x$label.left,
                 label.right = x$label.right,
                 #
                 method.tau = "DL")
    #
    if (x$sm == "PFT") {
      m$n.harmonic.mean <- x$n.harmonic.mean
      m$n.harmonic.mean.pooled <- x$n.harmonic.mean.pooled
    }
    else if (x$sm == "IRFT") {
      m$t.harmonic.mean <- x$t.harmonic.mean
      m$t.harmonic.mean.pooled <- x$t.harmonic.mean.pooled
    }
    #
    m$pval <- pval
    m$statistic <- statistic
    m$I2 <- I2
    m$tau2 <- tau2
    m$tau <- tau
    #
    if (print.tau2.ci)
      m$tau2.ci <- tau2.ci
    #
    if (print.tau.ci)
      m$tau.ci <- tau.ci
    #
    if (print.I2.ci)
      m$I2.ci <- I2.ci
    #
    if (print.cidprop)
      m$prop.cid.below.null <- x$prop.cid.below.null
    #
    if (print.cidprop)
      m$prop.cid.above.null <- x$prop.cid.above.null
    #
    type.study <- type
  }
  #
  pval.pooled <- formatPT(x$pval.pooled, digits = digits.pval, lab.NA = lab.NA)
  statistic.pooled <- formatN(x$statistic.pooled, digits = digits.stat,
                              lab.NA = lab.NA, big.mark = big.mark)
  #
  m$TE.common <- x$TE.pooled
  m$lower.common <- x$lower.pooled
  m$upper.common <- x$upper.pooled
  m$statistic.common <- statistic.pooled
  m$pval.common <- pval.pooled
  #
  m$TE.random <- x$TE.pooled
  m$lower.random <- x$lower.pooled
  m$upper.random <- x$upper.pooled
  m$statistic.random <- statistic.pooled
  m$pval.random <- pval.pooled
  #
  m$df.random <- x$df.random.pooled
  #
  if (any(prediction)) {
    m$lower.predict <- x$lower.predict.pooled
    m$upper.predict <- x$upper.predict.pooled
    m$method.predict <- x$method.predict
  }
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
  m$k <- x$k.pooled
  m$k.study <- x$k.study.pooled
  m$k.all <- x$k.all.pooled
  m$k <- x$k.TE.pooled
  #
  m$.text.details.methods <- text.details
  
    
  data.p <-
    data.frame(
      tau2 = formatPT(x$tau2.pooled, digits = digits.tau2, lab.NA = lab.NA),
      tau = formatPT(x$tau.pooled, digits = digits.tau, lab.NA = lab.NA),
      I2 = percent_label(x$I2.pooled, "", digits.I2, formatPT, lab.NA = lab.NA)
      )
  #
  if (print.tau2.ci) {
    data.p$tau2.ci <-
      formatCI(x$lower.tau2.pooled, x$upper.tau2.pooled,
               digits = digits.tau2, lab.NA = lab.NA, na.lab = TRUE)
  }
  #
  if (print.tau.ci) {
    data.p$tau.ci <-
      formatCI(x$lower.tau.pooled, x$upper.tau.pooled,
               digits = digits.tau, lab.NA = lab.NA, na.lab = TRUE)
  }
  #
  if (print.I2.ci) {
    data.p$I2.ci <-
      formatCI(x$lower.I2.pooled, x$upper.I2.pooled,
               digits = digits.I2, lab.NA = lab.NA, na.lab = TRUE)
  }
  #
  if (print.cidprop.below.null) {
    data.p$prop.cid.below.null <-
      percent_label(x$prop.cid.below.null.pooled, "", digits.percent, formatPT,
                    lab.NA = lab.NA)
  }
  #
  if (print.cidprop.above.null) {
    data.p$prop.cid.above.null <-
      percent_label(x$prop.cid.above.null.pooled, "", digits.percent, formatPT,
                    lab.NA = lab.NA)
  }
  
  
  class(m) <- c(class(m), class(x))
  m$classes <- x$classes
  #
  forest_args <-
    list(x = m,
         #
         layout = layout,
         #
         leftcols = leftcols, leftlabs = leftlabs,
         rightcols = rightcols, rightlabs = rightlabs,
         #
         test.overall = FALSE,
         test.overall.common = FALSE, test.overall.random = FALSE,
         overall.hetstat = FALSE,
         #
         type.study = type.study,
         weight.study = "same",
         lab.NA = lab.NA,
         #
         main = main,
         data.pooled = data.p,
         just.addcols = just.addcols,
         #
         backtransf = backtransf,
         #
         big.mark = big.mark,
         #
         digits = list(digits = digits,
                       pval = digits.pval, stat = digits.stat,
                       tau2 = digits.tau2, tau = digits.tau,
                       I2 = digits.I2),
         #
         col = list(study = col, square = col.bg, square.lines = col.border),
         #
         details = details,
         #
         # Output
         #
         filename = filename,
         path = path,
         device = device, device.args = device.args,
         units = units,
         height = height, width = width,
         dpi = dpi,
         dev.off = dev.off,
         new = new
    )
  #
  res <- do.call("forest.meta", c(forest_args, args))
  #
  invisible(res)
}
