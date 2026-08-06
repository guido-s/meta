#' Forest plot to display the result of a leave-one-out meta-analysis
#' 
#' @description
#' Draws a forest plot in the active graphics window (using grid
#' graphics system).
#' 
#' @aliases forest.metainf
#' 
#' @param x An object of class \code{\link{metainf}}.
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
#' @param lab.NA A character string to label missing values.
#' @param layout A character string specifying the layout of the
#'   forest plot (see \code{\link{forest.meta}}).
#' @param calcwidth.details A logical indicating whether the first line
#'   of meta-analysis details should be considered to calculate width
#'   of columns on the left side of the forest plot.
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
#' @param details A logical specifying whether details on statistical
#'   methods should be printed.
#' @param \dots Additional graphical arguments (passed on to
#'   \code{\link{forest.meta}}).
#' 
#' @details
#' A forest plot, also called confidence interval plot, is drawn in
#' the active graphics window. Internally, R function
#' \code{\link{forest.metacum}} is called to produce the forest plot.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{forest.metacum}}, \code{\link{forest.meta}},
#'   \code{\link{metainf}}, \code{\link{settings.meta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' data(Fleiss1993bin)
#' ma <- metabin(d.asp, n.asp, d.plac, n.plac,
#'   data = Fleiss1993bin, studlab = study, sm = "RR", method = "I")
#' ma
#' metainf(ma)
#' metainf(ma, pooled = "random")
#' 
#' forest(metainf(ma))
#' forest(metainf(ma, pooled = "random"))
#' forest(metainf(ma, pooled = "random", prediction = TRUE))
#'
#' @method forest metainf
#' @export

forest.metainf <- function(x,
                           #
                           prediction = x$prediction,
                           overall = x$overall,
                           just.addcols = "right",
                           smlab = "Leave-One-Out Meta-Analysis",
                           type = "square",
                           #
                           layout = gs("layout"),
                           lab.NA = ".",
                           #
                           calcwidth.details = TRUE,
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
                           col,
                           col.bg,
                           col.border,
                           col.bg.predict,
                           col.border.predict,
                           #
                           details = gs("forest.details"),
                           ...) {
  
  chkclass(x, "metainf")
  #
  x$prediction <- prediction
  x$overall <- overall
  x$backtransf <- backtransf
  #
  layout <- setchar(layout, c("meta", "BMJ", "RevMan5", "JAMA"))
  #
  # Colour schemes for layouts:
  # - colors[1] - vertical line for common effect or random effects model
  # - colors[2] - diamond for meta-analysis results
  # - colors[3] - outer lines of diamonds
  # - colors[4] - square or circle colour
  # - colors[5] - outer lines of squares or circles
  # - colors[6] - line color (axes, reference line, header lines)
  # - colors[7] - prediction interval
  # - colors[8] - outer lines of prediction intervals
  # - colors[9] - subgroups
  # - colors[10] - color within squares or circles
  #
  colors <- layout_colors(layout)
  #
  if (missing(col))
    col <- colors[1]
  #
  if (missing(col.bg))
    col.bg <- ifelse(type == "diamond", colors[2], colors[4])
  #
  if (missing(col.border))
    col.border <- ifelse(type == "diamond", colors[3], colors[5])
  #
  if (missing(col.bg.predict))
    col.bg.predict <- colors[7]
  if (missing(col.border.predict))
    col.border.predict = colors[8]
  #
  res <- forest.metacum(x,
                        #
                        just.addcols = just.addcols,
                        smlab = smlab,
                        type = type,
                        layout = layout,
                        lab.NA = lab.NA,
                        calcwidth.details = calcwidth.details,
                        #
                        big.mark = big.mark,
                        digits = digits,
                        digits.pval = digits.pval,
                        digits.tau2 = digits.tau2,
                        digits.tau = digits.tau,
                        digits.I2 = digits.I2,
                        digits.cid = digits.cid,
                        digits.percent = digits.percent,
                        #
                        col = col,
                        col.bg = col.bg,
                        col.border = col.border,
                        col.bg.predict = col.bg.predict,
                        col.border.predict = col.border.predict,
                        #
                        details = details,
                        ...)
  #
  return(invisible(res))
}
