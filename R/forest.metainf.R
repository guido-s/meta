#' Forest plot to display the result of a leave-one-out meta-analysis
#' 
#' @description
#' Draws a forest plot in the active graphics window (using grid
#' graphics system).
#' 
#' @aliases forest.metainf
#' 
#' @param x An object of class \code{\link{metainf}}.
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
#' @seealso \code{\link{forest.meta}}, \code{\link{metainf}},
#'   \code{\link{settings.meta}}
#' 
#' @keywords hplot
#' 
#' @examples
#' data(Fleiss1993bin)
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac,
#'   data = Fleiss1993bin, studlab = study, sm = "RR", method = "I")
#' m1
#' metainf(m1)
#' metainf(m1, pooled = "random")
#' 
#' forest(metainf(m1))
#' forest(metainf(m1, pooled = "random"))
#' forest(metainf(m1, pooled = "random", prediction = TRUE))
#'
#' @method forest metainf
#' @export

forest.metainf <- function(x,
                           #
                           leftcols = NULL, leftlabs = NULL,
                           rightcols = NULL, rightlabs = NULL,
                           #
                           prediction = x$prediction,
                           overall = x$overall,
                           just.addcols = "right",
                           smlab = "Leave-One-Out Meta-Analysis",
                           type = "square",
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
  
  chkclass(x, "metainf")
  #
  res <- forest.metacum(x,
                        leftcols = leftcols, leftlabs = leftlabs,
                        rightcols = rightcols,
                        rightlabs = rightlabs,
                        prediction = prediction,
                        overall = overall,
                        just.addcols = just.addcols,
                        smlab = smlab,
                        type = type,
                        lab.NA = lab.NA,
                        #
                        backtransf = backtransf,
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
                        addrows.below.overall = addrows.below.overall,
                        details = details,
                        ...)
  #
  return(invisible(res))
}
