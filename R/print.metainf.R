#' Print results of a leave-one-out meta-analysis
#' 
#' @description
#' Print results of a leave-one-out meta-analysis
#' 
#' @aliases print.metainf
#' 
#' @param x An object of class \code{\link{metainf}}.
#' @param prediction A logical indicating whether prediction
#'   intervals should be printed.
#' @param overall A logical indicating whether overall results should be
#'   printed.
#' @param backtransf A logical indicating whether printed results
#'   should be back transformed. If \code{backtransf=TRUE}, results
#'   for \code{sm="OR"} are printed as odds ratios rather than log
#'   odds ratios, for example.
#' @param header A logical indicating whether information on title of
#'   meta-analysis, comparison and outcome should be printed at the
#'   beginning of the printout.
#' @param lab.NA A character string to label missing values.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.stat Minimal number of significant digits for z- or
#'   t-value of test for overall effect, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance \eqn{\tau^2}, see \code{print.default}.
#' @param digits.tau Minimal number of significant digits for
#'   \eqn{\tau}, the square root of the between-study variance
#'   \eqn{\tau^2}.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   and Rb statistic, see \code{print.default}.
#' @param digits.cid Minimal number of significant digits for
#'   CID / decision thresholds, see \code{print.default}.
#' @param digits.percent Minimal number of significant digits for
#'   probabilities, printed as percentages, see \code{print.default}.
#' @param big.mark A character used as thousands separator.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   overall effect should be printed according to JAMA reporting
#'   standards.
#' @param print.stat A logical value indicating whether z- or t-value
#'   for test of treatment effect should be printed.
#' @param print.tau2 A logical specifying whether between-study
#'   variance \eqn{\tau^2} should be printed.
#' @param print.tau2.ci A logical value indicating whether to print
#'   the confidence interval of \eqn{\tau^2}.
#' @param print.tau A logical specifying whether \eqn{\tau}, the
#'   square root of the between-study variance \eqn{\tau^2}, should be
#'   printed.
#' @param print.tau.ci A logical value indicating whether to print the
#'   confidence interval of \eqn{\tau}.
#' @param print.I2 A logical specifying whether heterogeneity
#'   statistic I\eqn{^2} should be printed.
#' @param print.I2.ci A logical specifying whether confidence interval for
#'   heterogeneity statistic I\eqn{^2} should be printed.
#' @param print.prob A logical specifying whether to print probabilities
#'   of clinically important benefit or harm.
#' @param text.tau2 Text printed to identify between-study variance
#'   \eqn{\tau^2}.
#' @param text.tau Text printed to identify \eqn{\tau}, the square
#'   root of the between-study variance \eqn{\tau^2}.
#' @param text.I2 Text printed to identify heterogeneity statistic
#'   I\eqn{^2}.
#' @param details.methods A logical specifying whether details on
#'   statistical methods should be printed.
#' @param \dots Additional arguments (ignored).
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metainf}}, \code{\link{settings.meta}}
#' 
#' @examples
#' data(Fleiss1993bin)
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac,
#'   data = Fleiss1993bin, studlab = study, sm = "RR", method = "I")
#' m1
#' metainf(m1)
#' metainf(m1, pooled = "random")
#' metainf(m1, pooled = "random", prediction = TRUE)
#'
#' @method print metainf
#' @export

print.metainf <- function(x,
                          #
                          prediction = x$prediction,
                          overall = x$overall,
                          backtransf = x$backtransf,
                          header = TRUE,
                          #
                          lab.NA = ".",
                          #
                          digits = gs("digits"),
                          digits.stat = gs("digits.stat"),
                          digits.pval = gs("digits.pval"),
                          digits.tau2 = gs("digits.tau2"),
                          digits.tau = gs("digits.tau"),
                          digits.I2 = gs("digits.I2"),
                          digits.cid = gs("digits.cid"),
                          digits.percent = 1,
                          #
                          big.mark = gs("big.mark"),
                          scientific.pval = gs("scientific.pval"),
                          zero.pval = gs("zero.pval"),
                          JAMA.pval = gs("JAMA.pval"),
                          #
                          print.stat = FALSE,
                          print.tau2 = TRUE,
                          print.tau2.ci = FALSE,
                          print.tau = TRUE,
                          print.tau.ci = FALSE,
                          print.I2 = TRUE,
                          print.I2.ci = FALSE,
                          print.prob = TRUE,
                          #
                          text.tau2 = gs("text.tau2"),
                          text.tau = gs("text.tau"),
                          text.I2 = gs("text.I2"),
                          #
                          details.methods = gs("details"),
                          ...) {
  
  chkclass(x, "metainf")
  #
  print.metacum(x,
                #
                prediction = prediction,
                overall = overall,
                backtransf = backtransf,
                header = header,
                #
                digits = digits,
                digits.stat = digits.stat,
                digits.pval = digits.pval,
                digits.tau2 = digits.tau2,
                digits.tau = digits.tau,
                digits.I2 = digits.I2,
                digits.cid = digits.cid,
                digits.percent = digits.percent,
                #
                big.mark = big.mark,
                scientific.pval = scientific.pval,
                zero.pval = zero.pval,
                JAMA.pval = JAMA.pval,
                #
                print.stat = print.stat,
                print.tau2 = print.tau2,
                print.tau2.ci = print.tau2.ci,
                print.tau = print.tau,
                print.tau.ci = print.tau.ci,
                print.I2 = print.I2,
                print.I2.ci = print.I2.ci,
                print.prob = print.prob,
                #
                text.tau2 = text.tau2,
                text.tau = text.tau,
                text.I2 = text.I2,
                #
                details.methods = details.methods,
                ...)
  #
  invisible(NULL)
}
