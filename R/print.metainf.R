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
#' @param print A named list controlling displayed information; see
#'   \code{\link{meta-formatting}}. By default, \code{tau2}, \code{tau},
#'   \code{I2}, and \code{prob} are enabled; \code{stat} and heterogeneity
#'   confidence intervals are disabled.
#' @param digits A named list controlling rounding; see
#'   \code{\link{meta-formatting}}. A scalar value is also accepted.
#'   Defaults are the same as in \code{\link{print.metacum}}.
#' @param big.mark A character used as thousands separator.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   overall effect should be printed according to JAMA reporting
#'   standards.
#' @param text.tau2 Text printed to identify between-study variance
#'   \eqn{\tau^2}.
#' @param text.tau Text printed to identify \eqn{\tau}, the square
#'   root of the between-study variance \eqn{\tau^2}.
#' @param text.I2 Text printed to identify heterogeneity statistic
#'   I\eqn{^2}.
#' @param details.methods A logical specifying whether details on
#'   statistical methods should be printed.
#' @param \dots Additional formatting arguments, including individual
#'   \code{print.*} and \code{digits.*} arguments.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metainf}}, \code{\link{settings.meta}}
#' 
#' @examples
#' data(Fleiss1993bin)
#' ma <- metabin(d.asp, n.asp, d.plac, n.plac,
#'   data = Fleiss1993bin, studlab = study, sm = "RR", method = "I")
#' ma
#' metainf(ma)
#' metainf(ma, pooled = "random")
#' metainf(ma, pooled = "random", prediction = TRUE)
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
                          print = NULL,
                          digits = NULL,
                          #
                          big.mark = gs("big.mark"),
                          scientific.pval = gs("scientific.pval"),
                          zero.pval = gs("zero.pval"),
                          JAMA.pval = gs("JAMA.pval"),
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
                print = print,
                digits = digits,
                lab.NA = lab.NA,
                #
                big.mark = big.mark,
                scientific.pval = scientific.pval,
                zero.pval = zero.pval,
                JAMA.pval = JAMA.pval,
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
