#' Summary method for objects of class metabind
#'
#' @description
#' Summary method for objects of class \code{metabind}. Summaries are printed
#' using a wrapper for \code{\link{print.metabind}}, providing a standard
#' \code{summary()} interface for combined meta-analysis results.
#'
#' @param object An object of class \code{metabind}.
#' @param x An object of class \code{"summary.metabind"}.
#' @param common A logical indicating whether common effect estimates
#'   should be printed.
#' @param random A logical indicating whether random effects estimates
#'   should be printed.
#' @param prediction A logical indicating whether prediction intervals
#'   should be printed.
#' @param overall A logical indicating whether overall summaries
#'   should be printed. This argument is useful in a meta-analysis
#'   with subgroups if summaries should only be printed on group
#'   level.
#' @param print A named list specifying which information to print; see
#'   \code{\link{meta-formatting}}.
#' @param digits A named list specifying digits used for displayed results; a
#'   scalar value is also accepted. See \code{\link{meta-formatting}}.
#' @param text.tau2 Text printed to identify between-study variance
#'   \eqn{\tau^2}.
#' @param text.tau Text printed to identify \eqn{\tau}, the square
#'   root of the between-study variance \eqn{\tau^2}.
#' @param text.I2 Text printed to identify heterogeneity statistic
#'   I\eqn{^2}.
#' @param lab.NA A character string to label missing values.
#' @param big.mark A character used as thousands separator.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param zero.pval A logical specifying whether p-values should be
#'   printed with a leading zero.
#' @param JAMA.pval A logical specifying whether p-values for test of
#'   overall effect should be printed according to JAMA reporting
#'   standards.
#' @param details.methods A logical specifying whether details on
#'   statistical methods should be printed.
#' @param \dots Additional arguments, ignored by \code{summary.metabind}
#'   and passed to \code{\link{print.metabind}} by \code{print.summary.metabind}.
#'
#' @return
#' An object of class \code{"summary.metabind"} with a corresponding print
#' method.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' @seealso \code{\link{metabind}}, \code{\link{print.metabind}}
#'
#' @examples
#' data(Fleiss1993bin)
#' #
#' ma1 <- metabin(d.asp, n.asp, d.plac, n.plac,
#'   data = Fleiss1993bin, studlab = study, sm = "RR", method = "I")
#' ma2 <- update(ma1, method = "MH", random = FALSE)
#' #
#' mb <- metabind(ma1, ma2)
#' #
#' # Identical printouts
#' summary(mb)
#' mb
#'
#' @method summary metabind 
#' @export

summary.metabind <- function(object, ...) {
  chkclass(object, "metabind")
  #
  res <- object
  #
  class(res) <- c("summary.metabind", class(res))
  #
  res
}


#' @rdname summary.metabind
#' @method print summary.metabind 
#' @export

print.summary.metabind <- function(x,
                                   #
                                   common = !is.null(x$common),
                                   random = !is.null(x$random),
                                   prediction = !is.null(x$predict),
                                   #
                                   overall = common || random,
                                   #
                                   print = NULL,
                                   digits = NULL,
                                   #
                                   lab.NA = ".",
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
  
  chkclass(x, "summary.metabind")
  #
  class(x) <- "metabind"
  #
  print(x,
        #
        common = common,
        random = random,
        prediction = prediction,
        #
        overall = overall,
        #
        print = print,
        digits = digits,
        #
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
