#' Cochrane review: Test for funnel plot asymmetry
#' 
#' @description
#' Conduct a test for funnel plot asymmetry for all outcomes in a
#' Cochrane review
#' 
#' @param x An object of class \code{rm5}.
#' @param comp.no Comparison number.
#' @param outcome.no Outcome number.
#' @param method.bias A character string indicating which test for
#'   small-study effects is to be used for all outcomes. Either
#'   \code{"rank"}, \code{"linreg"}, or \code{"mm"}, can be
#'   abbreviated. See function \code{\link{metabias}}
#' @param method.bias.binary A character string indicating which test
#'   is to be used for binary outcomes. Either \code{"rank"},
#'   \code{"linreg"}, \code{"mm"}, \code{"count"}, \code{"score"}, or
#'   \code{"peters"}, can be abbreviated. See function
#'   \code{\link{metabias}}
#' @param method.bias.or A character string indicating which test is
#'   to be used for binary outcomes with odds ratio as summary
#'   measure. Either \code{"rank"}, \code{"linreg"}, \code{"mm"},
#'   \code{"count"}, \code{"score"}, or \code{"peters"}, can be
#'   abbreviated. See function \code{\link{metabias}}
#' @param k.min Minimum number of studies to perform test for
#'   small-study effects.
#' @param ... Additional arguments (ignored at the moment)
#' 
#' @details
#' This function can be used to conduct a test for funnel plot
#' asymmetry for all or selected meta-analyses in a Cochrane Review.
#' 
#' Review Manager 5 (RevMan 5) is the current software used for
#' preparing and maintaining Cochrane Reviews
#' (\url{http://community.cochrane.org/tools/review-production-tools/revman-5}).
#' In RevMan 5, subgroup analyses can be defined and data from a
#' Cochrane review can be imported to R using the function
#' \code{read.rm5}.
#' 
#' The R function \code{\link{metacr}} is called internally.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{metabias}}, \code{\link{metacr}},
#'   \code{\link{read.rm5}}, \code{\link{summary.rm5}}
#' 
#' @references
#' Higgins, J.P.T and S. Green (2011):
#' \emph{Cochrane Handbook for Systematic Reviews of Interventions
#'   Version 5.1.0 [Updated March 2011]}.
#' The Cochrane Library: http://www.cochrane-handbook.org
#' 
#' @keywords htest
#' 
#' @examples
#' # Locate export data file "Fleiss93_CR.csv" in sub-directory of
#' # package "meta"
#' #
#' filename <- system.file("extdata", "Fleiss93_CR.csv", package = "meta")
#' Fleiss93_CR <- read.rm5(filename)
#' 
#' # Print results for all tests of small-study effects
#' #
#' metabias(Fleiss93_CR, k.min = 5)
#' 
#' # Print result of test of small-study effects for second outcome in
#' # first comparison
#' #
#' metabias(Fleiss93_CR, comp.no = 1, outcome.no = 2, k.min = 5)
#' 
#' @export
#' @export metabias.rm5


metabias.rm5 <- function(x, comp.no, outcome.no,
                         method.bias = "linreg",
                         method.bias.binary = method.bias,
                         method.bias.or = "score",
                         k.min = 10, ...) {
  
  
  ##
  ##
  ## (1) Check for rm5 object
  ##
  ##
  chkclass(x, "rm5")
  
  
  if (missing(comp.no))
    comp.no <- unique(x$comp.no)
  
  
  n <- 1
  ##
  for (i in comp.no) {
    if (missing(outcome.no))
      jj <- unique(x$outcome.no[x$comp.no == i])
    else
      jj <- outcome.no
    for (j in jj) {
      ##
      m1 <- metacr(x, i, j)
      ##
      if (inherits(m1, "metabin")) {
        if (m1$sm == "OR")
          mb1 <- metabias(m1, k.min = k.min, method.bias = method.bias.or)
        else 
          mb1 <- metabias(m1, k.min = k.min, method.bias = method.bias.binary)
      }
      else
        mb1 <- metabias(m1, k.min = k.min, method.bias = method.bias)
      ##
      if (!is.null(mb1$estimate)) {
        if (n > 1)
          cat("\n*****\n\n")
        print(mb1)
        ##
        n <- n + 1
      }
    }
  }
  
  invisible(NULL)
}
