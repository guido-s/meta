#' Cochrane review: Test for funnel plot asymmetry
#' 
#' @description
#' Conduct a test for funnel plot asymmetry for all outcomes in a
#' Cochrane review of intervention studies
#' 
#' @aliases metabias.rm5 metabias.cdir 
#' 
#' @param x An object of class \code{rm5} or \code{cdir}.
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
#' asymmetry for all or selected meta-analyses in a Cochrane review of
#' intervention studies (Higgins et al, 2023).
#' 
#' The R function \code{\link{metacr}} is called internally.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metabias}}, \code{\link{metacr}},
#'   \code{\link{read.rm5}}, \code{\link{read.cdir}},
#'   \code{\link{summary.rm5}}, \code{\link{summary.cdir}}
#' 
#' @references
#' Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch
#' VA (editors) (2023):
#' \emph{Cochrane Handbook for Systematic Reviews of Interventions
#'   Version 6.4 (updated August 2023)}.
#' Available from \url{https://www.cochrane.org/authors/handbooks-and-manuals/handbook}
#' 
#' @keywords htest
#' 
#' @examples
#' # Locate export data file "Fleiss1993_CR.csv" in sub-directory of
#' # package "meta"
#' #
#' filename <- system.file("extdata", "Fleiss1993_CR.csv", package = "meta")
#' Fleiss1993_CR <- read.rm5(filename)
#' 
#' # Print results for all tests of small-study effects
#' #
#' metabias(Fleiss1993_CR, k.min = 5)
#' 
#' # Print result of test of small-study effects for second outcome in
#' # first comparison
#' #
#' metabias(Fleiss1993_CR, comp.no = 1, outcome.no = 2, k.min = 5)
#' 
#' @method metabias rm5
#' @export


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





#' @rdname metabias.rm5
#' @method metabias cdir
#' @export


metabias.cdir <- function(x, comp.no, outcome.no,
                          method.bias = "linreg",
                          method.bias.binary = method.bias,
                          method.bias.or = "score",
                          k.min = 10, ...) {
  
  
  ##
  ##
  ## (1) Check for cdir object
  ##
  ##
  chkclass(x, "cdir")
  
  
  if (missing(comp.no))
    comp.no <- unique(x$data$comp.no)
  
  
  n <- 1
  ##
  for (i in comp.no) {
    if (missing(outcome.no))
      jj <- unique(x$data$outcome.no[x$data$comp.no == i])
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
