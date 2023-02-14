#' Cochrane review: trim-and-fill method
#' 
#' @description
#' 
#' Conduct trim-and-fill analysis for all meta-analyses in a Cochrane
#' review.
#' 
#' @param x An object of class \code{rm5} or \code{trimfill.rm5}.
#' @param comp.no Comparison number.
#' @param outcome.no Outcome number.
#' @param ... Additional arguments (passed on to \code{metacr}).
#' 
#' @details
#' This function can be used to conduct a trim-and-fill analysis for
#' all or selected meta-analyses in a Cochrane review.
#' 
#' Review Manager 5 (RevMan 5) was the software used for preparing and
#' maintaining Cochrane Reviews
#' (\url{https://training.cochrane.org/online-learning/core-software/revman}).
#' Data from a Cochrane review can be imported to R using the function
#' \code{read.rm5}.
#' 
#' The R function \code{\link{metacr}} is called internally.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{summary.meta}}, \code{\link{metacr}},
#'   \code{\link{read.rm5}}, \code{\link{metabias.rm5}}
#' 
#' @references
#' Higgins, J.P.T and S. Green (2011):
#' \emph{Cochrane Handbook for Systematic Reviews of Interventions
#'   Version 5.1.0 [Updated March 2011]}.
#' The Cochrane Library: http://www.cochrane-handbook.org
#' 
#' @examples
#' # Locate export data file "Fleiss1993_CR.csv"
#' # in sub-directory of package "meta"
#' #
#' filename <- system.file("extdata", "Fleiss1993_CR.csv", package = "meta")
#' Fleiss1993_CR <- read.rm5(filename)
#' 
#' # Conduct trim-and-fill analysis
#' #
#' trimfill(Fleiss1993_CR)
#' 
#' # Conduct trim-and-fill analysis only for second outcome of first
#' # comparison
#' #
#' trimfill(Fleiss1993_CR, comp.no = 1, outcome.no = 2)
#' 
#' @method trimfill rm5
#' @export


trimfill.rm5 <- function(x, comp.no, outcome.no, ...) {
  
  chkclass(x, "rm5")
  ##
  if (missing(comp.no))
    comp.no <- unique(x$comp.no)
  ##
  res <- list()
  ##
  n <- 1
  ##
  for (i in comp.no) {
    if (missing(outcome.no))
      jj <- unique(x$outcome.no[x$comp.no == i])
    else
      jj <- outcome.no
    for (j in jj) {
      ##
      m1 <- metacr(x, i, j, ...)
      ##
      if (!is.na(m1$TE.common) & is.null(m1$subgroup)) {
        tf1 <- res[[n]] <- trimfill(m1)
        ##
        n <- n + 1
      }
    }
  }
  ##
  class(res) <- "trimfill.rm5"
  
  res
}





#' @rdname trimfill.rm5
#' @method print trimfill.rm5
#' @export


print.trimfill.rm5 <- function(x, ...) {
  
  chkclass(x, "trimfill.rm5")
  ##
  n <- 1
  ##
  for (i in 1:length(x)) {
    if (n > 1)
      cat("\n*****\n\n")
    ##
    print(x[[i]])
    ##
    n <- n + 1
  }
  ##
  invisible(NULL)
}
