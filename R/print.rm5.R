#' Cochrane review: summary of meta-analyses
#' 
#' @description
#' Calculate and print a summary of all meta-analyses in a
#' Cochrane review.
#' 
#' @param x An object of class \code{rm5}.
#' @param comp.no Comparison number.
#' @param outcome.no Outcome number.
#' @param ... Additional arguments (passed on to \code{metacr}).
#' 
#' @details
#' This function can be used to redo all or selected meta-analyses of
#' a Cochrane Review of interventions (Higgins et al., 2023).
#' 
#' Review Manager 5 (RevMan 5) was the software used for preparing and
#' maintaining Cochrane Reviews. In RevMan 5, subgroup analyses can be defined
#' and data from a Cochrane review can be imported to R using the function
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
#' Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch
#' VA (editors) (2023):
#' \emph{Cochrane Handbook for Systematic Reviews of Interventions
#'   Version 6.4 (updated August 2023)}.
#' Available from \url{https://www.cochrane.org/authors/handbooks-and-manuals/handbook}
#' 
#' @examples
#' # Locate export data file "Fleiss1993_CR.csv"
#' # in sub-directory of package "meta"
#' #
#' filename <- system.file("extdata", "Fleiss1993_CR.csv", package = "meta")
#' Fleiss1993_CR <- read.rm5(filename)
#' 
#' # Print results for all meta-analysis
#' #
#' Fleiss1993_CR
#' 
#' # Print results only for second outcome of first comparison
#' #
#' print(Fleiss1993_CR, comp.no = 1, outcome.no = 2)
#' 
#' @method print rm5
#' @export


print.rm5 <- function(x, comp.no, outcome.no, ...) {
  
  
  ##
  ##
  ## (1) Check for rm5 object
  ##
  ##
  chkclass(x, "rm5")
  
  
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
      res[[n]] <- metacr(x, i, j, ...)
      ##
      n <- n + 1
    }
  }
  ##
  n <- 1
  ##
  for (i in seq_len(length(res))) {
    if (n > 1)
      cat("\n*****\n\n")
    print(res[[n]])
    n <- n + 1
  }
  
  
  invisible(NULL)
}
