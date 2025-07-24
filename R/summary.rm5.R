#' Cochrane review: detailed summary of meta-analyses
#' 
#' @description
#' Calculate and print a detailed summary of all meta-analyses in a
#' Cochrane review.
#' 
#' @aliases summary.rm5 summary.cdir 
#' 
#' @param object An object of class \code{rm5} or \code{cdir}.
#' @param x An object of class \code{summary.rm5} or
#'   \code{summary.cdir}.
#' @param comp.no Comparison number.
#' @param outcome.no Outcome number.
#' @param ... Additional arguments (passed on to \code{metacr}).
#' 
#' @details
#' This function can be used to rerun all or selected meta-analyses of
#' a Cochrane Review of interventions (Higgins et al., 2023).
#' 
#' The R function \code{\link{metacr}} is called internally.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{summary.meta}}, \code{\link{metacr}},
#'   \code{\link{read.rm5}}, \code{\link{read.cdir}},
#'   \code{\link{metabias.rm5}}, \code{\link{metabias.cdir}}
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
#' # Print summary results for all meta-analysis
#' #
#' summary(Fleiss1993_CR)
#' 
#' # Print summary results only for second outcome of first comparison
#' #
#' summary(Fleiss1993_CR, comp.no = 1, outcome.no = 2)
#' 
#' @method summary rm5
#' @export


summary.rm5 <- function(object, comp.no, outcome.no, ...) {
  
  chkclass(object, "rm5")
  ##
  if (missing(comp.no))
    comp.no <- unique(object$comp.no)
  ##
  res <- list()
  ##
  n <- 1
  ##
  for (i in comp.no) {
    if (missing(outcome.no))
      jj <- unique(object$outcome.no[object$comp.no == i])
    else
      jj <- outcome.no
    ##
    for (j in jj) {
      res[[n]] <- summary(metacr(object, i, j, ...))
      ##
      n <- n + 1
    }
  }
  ##
  class(res) <- "summary.rm5"
  
  res
}





#' @rdname summary.rm5
#' @method summary cdir
#' @export


summary.cdir <- function(object, comp.no, outcome.no, ...) {
  
  chkclass(object, "cdir")
  ##
  if (missing(comp.no))
    comp.no <- unique(object$data$comp.no)
  ##
  res <- list()
  ##
  n <- 1
  ##
  for (i in comp.no) {
    if (missing(outcome.no))
      jj <- unique(object$data$outcome.no[object$data$comp.no == i])
    else
      jj <- outcome.no
    ##
    for (j in jj) {
      res[[n]] <- summary(metacr(object, i, j, ...))
      ##
      n <- n + 1
    }
  }
  ##
  class(res) <- "summary.cdir"
  
  res
}





#' @rdname summary.rm5
#' @method print summary.rm5
#' @export


print.summary.rm5 <- function(x, ...) {
  
  chkclass(x, "summary.rm5")
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
  
  invisible(NULL)
}





#' @rdname summary.rm5
#' @method print summary.cdir
#' @export


print.summary.cdir <- function(x, ...) {
  
  chkclass(x, "summary.cdir")
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
  
  invisible(NULL)
}
