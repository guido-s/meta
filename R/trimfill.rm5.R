#' Cochrane review: trim-and-fill method
#' 
#' @description
#' Conduct trim-and-fill analysis for all meta-analyses in a Cochrane
#' review.
#' 
#' @aliases trimfill.rm5 trimfill.cdir 
#' 
#' @param x An object of class \code{rm5}, \code{trimfill.rm5},
#'   \code{cdir} or \code{trimfill.cdir}.
#' @param comp.no Comparison number.
#' @param outcome.no Outcome number.
#' @param ... Additional arguments (passed on to \code{metacr}).
#' 
#' @details
#' This function can be used to conduct a trim-and-fill analysis for
#' all or selected meta-analyses in a Cochrane review of intervention
#' studies (Higgins et al., 2023).
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
#' @method trimfill cdir
#' @export


trimfill.cdir <- function(x, comp.no, outcome.no, ...) {
  
  chkclass(x, "cdir")
  ##
  if (missing(comp.no))
    comp.no <- unique(x$data$comp.no)
  ##
  res <- list()
  ##
  n <- 1
  ##
  for (i in comp.no) {
    if (missing(outcome.no))
      jj <- unique(x$data$outcome.no[x$data$comp.no == i])
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





#' @rdname trimfill.rm5
#' @method print trimfill.cdir
#' @export


print.trimfill.cdir <- function(x, ...) {
  
  chkclass(x, "trimfill.cdir")
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
