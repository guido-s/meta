#' Coerce to a data frame
#' 
#' @description
#' The \code{as.data.frame} method returns a data frame containing
#' information on individual studies, e.g., estimated treatment effect
#' and its standard error.
#' 
#' @param x An object of class \code{meta}.
#' @param row.names \code{NULL} or a character vector giving the row
#'   names for the data frame.
#' @param optional logical. If \code{TRUE}, setting row names and
#'   converting column names (to syntactic names) is optional.
#' @param \dots other arguments
#' 
#' @return A data frame is returned by the function \code{as.data.frame}.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{metacont}},
#'   \code{\link{metagen}}, \code{\link{forest.meta}}
#' 
#' @examples
#' data(Fleiss1993cont)
#' #
#' # Generate additional variable with grouping information
#' #
#' Fleiss1993cont$group <- c(1, 2, 1, 1, 2)
#' #
#' # Do meta-analysis without grouping information
#' #
#' m1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "SMD", studlab = paste(study, year))
#' #
#' # Update meta-analysis object and do subgroup analyses
#' #
#' update(m1, subgroup = group)
#' 
#' # Same result using metacont function directly
#' #
#' m2 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "SMD", studlab = paste(study, year),
#'   subgroup = group)
#' m2
#' 
#' # Compare printout of the following two commands
#' #
#' as.data.frame(m1)
#' m1$data
#'
#' @method as.data.frame meta
#' @export


as.data.frame.meta <- function(x, row.names = NULL, optional = FALSE, ...) {
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(x, "meta")
  
  
  ## Remove element 'call' from object of class meta to get rid
  ## of an error message in meta-analyses with six studies:
  ## 'Error: evaluation nested too deeply: infinite recursion ...'
  ##
  ## NB: Element 'call' which is of length six contains information
  ##     on the function call.
  ##
  x$call <- NULL
  x$call.object <- NULL
  
  ## Remove data set from output
  ##
  x$data <- NULL
  
  ## Remove risk of bias table from output
  ##
  x$rob <- NULL
  
  ## Remove list information from metabind()
  ##
  x$list <- NULL
  
  ## Remove debug information
  ##
  x$debug <- x$tau2.calc <- x$rma.three.level <- NULL
  
  if (!is.null(x$approx.TE) && all(x$approx.TE == ""))
    x$approx.TE <- NULL
  ##
  if (!is.null(x$approx.seTE) && all(x$approx.seTE == ""))
    x$approx.seTE <- NULL
  
  sel <- as.vector(lapply(x, length) == length(x$TE))
  
  res <- as.data.frame(x[names(x)[sel]], ...)
  ##
  if (!is.null(row.names)) {
    if (length(row.names) != nrow(res))
      stop("Argument 'row.names' must be of length ", nrow(res))
    row.names(res) <- row.names
  }
  else
    row.names(res) <- seq_len(nrow(res))
  
  attr(res, "version") <- packageDescription("meta")$Version
  
  res
}
