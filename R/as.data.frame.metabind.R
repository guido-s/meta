#' Coerce to a data frame
#' 
#' @description
#' The \code{as.data.frame} method returns a data frame containing
#' information on meta-analysis results.
#' 
#' @param x An object of class \code{metabind}.
#' @param row.names \code{NULL} or a character vector giving the row
#'   names for the data frame.
#' @param optional Argument is ignored.
#' @param \dots Other arguments (ignored).
#' 
#' @return A data frame is returned by the function \code{as.data.frame}.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metabind}}
#' 
#' @examples
#' data(Fleiss1993cont)
#' #
#' ma1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "SMD", studlab = paste(study, year))
#' #
#' ma2 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "SMD", studlab = paste(study, year),
#'   method.smd = "glass")
#' 
#' mb12 <-
#'   metabind(ma1, ma2, name = c("Hedges' g", "Glass' delta"), common = FALSE)
#' 
#' as.data.frame(mb12)
#'
#' @method as.data.frame metabind
#' @export

as.data.frame.metabind <- function(x, row.names = NULL, optional = FALSE, ...) {
  
  chkclass(x, "metabind")
  x <- updateversion(x)
  
  full.results <- bind_rows(x$common, x$random, x$predict)
  #
  res <- left_join(full.results, x$overall, by = "id")
  #
  if (!is.null(row.names)) {
    if (length(row.names) != nrow(res))
      stop("Argument 'row.names' must be of length ", nrow(res))
    row.names(res) <- row.names
  }
  else
    row.names(res) <- seq_len(nrow(res))
  #
  attr(res, "version") <- packageDescription("meta")$Version
  
  res
}
