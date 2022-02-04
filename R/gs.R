#' Get default for a meta-analysis setting.
#' 
#' @description
#' Get default for a meta-analysis setting in R package \bold{meta}.
#' 
#' @param x A character string holding a settings name.
#' 
#' @details
#' This function can be used to get the default for a meta-analysis
#' setting defined using \code{\link{settings.meta}}.
#' 
#' This function is primarily used to define default settings in
#' meta-analysis functions, e.g., \code{\link{metabin}} or
#' \code{\link{metacont}}. A list of all arguments with current
#' settings is printed using the command
#' \code{settings.meta("print")}.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{settings.meta}}
#' 
#' @examples
#' # Get default setting for Hartung-Knapp method
#' #
#' gs("hakn")
#' 
#' # Get default setting for summary measure in metabin()
#' #
#' gs("smbin")
#' 
#' @export gs


gs <- function(x) {
  
  if (missing(x))
    return(NULL)
  
  chkchar(x)
  ##
  x.nam <- deparse(substitute(x))
  x.nam <- substring(x.nam, 2, nchar(x.nam) - 1)
  ##
  if (!(x %in% .settings$argslist.internal))
    x <- setchar(x, .settings$argslist, "unmatched", name = x.nam,
                 stop.at.error = FALSE)
  
  if (!is.null(x))
    x <- settings.meta(quietly = TRUE)[[x]]
  
  x
}
