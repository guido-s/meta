#' Return subset of longarm object
#' 
#' @description
#' The \code{subset} method returns a subset of a longarm object.
#' 
#' @param x An object of class \code{longarm}.
#' @param subset A logical expression indicating elements or rows to
#'   keep: missing values are taken as false.
#' @param \dots Additional arguments.
#' @return A \code{\link{longarm}} object is returned.
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' @seealso \code{\link{longarm}}
#' 
#' @examples
#' # Artificial example with three studies
#' m <- metabin(1:3, 100:102, 4:6, 200:202, studlab = LETTERS[1:3])
#' # Transform data to long arm-based format
#' l1 <- longarm(m)
#' l1
#'
#' # Subset without Study B
#' subset(l1, studlab != "B")
#'
#' @method subset longarm
#' @export


subset.longarm <- function(x, subset, ...){
  
  chkclass(x, "longarm")
  
  
  ## Catch 'subset'
  ##
  if (missing(subset))
    return(x)
  ##
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  subset <- catch("subset", mc, x, sfsp)
  if (is.null(subset))
    subset <- catch("subset", mc, x$data, sfsp)
  
  
  ## Select subset
  ##
  x.df <- x
  class(x.df) <- "data.frame"
  ##
  res <- subset(x.df, subset, ...)
  
  
  ## Return subset
  ## 
  attr(res, "type") <- attr(x, "type")
  #
  class(res) <- unique(c("longarm", class(res)))
  #
  res
}
