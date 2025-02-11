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
#' # Transform data from wide arm-based format to contrast-based format
#' pw1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
#'   n = list(n1, n2, n3),
#'   mean = list(y1, y2, y3), sd = list(sd1, sd2, sd3),
#'   data = dat.franchini2012, studlab = Study)
#' 
#' # Transform data from contrast-based to long arm-based format
#' # and only keep the main variables
#' la1 <- longarm(pw1, append = FALSE)
#' head(la1)
#'
#' # Subset without Lieberman studies
#' subset(la1, grepl("Lieberman", studlab))
#'
#' @method subset longarm
#' @export


subset.longarm <- function(x, subset, ...){
  
  chkclass(x, "longarm")
  #
  attribs <- attributes(x)
  #
  attribs$names <- attribs$row.names <- NULL
  
  
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
  for (i in names(attribs))
    attr(res, i) <- attr(x, i)
  #
  class(res) <- unique(c("longarm", class(res)))
  #
  res
}
