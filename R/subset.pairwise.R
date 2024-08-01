#' Return subset of pairwise object
#' 
#' @description
#' The \code{subset} method returns a subset of a pairwise object.
#' 
#' @param x An object of class \code{pairwise}.
#' @param subset A logical expression indicating elements or rows to
#'   keep: missing values are taken as false.
#' @param \dots Additional arguments.
#' @return A pairwise object is returned.
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' @seealso \code{\link{pairwise}}
#' 
#' @examples
#' \dontrun{
#' # Transform data from arm-based format to contrast-based format
#' data(Franchini2012, package = "netmeta")
#' p1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
#'   n = list(n1, n2, n3),
#'   mean = list(y1, y2, y3), sd = list(sd1, sd2, sd3),
#'   data = Franchini2012, studlab = Study)
#' p1[, 1:5]
#'
#' # Subset without Lieberman studies
#' subset(p1, !grepl("Lieberman", studlab))[, 1:5]
#' }
#'
#' @method subset pairwise
#' @export


subset.pairwise <- function(x, subset, ...){
  
  chkclass(x, "pairwise")
  
  
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
  attr(res, "sm") <- attr(x, "sm")
  attr(res, "method") <- attr(x, "method")
  attr(res, "incr") <- attr(x, "incr")
  attr(res, "allincr") <- attr(x, "allincr")
  attr(res, "addincr") <- attr(x, "addincr")
  attr(res, "addstudies") <- attr(x, "addstudies")
  attr(res, "pairwise") <- attr(x, "pairwise")
  attr(res, "reference.group") <- attr(x, "reference.group")
  attr(res, "keep.all.comparisons") <- attr(x, "keep.all.comparisons")
  attr(res, "type") <- attr(x, "type")
  attr(res, "varnames") <- attr(x, "varnames")
  attr(res, "version") <- attr(x, "version")
  ##
  class(res) <- c("pairwise", class(res))
  ##
  res
}
