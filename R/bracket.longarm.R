#' Extract parts of longarm object
#' 
#' @description
#' Auxiliary function to extract parts of \code{\link{longarm}} object.
#' 
#' @param x An object of class \code{\link{longarm}}.
#' @param \dots Additional arguments (passed on to [.data.frame).
#' 
#' @aliases [,longarm
#' 
#' @name [.longarm
#' 
#' @usage \method{[}{longarm}(x, ...)
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{longarm}}, \code{\link{subset.longarm}},
#'   \code{\link[metadat]{dat.franchini2012}}
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
#' la1[la1$studlab == "Lieberman 1998", ]
#'
#' @method [ longarm
#' @export

`[.longarm` <- function(x, ...) {
  
  chkclass(x, "longarm")
  #
  attribs <- attributes(x)
  #
  attribs$names <- attribs$row.names <- NULL
  
  res <- x
  class(res) <- "data.frame"
  res <- res[...]
  
  # Check whether 'reference.group' is available in the subset
  #
  reference.group <- attr(x, "reference.group")
  #
  trts <- sort(unique(c(res$treat1, res$treat2)))
  if (!is.null(trts) && !(reference.group %in% trts))
    reference.group <- trts[1]
  
  # Return results
  # 
  for (i in names(attribs))
    attr(res, i) <- attr(x, i)
  #
  class(res) <- unique(c("longarm", class(res)))
  #
  res
}
