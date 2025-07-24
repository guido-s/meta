#' Extract parts of pairwise object
#' 
#' @description
#' Auxiliary function to extract parts of \code{\link{pairwise}} object.
#' 
#' @param x An object of class \code{\link{pairwise}}.
#' @param \dots Additional arguments (passed on to [.data.frame).
#' 
#' @aliases [,pairwise
#' 
#' @name [.pairwise
#' 
#' @usage \method{[}{pairwise}(x, ...)
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{pairwise}}, \code{\link{subset.pairwise}},
#'   \code{\link[metadat]{dat.franchini2012}}
#' 
#' @examples
#' # Transform data from arm-based format to contrast-based format
#' pw1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
#'   n = list(n1, n2, n3),
#'   mean = list(y1, y2, y3), sd = list(sd1, sd2, sd3),
#'   data = dat.franchini2012, studlab = Study)
#' 
#' pw1[, 1:5]
#' pw1[!grepl("Lieberman", pw1$studlab), 1:5]
#'
#' @method [ pairwise
#' @export

`[.pairwise` <- function(x, ...) {
  
  chkclass(x, "pairwise")
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
  class(res) <- unique(c("pairwise", class(res)))
  #
  res
}
