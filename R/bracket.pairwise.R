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
#' p1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
#'   n = list(n1, n2, n3),
#'   mean = list(y1, y2, y3), sd = list(sd1, sd2, sd3),
#'   data = dat.franchini2012, studlab = Study)
#' 
#' p1[, 1:5]
#' p1[!grepl("Lieberman", p1$studlab), 1:5]
#'
#' @method [ pairwise
#' @export

`[.pairwise` <- function(x, ...) {
  
  chkclass(x, "pairwise")
  
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
  attr(res, "sm") <- attr(x, "sm")
  attr(res, "method") <- attr(x, "method")
  attr(res, "incr") <- attr(x, "incr")
  attr(res, "allincr") <- attr(x, "allincr")
  attr(res, "addincr") <- attr(x, "addincr")
  attr(res, "addstudies") <- attr(x, "addstudies")
  attr(res, "pairwise") <- attr(x, "pairwise")
  attr(res, "reference.group") <- reference.group
  attr(res, "keep.all.comparisons") <- attr(x, "keep.all.comparisons")
  attr(res, "type") <- attr(x, "type")
  attr(res, "varnames") <- attr(x, "varnames")
  attr(res, "version") <- attr(x, "version")
  #
  class(res) <- c("pairwise", class(res))
  #
  res
}
