#' Return subset of pairwise object
#' 
#' @description
#' The \code{subset} method returns a subset of a pairwise object.
#' 
#' @param x An object of class \code{pairwise}.
#' @param subset A logical expression indicating rows to keep; missing values
#'   are taken as false.
#' @param \dots Additional arguments.
#' @return A pairwise object is returned.
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' @seealso \code{\link{pairwise}}, \code{\link[metadat]{dat.franchini2012}}
#' 
#' @examples
#' # Transform data from arm-based format to contrast-based format
#' pw1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
#'   n = list(n1, n2, n3),
#'   mean = list(y1, y2, y3), sd = list(sd1, sd2, sd3),
#'   data = dat.franchini2012, studlab = Study)
#' # Subset without Lieberman studies
#' subset(pw1, !grepl("Lieberman", studlab))[, 1:5]
#'
#' @method subset pairwise
#' @export

subset.pairwise <- function(x, subset, ...){
  
  chkclass(x, "pairwise")
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
  
  # Check whether 'reference.group' is available in the subset
  #
  reference.group <- attr(x, "reference.group")
  #
  trts <- sort(unique(c(res$treat1, res$treat2)))
  if (!(reference.group %in% trts))
    reference.group <- trts[1]
  
  
  # Return subset
  #
  for (i in names(attribs))
    attr(res, i) <- attr(x, i)
  #
  class(res) <- unique(c("pairwise", class(res)))
  #
  res
}
