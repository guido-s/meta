#' Auxiliary functions to calculate the dimensions of a forest plot
#'
#' @description
#' Extracts the exact width and height of a `meta` forest plot from its internal
#' grid layout. This is primarily used to determine the optimal dimensions for
#' saving the forest plot to a file.
#' 
#' @details
#' Rather than using guesswork or manual row counting, `forest_dims()` captures
#' the forest plot as a true graphics object and extracts the precise width and
#' height directly from its underlying structure.
#'
#' Because it mathematically measures the actual rendered components, it is
#' highly robust. It works seamlessly with any forest plot configuration.
#'
#' @param x An object of class \code{meta} (e.g., from \code{\link{metacont}},
#'   \code{\link{metabin}}, or \code{\link{metagen}}).
#' @param \dots Additional arguments passed on to the underlying forest plot
#'   method (e.g., \code{\link{forest.meta}}, \code{\link{forest.metabind}},
#'    \code{\link{forest.metacum}}, or \code{\link{forest.metainf}}).
#' @param units Units of the returned `width` and `height`. One of
#'   \code{"in"} (for inches, default), \code{"cm"}, or \code{"mm"}, can be
#'   abbreviated. In \code{inches2units} and \code{units2inches}, \code{"px"}
#'   (for pixels) is also admissible.
#' @param dpi Plot resolution.
#'
#' @return A named list with elements:
#' \item{width}{Forest plot width.}
#' \item{height}{Forest plot height.}
#' \item{units}{Units of `width` and `height`.}
#' 
#' @author Nour Edin Darwish \email{nouredindarwish@gmail.com}
#'
#' @examples
#' # Create a simple meta-analysis object
#' ma <- metagen(TE = c(0.5, 0.8, 0.3), seTE = c(0.2, 0.3, 0.15),
#'   studlab = paste("Study", LETTERS[1:3]))
#'
#' # Get dimensions in inches (default)
#' forest_dims(ma)
#'
#' # Get dimensions in centimetres
#' forest_dims(ma, units = "cm")
#'
#' # Forest plot with details on meta-analysis methods
#' forest_dims(ma, details = TRUE, units = "cm")
#'
#' @export

forest_dims <- function(x, ..., units = "in") {
  
  chklength(units, 1, text = "Argument 'units' must be of length 1.")
  units <- setchar(units, c("in", "cm", "mm"), pre = "either ")
  
  # Capture ... as unevaluated expressions so that NSE arguments (e.g.,
  # sortvar = TE) are forwarded as-is to forest(), which resolves
  # them via match.call() + catch().
  #
  dots_exprs <- enexprs(...)
  user_env <- parent.frame()
  filename <- tempfile(fileext = ".pdf")
  on.exit(unlink(filename), add = TRUE)
  #
  call_expr <-
    expr(forest(!!x, !!!dots_exprs, filename = !!filename, units = !!units))
  res <- eval(call_expr, envir = user_env)
  res <- res[c("width", "height", "units")]
  #
  res
}


#' @rdname forest_dims
#' @export

units2inches <- function(x, units = "in", dpi = 300) {
  chklength(units, 1, text = "Argument 'units' must be of length 1.")
  units <- setchar(units, c("in", "cm", "mm", "px"), pre = "either ")
  #
  x / c("in" = 1, cm = 2.54, mm = 2.54 * 10, px = dpi)[[units]]
}


#' @rdname forest_dims
#' @export

inches2units <- function(x, units = "in", dpi = 300) {
  chklength(units, 1, text = "Argument 'units' must be of length 1.")
  units <- setchar(units, c("in", "cm", "mm", "px"), pre = "either ")
  #
  x * c("in" = 1, cm = 2.54, mm = 2.54 * 10, px = dpi)[[units]]
}
