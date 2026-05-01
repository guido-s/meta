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
#'   \code{\link{metabin}}, or \code{\link{metagen}}.
#' @param \dots Additional arguments passed on to the underlying forest plot
#'   method (e.g., \code{\link{forest.meta}}, \code{\link{forest.metabind}},
#'    \code{\link{forest.metacum}}, or \code{\link{forest.metainf}}.
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
  #
  grid_unit <- c("in" = "inches", cm = "cm", mm = "mm")[[units]]
  
  # Capture ... as unevaluated expressions so that NSE arguments (e.g.,
  # sortvar = TE) are forwarded as-is to forest(), which resolves
  # them via match.call() + catch().
  #
  dots_exprs <- enexprs(...)
  user_env <- caller_env()
  call_expr <- expr(forest(!!x, !!!dots_exprs))
  
  # convertWidth()/convertHeight() below require an open graphics device to
  # resolve grid unit conversions, so we open a null PDF device (writes
  # nowhere) for that purpose.
  #
  # Crucially, this device must be opened BEFORE grid.grabExpr(), not after.
  # This placement is entirely unrelated to convertWidth() — it's to fix a
  # separate bug where grid.grabExpr() creates an unwanted "Rplots.pdf" file
  # when running tests via the RStudio "Test" button or R CMD check.
  #
  # Here is how that bug happens:
  # grid.grabExpr() saves the current device at start (cd <- dev.cur()) and
  # restores it on exit via dev.set(cd). When no device is open before the
  # call, cd = 1 (the "null device"). grid.grabExpr() opens its own offscreen
  # device, evaluates, grabs the display list, then on cleanup: closes the
  # offscreen device and calls dev.set(1). At this exact moment there are ZERO
  # open devices. Calling dev.set(1) with no devices open triggers R to
  # automatically open a new device via getOption("device").
  #
  # Which device getOption("device") auto-opens depends on the environment.
  # In an interactive RStudio console the default device is "RStudioGD" (the
  # plot pane) — harmless, no file is written. In a non-interactive subprocess
  # (RStudio "Test" button) the default device is pdf(), which without a file
  # argument writes to "Rplots.pdf" — an unwanted side-effect.
  #
  # Opening pdf(file = NULL) here first ensures grid.grabExpr() sees
  # dev.cur() >= 2 (our null PDF), not 1. On exit it calls dev.set(2), which
  # harmlessly restores our null PDF instead of triggering the auto-open.
  #
  # Note: our own on.exit below does NOT suffer from this problem because the
  # `if (old_dev > 1)` guard skips dev.set() when old_dev is 1 (the null
  # device), so we never call dev.set(1) ourselves.
  #
  old_dev <- dev.cur()
  pdf(file = NULL)
  on.exit({
    dev.off()
    if (old_dev > 1) dev.set(old_dev)
    },
    add = TRUE)
  
  # Suppress grid.newpage hooks so that external hooks (e.g., the hook
  # registered by R CMD check to annotate each plot page with help("topic")
  # labels) don't inject extra viewports into the captured gTree.
  #
  old_hooks <- getHook("grid.newpage")
  setHook("grid.newpage", NULL, "replace")
  on.exit(setHook("grid.newpage", old_hooks, "replace"), add = TRUE)
  
  gtree <- grid.grabExpr(eval(call_expr, envir = user_env))
  
  # The main viewport's layout sits at the vpTree parent
  #
  layout <- gtree$childrenvp[[1]]$parent$layout
  
  # Widths: exact per-column units, sum directly
  #
  width <- convertWidth(sum(layout$widths), grid_unit, valueOnly = TRUE)
  
  # Heights: single unit recycled across nrow — expand then sum
  #
  height <-
    convertHeight(sum(rep(layout$heights, layout$nrow)),
                  grid_unit, valueOnly = TRUE)
  
  # Hardcoded padding to prevent clipping of axis labels
  #
  extra_width <- inches2units(0.3, units)
  extra_height <- inches2units(0.8, units)
  #
  res <- list(width = width + extra_width, height = height + extra_height,
              units = units)
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
