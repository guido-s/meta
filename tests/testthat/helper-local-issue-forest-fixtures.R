issue_forest_draw <- function(x, ...) {
  observed <- NULL
  draw <- meta:::forest_meta_internal
  testthat::local_mocked_bindings(
    forest_meta_internal = function(...) {
      observed <<- list(...)
      draw(...)
    }, .package = "meta")
  grDevices::pdf(NULL, width = 24, height = 16)
  on.exit(grDevices::dev.off())
  forest(x, ...)
  observed$widths.mm <- grid::convertWidth(observed$x1, "mm", valueOnly = TRUE)
  observed$label.widths.mm <- vapply(observed$cols$col.studlab$labels,
    function(g) grid::convertWidth(grid::grobWidth(g), "mm", valueOnly = TRUE),
    numeric(1))
  observed
}
