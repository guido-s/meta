# Auxiliary function for forest_dims()
#
# Package: meta
# Author:  Nour Edin Darwish <nouredindarwish@gmail.com>
# License: GPL (>= 2)
#

forest_dims_internal <- function(x, units = "in") {
  
  chklength(units, 1, text = "Argument 'units' must be of length 1.")
  units <- setchar(units, c("in", "cm", "mm"), pre = "either ")
  
  grid_unit <- c("in" = "inches", cm = "cm", mm = "mm")[[units]]
  old_dev <- dev.cur()
  pdf(file = NULL)
  on.exit(
    {
      dev.off()
      if (old_dev > 1)
        dev.set(old_dev)
    },
    add = TRUE)
  
  # Suppress grid.newpage hooks so that external hooks (e.g., the hook
  # registered by R CMD check to annotate each plot page with help("topic")
  # labels) don't inject extra viewports into the captured gTree.
  #
  old_hooks <- getHook("grid.newpage")
  setHook("grid.newpage", NULL, "replace")
  on.exit(setHook("grid.newpage", old_hooks, "replace"), add = TRUE)
  #
  gtree <- grid.grabExpr(do.call(forest_meta_internal, x))
  
  # The main viewport's layout sits at the vpTree parent
  #
  layout <- gtree$childrenvp[[1]]$parent$layout
  
  # Widths: exact per-column units, sum directly
  #
  width <- convertWidth(sum(layout$widths), grid_unit, valueOnly = TRUE)
  
  # Heights: exact per-row units, sum directly
  #
  height <- convertHeight(sum(layout$heights), grid_unit, valueOnly = TRUE)
  
  # Padding to prevent clipping of axis labels and text below the x-axis.
  #
  n_lines <- function(x)
    if (length(x) == 0 || all(x == "")) 0 else length(x)
  #
  bottom_positions <-
    if (n_lines(x$grob.xlab) == 0)
      numeric(0)
    else
      x$xlab.ypos - seq_len(n_lines(x$grob.xlab)) + 1
  #
  if (isTRUE(x$print.label) && isTRUE(x$bottom.lr))
    bottom_positions <-
      c(bottom_positions,
        x$y.bottom.lr - seq_len(max(n_lines(x$grob.label.left),
                                     n_lines(x$grob.label.right))) + 1)
  #
  bottom_depth <-
    if (length(bottom_positions) == 0)
      0
    else if (min(bottom_positions) >= 0)
      0
    else
      abs(min(bottom_positions)) + 0.5
  #
  extra_width <- inches2units(0.3, units)
  extra_height.default <- inches2units(0.3, units)
  extra_height <-
    if (bottom_depth == 0)
      if (isTRUE(attr(x, "details"))) 0 else extra_height.default
    else
      max(extra_height.default,
          convertHeight(unit(bottom_depth, "lines"),
                        grid_unit, valueOnly = TRUE))
  extra_height <- extra_height +
    convertHeight(unit(x$top.pad.lines, "lines"), grid_unit, valueOnly = TRUE)
  #
  res <- list(width = width + extra_width, height = height + extra_height,
              units = units)
  #
  res
}
