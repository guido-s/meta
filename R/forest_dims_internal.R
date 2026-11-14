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
  
  extra_width <- inches2units(0.3, units)
  extra_height.default <- inches2units(0.3, units)
  overhang <- x$top.depth.lines + x$bottom.depth.lines
  extra_height <-
    if (overhang == 0)
      if (isTRUE(attr(x, "details"))) 0 else extra_height.default
    else
      max(extra_height.default,
          convertHeight(unit(overhang, "lines"),
                        grid_unit, valueOnly = TRUE))
  #
  res <- list(width = width + extra_width, height = height + extra_height,
              units = units)
  #
  res
}
