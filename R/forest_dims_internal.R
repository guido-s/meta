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
