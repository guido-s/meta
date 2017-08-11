radial.default <- function(x, y,
                           xlim = NULL, ylim = NULL,
                           xlab = "Inverse of standard error",
                           ylab = "Standardised treatment effect (z-score)",
                           comb.fixed = TRUE, axes = TRUE,
                           pch = 1, text = NULL, cex = 1, col = NULL,
                           level = NULL, ...) {
  
  
  ##
  ##
  ## (1) Check essential arguments
  ##
  ##
  k.All <- length(x)
  ##
  chknumeric(x)
  chknumeric(y)
  ##
  fun <- "radial"
  chklength(y, k.All, fun)
  
  
  ##
  ##
  ## (2) Do meta-analysis
  ##
  ##
  m <- metagen(x, y)
  
  
  ##
  ##
  ## (3) Produce radial plot
  ##
  ##
  radial(m, xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab,
         comb.fixed = comb.fixed, axes = axes,
         pch = pch, text = text,
         cex = cex, col = col, level = level, ...)
  
  
  invisible(NULL)
}
