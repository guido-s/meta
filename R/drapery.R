#' Drapery plot
#' 
#' @description
#' Draw a drapery plot with p-value curves for individual studies and
#' meta-analysis estimates.
#' 
#' @aliases drapery
#' 
#' @param x An object of class \code{meta}.
#' @param type A character string indicating whether to plot p-values
#'   (\code{"pvalue"}) or critical values (\code{"cvalue"}), can be
#'   abbreviated.
#' @param layout A character string for the line layout of individual
#'   studies: \code{"equal"}, \code{"grayscale"}, or
#'   \code{"linewidth"} (see Details), can be abbreviated.
#' @param lty.study Line type for individual studies.
#' @param lwd.study Line width for individual studies.
#' @param col.study Colour of lines for individual studies.
#' @param comb.fixed A logical indicating whether to show result for
#'   the fixed effect model.
#' @param comb.random A logical indicating whether to show result for
#'   the random effects model.
#' @param lty.fixed Line type for fixed effect meta-analysis.
#' @param lwd.fixed Line width for fixed effect meta-analysis.
#' @param col.fixed Colour of lines for fixed effect meta-analysis.
#' @param lty.random Line type for random effects meta-analysis.
#' @param lwd.random Line width for random effects meta-analysis.
#' @param col.random Colour of lines for random effects meta-analysis.
#' @param prediction A logical indicating whether to show prediction
#'   region.
#' @param col.predict Colour of prediction region
#' @param alpha Horizonal lines are printed for the specified alpha
#'   values.
#' @param col.alpha Colour of horizonal lines for alpha values.
#' @param col.null.effect Colour of vertical line indicating null
#'   effect.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param pos.legend Position of legend (see \code{\link{legend}}).
#' @param bg Background colour of legend (see \code{\link{legend}}).
#' @param bty Type of the box around the legend; either \code{"o"} or
#'   \code{"n"} (see \code{\link{legend}}).
#' @param backtransf A logical indicating whether results should be
#'   back transformed on the x-axis. For example, if \code{backtransf
#'   = FALSE}, log odds ratios instead of odds ratios are shown on the
#'   x-axis.
#' @param xlab A label for the x-axis.
#' @param ylab A label for the y-axis.
#' @param xlim The x limits (min, max) of the plot.
#' @param ylim The y limits (min, max) of the plot.
#' @param lwd.max The maximum line width (only considered if argument
#'   \code{layout} is equal to \code{"linewidth"}).
#' @param lwd.study.weight A character string indicating whether to
#'   determine line width for individual studies using weights from
#'   fixed effect (\code{"fixed"}) or random effects model
#'   (\code{"random"}), can be abbreviated (only considered if
#'   argument \code{layout} is equal to \code{"linewidth"}).
#' @param \dots Graphical arguments as in \code{par} may also be
#'   passed as arguments.
#' 
#' @details
#' The concept of a p-value function also called confidence curve goes
#' back to Birnbaum (1961). A drapery plot is showing p-value
#' functions for individual studies as well as meta-analysis estimates
#' is drawn in the active graphics window. Furthermore, a prediction
#' region for a single future study is shown as a shaded area. In
#' contrast to a forest plot, a drapery plot does not provide
#' information for a single confidence level however any confidence
#' level.
#'
#' Instead of p-value functions, curves with critical values can be
#' shown using argument \code{type = "cvalue"}.
#' 
#' Argument \code{layout} determines how curves for individual studies
#' are presented:
#' \itemize{
#' \item equal lines (\code{layout = "equal"})
#' \item darker gray tones with increasing precision (\code{layout = "grayscale"})
#' \item thicker lines with increasing precision (\code{layout = "linewidth"})
#' }
#' 
#' @author Gerta Rücker \email{sc@@imbi.uni-freiburg.de}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{forest}}, \code{\link{radial}}
#' 
#' @references
#' 
#' Birnbaum A (1961):
#' Confidence Curves: An Omnibus Technique for Estimation and Testing
#' Statistical Hypotheses.
#' \emph{Journal of the American Statistical Association},
#' \bold{56}, 246--9
#' 
#' @keywords hplot
#' 
#' @examples
#' data("lungcancer")
#' m1 <- metainc(d.smokers, py.smokers,
#'               d.nonsmokers, py.nonsmokers,
#'               data = lungcancer, studlab = study)
#' 
#' # Drapery plot
#' #
#' drapery(m1, xlim = c(0.5, 50))
#' 
#' @rdname drapery
#' @export drapery


drapery <- function(x, type = "pvalue", layout = "equal",
                    ##
                    lty.study = 2, lwd.study = 1, col.study = "black",
                    ##
                    comb.fixed = x$comb.fixed, comb.random = x$comb.random,
                    lty.fixed = 1, lwd.fixed = max(3, lwd.study),
                    col.fixed = "blue",
                    lty.random = 1, lwd.random = max(3, lwd.study),
                    col.random = "red",
                    ##
                    prediction = comb.random, col.predict = "lightblue",
                    ##
                    alpha = c(0.01, 0.05, 0.1),
                    col.alpha = "darkgray",
                    col.null.effect = "black",
                    ##
                    legend = TRUE, pos.legend = "topleft",
                    bg = "white", bty = "o",
                    ##
                    backtransf = x$backtransf,
                    xlab,
                    ylab =
                      if (type == "pvalue") "P-value" else "Critical value",
                    xlim, ylim,
                    lwd.max = 2.5,
                    lwd.study.weight = if (comb.random) "random" else "fixed",
                    ...) {
  
  
  pvf <- function(x, TE, seTE)
    2 * pnorm(abs(TE - x) / seTE, lower.tail = FALSE)
  

  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chkclass(x, "meta")
  ##
  type <- setchar(type, c("pvalue", "cvalue"))
  layout <- setchar(layout, c("equal", "grayscale", "linewidth"))
  ##
  chknumeric(lty.study, min = 0, zero = TRUE, single = TRUE)
  chknumeric(lwd.study, min = 0, zero = TRUE, single = TRUE)
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(prediction)
  ##
  if (!all(is.na(alpha)))
    chklevel(alpha, single = FALSE)
  ##
  chklogical(legend)
  bty <- setchar(bty, c("o", "n"))
  chklogical(backtransf)
  ##
  if (missing(xlab)) {
    xlab <- ""
    xlab <- xlab(x$sm, backtransf = backtransf)
  }
  else
    chkchar(xlab)
  chkchar(ylab)
  ##
  chknumeric(lwd.max, min = 0, zero = TRUE, single = TRUE)
  lwd.study.weight <- setchar(lwd.study.weight, c("random", "fixed"))
  
  
  ##
  ##
  ## (2) Set variables for plotting
  ##
  ##
  x.backtransf <- x$sm %in% c("RR", "OR", "HR", "IRR", "ROM") & backtransf
  ##
  if (missing(xlim)) {
    mn <- min(x$lower)
    mx <- max(x$upper)
    xlim <- c(mn, mx)
  }
  else {
    if (x.backtransf)
      xlim <- log(xlim)
    ##
    mn <- min(xlim)
    mx <- max(xlim)
  }
  ##
  grid <- sort(c(seq(mn, mx, length.out = 1000), x$TE, x$TE.fixed, x$TE.random))
  ##
  if (x.backtransf) {
    mn <- exp(mn)
    x.grid <- exp(grid)
    xlim <- exp(xlim)
    null.effect <- exp(x$null.effect)
  }
  else {
    x.grid <- grid
    null.effect <- x$null.effect
  }
  ##
  if (missing(ylim))
    ylim <- if (type == "pvalue") c(0, 1) else c(-3, 0)
  ##
  o <- order(x$seTE)
  TE <- x$TE[o]
  seTE <- x$seTE[o]
  w.fixed <- x$w.fixed[o]
  w.random <- x$w.random[o]
  ##
  k <- x$k
  TE.fixed <- x$TE.fixed
  seTE.fixed <- x$seTE.fixed
  TE.random <- x$TE.random
  seTE.random <- x$seTE.random
  seTE.predict <- x$seTE.predict
  t.quantile <- (x$upper.predict - x$lower.predict) / seTE.predict / 2
  ##
  lty.study <- rep(lty.study, k)
  lwd.study <- rep(lwd.study, k)
  col.study <- rep(col.study, k)
  ##
  if (layout == "grayscale")
    col.study <- gray.colors(k)
  ##
  else if (layout == "linewidth") {
    if (lwd.study.weight == "fixed")
      lwd.study <- lwd.max * w.fixed / max(w.fixed)
    else
      lwd.study <- lwd.max * w.random / max(w.random)
  }
  ##  
  y.1 <- pvf(grid, TE[1], seTE[1])
  y.fixed <- pvf(grid, TE.fixed, seTE.fixed)
  y.random <- pvf(grid, TE.random, seTE.random)
  y.predict <- pvf(grid, TE.random, t.quantile / qnorm(0.975) * seTE.predict)
  y.alpha <- alpha
  ##
  if (type == "cvalue") {
    y.1 <- qnorm(y.1 / 2)
    y.fixed <- qnorm(y.fixed / 2)
    y.random <- qnorm(y.random / 2)
    y.predict <- qnorm(y.predict / 2)
    y.alpha <- qnorm(y.alpha / 2)
  }
  ##
  sel.fixed <- y.fixed >= min(ylim)
  sel.random <- y.random >= min(ylim)
  sel.predict <- y.predict >= min(ylim)
  
  
  ##
  ##
  ## (3) Generate drapery plot
  ##
  ##
  plot(x.grid, y.1,
       xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
       type = "n", las = 1, if (x.backtransf) log = "x" else log = "", ...)
  ##
  ## Add prediction region
  ##
  if (prediction) {
    polygon(x.grid[sel.predict], y.predict[sel.predict],
            col = col.predict, border = NA)
    ##
    if (all(y.predict > min(ylim))) {
      x.min <- min(x.grid)
      x.max <- max(x.grid)
      y.min <- min(ylim)
      y.max <- min(y.predict)
      polygon(c(x.min, x.min, x.max, x.max),
              c(y.min, y.max, y.max, y.min),
              col = col.predict, border = NA)
    }
  }
  ##
  ## Add p-value lines and null effect line
  ##
  for (alpha.i in y.alpha)
    abline(h = alpha.i, col = col.alpha)
  ##
  abline(v = null.effect, col = col.null.effect)
  ##
  for (i in k:1) {
    y.i <- pvf(grid, TE[i], seTE[i])
    if (type == "cvalue")
      y.i <- qnorm(y.i / 2)
    sel.i <- y.i >= min(ylim)
    lines(x.grid[sel.i], y.i[sel.i],
          lty = lty.study[i], lwd = lwd.study[i], col = col.study[i])
  }
  ##
  ## Add fixed effect lines
  ##
  if (comb.fixed)
    lines(x.grid[sel.fixed], y.fixed[sel.fixed],
          lty = lty.fixed, lwd = lwd.fixed, col = col.fixed)
  ##
  ## Add random effects lines
  ##
  if (comb.random)
    lines(x.grid[sel.random], y.random[sel.random],
          lty = lty.random, lwd = lwd.random, col = col.random)
  ##
  ## Add text for alpha levels
  ##
  for (i in seq(along = y.alpha))
    text(mn, y.alpha[i] + (ylim[2] - ylim[1]) / 200,
         paste("p =", alpha[i]),
         cex = 0.7, col = col.alpha, adj = c(0, 0))
  ##
  ## Add legend
  ##
  if (legend & any(c(comb.fixed, comb.random, prediction)))
    legend(pos.legend,
           col = c(if (comb.fixed) col.fixed,
                   if (comb.random) col.random,
                   if (prediction) col.predict),
           lwd = c(if (comb.fixed) lwd.fixed,
                   if (comb.random) lwd.random,
                   if (prediction) lwd.random),
           bty = bty, bg = bg,
           c(if (comb.fixed) "Fixed effect model",
             if (comb.random)"Random effects model",
             if (prediction) "Range of prediction"))
  
  
  invisible(NULL)
}
