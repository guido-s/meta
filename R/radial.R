#' Radial plot
#' 
#' @description
#' Draw a radial plot (also called Galbraith plot) which can be used
#' to assess bias in meta-analysis.
#' 
#' @aliases radial radial.default radial.meta
#' 
#' @param x An object of class \code{meta}, or estimated treatment
#'   effect in individual studies.
#' @param y Standard error of estimated treatment effect.
#' @param xlim The x limits (min, max) of the plot.
#' @param ylim The y limits (min, max) of the plot.
#' @param xlab A label for the x-axis.
#' @param ylab A label for the y-axis.
#' @param common A logical indicating whether the pooled common effect
#'   estimate should be plotted.
#' @param axes A logical indicating whether axes should be drawn on
#'   the plot.
#' @param pch The plotting symbol used for individual studies.
#' @param text A character vector specifying the text to be used
#'   instead of plotting symbol.
#' @param cex The magnification to be used for plotting symbol.
#' @param col A vector with colour of plotting symbols.
#' @param level The confidence level utilised in the plot.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param fixed Deprecated argument (replaced by 'common').
#' @param \dots Graphical arguments as in \code{par} may also be
#'   passed as arguments.
#' 
#' @details
#' A radial plot (Galbraith 1988a,b), also called Galbraith plot, is
#' drawn in the active graphics window. If \code{common} is TRUE, the
#' pooled estimate of the common effect model is plotted. If
#' \code{level} is not NULL, the corresponding confidence limits are
#' drawn.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metabias}}, \code{\link{metabin}},
#'   \code{\link{metagen}}, \code{\link{funnel}}
#' 
#' @references
#' Galbraith RF (1988a):
#' Graphical display of estimates having differing standard errors.
#' \emph{Technometrics},
#' \bold{30}, 271--81
#' 
#' Galbraith RF (1988b):
#' A note on graphical presentation of estimated odds ratios from
#' several clinical trials.
#' \emph{Statistics in Medicine},
#' \bold{7}, 889--94
#' 
#' @keywords hplot
#' 
#' @examples
#' data(Olkin1995)
#' m1 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'   data = Olkin1995, subset = c(41, 47, 51, 59),
#'   studlab = paste(author, year), sm = "RR", method = "I")
#' 
#' # Radial plot
#' #
#' radial(m1, level = 0.95)
#' 
#' @rdname radial
#' @method radial meta
#' @export


radial.meta <- function(x,
                        ##
                        xlim = NULL, ylim = NULL,
                        xlab = "Inverse of standard error",
                        ylab = "Standardised treatment effect (z-score)",
                        ##
                        common = TRUE,
                        ##
                        axes = TRUE,
                        pch = 1, text = NULL, cex = 1, col = NULL,
                        ##
                        level = NULL,
                        warn.deprecated = gs("warn.deprecated"),
                        fixed,
                        ...) {
  
  
  ##
  ##
  ## (1) Check for meta object
  ##
  ##
  chkclass(x, "meta")
  chksuitable(x, "Radial plot", "metamerge", check.mlm = FALSE)
  ##
  x <- updateversion(x)
  ##
  TE <- x$TE
  seTE <- x$seTE
  
  
  ##
  ##
  ## (2) Check arguments
  ##
  ##
  if (!is.null(xlim))
    chknumeric(xlim, length = 2)
  if (!is.null(ylim))
    chknumeric(ylim, length = 2)
  chkchar(xlab, length = 1)
  chkchar(ylab, length = 1)
  ##
  chklogical(warn.deprecated)
  ##
  common <-
    deprecated2(common, missing(common), fixed, missing(fixed),
                warn.deprecated)
  chklogical(common)
  ##
  chklogical(axes)
  if (!is.null(level))
    chklevel(level)
  ##
  if (!is.null(text))
    chklength(text, length(TE),
              text = paste(
                "Argument 'text' must have same length as number of",
                "studies in meta-analysis."))
  
  
  ##
  ##
  ## (3) Produce radial plot
  ##
  ## 
  if (!is.null(x$exclude)) {
    TE <- TE[!x$exclude]
    seTE <- seTE[!x$exclude]
    if (!is.null(text))
      text <- text[!x$exclude]
  }
  ##
  zscore <- TE / seTE
  ##
  if (is.null(xlim)) xlim <- c(0, max(1 / seTE, na.rm = TRUE))
  if (is.null(ylim)) ylim <- c(min(c(-2, zscore), na.rm = TRUE),
                               max(c( 2, zscore), na.rm = TRUE))
  ##
  plot(1 / seTE, zscore,
       xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
       axes = axes, type = "n", ...)

  if (is.null(text)) {
    if (is.null(col))
      points(1 / seTE, zscore, pch = pch, cex = cex)
    else
      points(1 / seTE, zscore, pch = pch, cex = cex, col = col)
  }
  else {
    if (is.null(col))
      text(1 / seTE, zscore, labels = text, cex = cex)
    else
      text(1 / seTE, zscore, labels = text, cex = cex, col = col)
  }

  if (common) {
    lmcomb <- lm(zscore ~ I(1 / seTE) - 1)
    abline(lmcomb, lty = 4)
    if (!is.null(level)) {
      alpha <- 1 - level
      abline(qnorm(1 - alpha / 2), coef(lmcomb), lty = 2)
      abline(qnorm(    alpha / 2), coef(lmcomb), lty = 2)
    }
  }
  
  invisible(NULL)
}





#' @rdname radial
#' @method radial default
#' @export


radial.default <- function(x, y,
                           ##
                           xlim = NULL, ylim = NULL,
                           xlab = "Inverse of standard error",
                           ylab = "Standardised treatment effect (z-score)",
                           ##
                           common = TRUE,
                           ##
                           axes = TRUE,
                           ##
                           pch = 1, text = NULL, cex = 1, col = NULL,
                           ##
                           level = NULL, ...) {
  
  
  ##
  ##
  ## (1) Check essential arguments
  ##
  ##
  chknumeric(x)
  chknumeric(y)
  ##
  chklength(y, length(x), "x")
  
  
  ##
  ##
  ## (2) Do meta-analysis
  ##
  ##
  m <- metagen(x, y, method.tau.ci = "")
  
  
  ##
  ##
  ## (3) Produce radial plot
  ##
  ##
  radial(m, xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab,
         common = common, axes = axes,
         pch = pch, text = text,
         cex = cex, col = col, level = level, ...)
  
  
  invisible(NULL)
}
