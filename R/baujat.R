#' Baujat plot to explore heterogeneity in meta-analysis
#' 
#' @description
#' Draw a Baujat plot to explore heterogeneity in meta-analysis.
#' 
#' @aliases baujat baujat.meta
#' 
#' @param x An object of class \code{meta}.
#' @param yscale Scaling factor for values on y-axis.
#' @param xlim The x limits (min,max) of the plot.
#' @param ylim The y limits (min,max) of the plot.
#' @param xlab A label for the x-axis.
#' @param ylab A label for the y-axis.
#' @param pch The plotting symbol used for individual studies.
#' @param cex The magnification to be used for plotting symbol.
#' @param col A vector with colour of plotting symbols.
#' @param bg A vector with background colour of plotting symbols (only
#'   used if \code{pch} in \code{21:25}).
#' @param studlab A logical indicating whether study labels should be
#'   printed in the graph. A vector with study labels can also be
#'   provided (must be of same length as \code{x$TE} then).
#' @param cex.studlab The magnification for study labels.
#' @param pos.studlab Position of study labels, see argument
#'   \code{pos} in \code{\link{text}}.
#' @param offset Offset for study labels (see \code{\link{text}}).
#' @param xmin A numeric specifying minimal value to print study
#'   labels (on x-axis).
#' @param ymin A numeric specifying minimal value to print study
#'   labels (on y-axis).
#' @param grid A logical indicating whether a grid is printed in the
#'   plot.
#' @param col.grid Colour for grid lines.
#' @param lty.grid The line type for grid lines.
#' @param lwd.grid The line width for grid lines.
#' @param pty A character specifying type of plot region (see
#'   \code{\link{par}}).
#' @param pooled A character string indicating whether a common effect
#'   or random effects model is used for pooling. Either missing (see
#'   Details), \code{"common"} or \code{"random"}, can be abbreviated.
#' @param \dots Graphical arguments as in \code{par} may also be
#'   passed as arguments.
#'
#' @details
#' Baujat et al. (2002) introduced a scatter plot to explore
#' heterogeneity in meta-analysis. On the x-axis the contribution of
#' each study to the overall heterogeneity statistic (see list object
#' \code{Q} of the meta-analysis object \code{x}) is plotted. On the
#' y-axis the standardised difference of the overall treatment effect
#' with and without each study is plotted; this quantity describes the
#' influence of each study on the overal treatment effect.
#' 
#' Information from object \code{x} is utilised if argument
#' \code{pooled} is missing. A common effect model is assumed
#' (\code{pooled="common"}) if argument \code{x$common} is
#' \code{TRUE}; a random effects model is assumed
#' (\code{pooled="random"}) if argument \code{x$random} is
#' \code{TRUE} and \code{x$common} is \code{FALSE}.
#' 
#' Internally, the \code{\link{metainf}} function is used to calculate
#' the values on the y-axis.
#' 
#' @return
#' A data.frame with the following variables:
#' \item{x}{Coordinate on x-axis (contribution to heterogeneity
#'   statistic)}
#' \item{y}{Coordinate on y-axis (influence on overall treatment
#'   effect)}
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{metagen}}, \code{\link{metainf}}
#' 
#' @references
#' Baujat B, Mahé C, Pignon JP, Hill C (2002):
#' A graphical method for exploring heterogeneity in meta-analyses:
#' Application to a meta-analysis of 65 trials.
#' \emph{Statistics in Medicine},
#' \bold{30}, 2641--52
#' 
#' @keywords hplot
#' 
#' @examples
#' data(Olkin1995)
#'
#' # Only consider first ten studies
#' m1 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'   data = Olkin1995, sm = "OR", method = "I", studlab = paste(author, year),
#'   subset = 1:10)
#' 
#' # Generate Baujat plot
#' baujat(m1)
#'
#' \dontrun{
#' m1 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'   data = Olkin1995, sm = "OR", method = "I", studlab = paste(author, year))
#' 
#' # Do not print study labels if the x-value is smaller than 4 and
#' # the y-value is smaller than 1
#' baujat(m1, yscale = 10, xmin = 4, ymin = 1)
#' 
#' # Change position of study labels
#' baujat(m1, yscale = 10, xmin = 4, ymin = 1,
#'        pos = 1, xlim = c(0, 6.5))
#' 
#' # Generate Baujat plot and assign x- and y- coordinates to R object
#' # b1
#' b1 <- baujat(m1)
#' 
#' # Calculate overall heterogeneity statistic
#' sum(b1$x)
#' m1$Q
#' }
#' 
#' @method baujat meta
#' @export


baujat.meta <- function(x,
                        yscale = 1,
                        xlim, ylim,
                        xlab = "Contribution to overall heterogeneity",
                        ylab = "Influence on overall result",
                        pch = 21, cex = 1, col = "black", bg = "darkgray",
                        studlab = TRUE, cex.studlab = 0.8,
                        pos.studlab, offset = 0.5,
                        xmin = 0, ymin = 0,
                        grid = TRUE, col.grid = "lightgray",
                        lty.grid = "dotted", lwd.grid = par("lwd"),
                        pty = "s",
                        pooled,
                        ...) {
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(x, "meta")
  chksuitable(x, "Baujat plot",
              c("trimfill", "metacum", "metainf", "metamerge", "netpairwise"))
  ##
  x <- updateversion(x)
  
  
  ##
  ##
  ## (2) Check arguments
  ##
  ##
  chknumeric(yscale)
  chknumeric(cex)
  chknumeric(cex.studlab)
  missing.pos.studlab <- missing(pos.studlab)
  if (!missing.pos.studlab) {
    pos.studlab <- as.numeric(setchar(pos.studlab, as.character(1:4)))
    chknumeric(pos.studlab, min = 1, max = 4)
  }
  ##
  if (!missing(pooled)) {
    pooled <- setchar(pooled, c("common", "random", "fixed"))
    pooled[pooled == "fixed"] <- "common"
  }
  else
    if (!x$common & x$random)
      pooled <- "random"
    else
      pooled <- "common"
  chknumeric(offset)
  chknumeric(xmin)
  chknumeric(ymin)
  chklogical(grid)
  chknumeric(lwd.grid)
  
  
  oldpar <- par(pty = pty)
  on.exit(par(oldpar))
  
  
  TE <- x$TE
  seTE <- x$seTE
  ##
  if (pooled == "random")
    TE.s <- x$TE.random
  else
    TE.s <- x$TE.common
  ##
  k <- x$k
  ##
  if (is.logical(studlab)) {
    if (studlab)
      labels.studlab <- x$studlab
    else
      labels.studlab <- rep("", length(TE))
  }
  else if (is.numeric(studlab)) {
    labels.studlab <- x$studlab
  }
  else {
    labels.studlab <- as.character(studlab)
    if (length(labels.studlab) != length(TE))
      stop("Length of argument 'studlab' must be the same as ",
           "number of studies in meta-analysis.")
  }
  ##
  if (!missing.pos.studlab) {
    if (length(pos.studlab) > 1)
      chklength(pos.studlab, TE,
                text = paste("Length of argument 'pos.studlab' must be",
                             "the same as number of studies",
                             "in meta-analysis."))
    else
      pos.studlab <- rep(pos.studlab, length(TE))
  }
  
  
  m.inf <- metainf(x, pooled = pooled)
  TE.inf <- m.inf$TE
  seTE.inf <- m.inf$seTE
  ##
  ys <- (TE.inf - TE.s)^2 / seTE.inf^2
  ys <- ys * yscale
  ##  
  xs <- (TE - TE.s)^2 / seTE^2
  ##
  if (!is.null(x$exclude))
    xs[x$exclude] <- NA
  
  
  if (missing(xlim))
    xlim <- c(0, max(xs, na.rm = TRUE))
  ##
  if (missing(ylim))
    ylim <- c(0, max(ys, na.rm = TRUE))
  
  
  ##
  ## Do not print labels for studies with x and/or y values below
  ## limits
  ##
  if (!missing(xmin) & !missing(ymin))
    labels.studlab[xs < xmin & ys < ymin] <- ""
  else if (!missing(xmin) & missing(ymin))
    labels.studlab[xs < xmin] <- ""
  else if (missing(xmin) & !missing(ymin))
    labels.studlab[ys < ymin] <- ""
  ##
  if (is.numeric(studlab)) {
    if (length(studlab) == 1)
      labels.studlab[xs < studlab] <- ""
    else if (length(studlab) == 2) {
      labels.studlab[xs < studlab[1]] <- ""
      labels.studlab[ys < studlab[2]] <- ""
    }
  }
  
  
  plot(xs, ys,
       xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       type = "n")
  ##
  if (grid)
    grid(col = col.grid, lty = lty.grid, lwd = lwd.grid)
  ##
  points(xs, ys, pch = pch, cex = cex, col = col, bg = bg)
  ##
  if (missing.pos.studlab)
    pos.studlab <- ifelse(xs > mean(xlim), 2, 4)
  text(xs, ys, labels = labels.studlab, cex = cex.studlab,
       pos = pos.studlab, offset = offset)
  
  
  res <- data.frame(studlab = x$studlab, x = xs, y = ys)
  attr(res, "pooled") <- pooled
  
  
  invisible(res)
}
