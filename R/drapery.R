#' Drapery plot
#' 
#' @description
#' Draw a drapery plot with (scaled) p-value curves for individual
#' studies and meta-analysis estimates.
#' 
#' @aliases drapery
#' 
#' @param x An object of class \code{meta}.
#' @param type A character string indicating whether to plot test
#'   statistics (\code{"zvalue"}) or p-values (\code{"pvalue"}), can
#'   be abbreviated.
#' @param layout A character string for the line layout of individual
#'   studies: \code{"grayscale"}, \code{"equal"}, or
#'   \code{"linewidth"} (see Details), can be abbreviated.
#' @param study.results A logical indicating whether results for
#'   individual studies should be shown in the figure.
#' @param lty.study Line type for individual studies.
#' @param lwd.study Line width for individual studies.
#' @param col.study Colour of lines for individual studies.
#' @param labels A logical or character string indicating whether
#'   study labels should be shown at the top of the drapery plot;
#'   either \code{FALSE}, \code{"id"}, or \code{"studlab"}; see
#'   Details.
#' @param col.labels Colour of study labels.
#' @param cex.labels The magnification for study labels.
#' @param subset.labels A vector specifying which study labels should
#'   be shown in the drapery plot.
#' @param srt.labels A numerical vector or single numeric (between 0
#'   and 90) specifying the angle to rotate study labels; see Details.
#' @param common A logical indicating whether to show result for the
#'   common effect model.
#' @param random A logical indicating whether to show result for the
#'   random effects model.
#' @param lty.common Line type for common effect meta-analysis.
#' @param lwd.common Line width for common effect meta-analysis.
#' @param col.common Colour of lines for common effect meta-analysis.
#' @param lty.random Line type for random effects meta-analysis.
#' @param lwd.random Line width for random effects meta-analysis.
#' @param col.random Colour of lines for random effects meta-analysis.
#' @param prediction A logical indicating whether to show prediction
#'   region.
#' @param col.predict Colour of prediction region
#' @param sign Significance level used to highlight significant values
#'   in curves.
#' @param lty.sign Line type for significant values.
#' @param lwd.sign Line width for significant values.
#' @param col.sign Line colour for significant values.
#' @param alpha Horizonal lines are printed for the specified alpha
#'   values.
#' @param lty.alpha Line type of horizonal lines for alpha values.
#' @param lwd.alpha Line width of horizonal lines for alpha values.
#' @param col.alpha Colour of horizonal lines for alpha values.
#' @param cex.alpha The magnification for the text of the alpha
#' @param col.null.effect Colour of vertical line indicating null
#'   effect.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param pos.legend A character string with position of legend (see
#'   \code{\link{legend}}).
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
#' @param ylim The y limits (min, max) of the plot (ignored if
#'   \code{type = "pvalue"}).
#' @param lwd.max The maximum line width (only considered if argument
#'   \code{layout} is equal to \code{"linewidth"}).
#' @param lwd.study.weight A character string indicating whether to
#'   determine line width for individual studies using weights from
#'   common effect (\code{"common"}) or random effects model
#'   (\code{"random"}), can be abbreviated (only considered if
#'   argument \code{layout} is equal to \code{"linewidth"}).
#' @param at Points at which tick-marks are to be drawn on the x-axis.
#' @param n.grid The number of grid points to calculate the p-value or
#'   test statistic functions.
#' @param mar Physical plot margin, see \code{\link{par}}.
#' @param plot A logical indicating whether to generate a figure.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param fixed Deprecated argument (replaced by 'common').
#' @param lwd.fixed Deprecated argument (replaced by 'lwd.common').
#' @param lty.fixed Deprecated argument (replaced by 'lty.common').
#' @param col.fixed Deprecated argument (replaced by 'col.common').
#' @param \dots Graphical arguments as in \code{par} may also be
#'   passed as arguments.
#' 
#' @details
#' The concept of a p-value function, also called confidence curve,
#' goes back to Birnbaum (1961). A drapery plot, showing p-value
#' functions (or a scaled version based on the corresponding test
#' statistics) for individual studies as well as meta-analysis
#' estimates, is drawn in the active graphics window. Furthermore, a
#' prediction region for a single future study is shown as a shaded
#' area. In contrast to a forest plot, a drapery plot does not provide
#' information for a single confidence level however for any
#' confidence level.
#'
#' Argument \code{type} can be used to either show p-value functions
#' (Birnbaum, 1961) or a scaled version (Infanger, 2019) with test
#' statistics (default).
#' 
#' Argument \code{layout} determines how curves for individual studies
#' are presented:
#' \itemize{
#' \item darker gray tones with increasing precision (\code{layout =
#'   "grayscale"})
#' \item thicker lines with increasing precision (\code{layout =
#'   "linewidth"})
#' \item equal lines (\code{layout = "equal"})
#' }
#'
#' Argument \code{labels} determines how curves of individual studies
#' are labelled:
#' \itemize{
#' \item number of the study in the (unsorted) forest plot / printout
#'   of a meta-analysis (\code{labels = "id"})
#' \item study labels provided by argument \code{studlab} in
#'   meta-analysis functions (\code{labels = "studlab"})
#' \item no study labels (\code{labels = FALSE})
#' }
#' By default, study labels are used (\code{labels = "studlab"}) if no
#' label has more than three characters; otherwise IDs are used
#' (\code{labels = "id"}). The connection between IDs and study labels
#' (among other information) is part of a data frame which is
#' invisibly returned (if argument \code{ study.results = TRUE}).
#'
#' Argument \code{srt.labels} can be used to change the rotation of
#' IDs or study labels. By default, study labels are rotated by +/- 45
#' degrees if at least one study label has more than three characters;
#' otherwise labels are not rotated.
#' 
#' If \code{labels = "studlab"}, labels are rotated by -45 degrees for
#' studies with a treatment estimate below the common effect estimate
#' and otherwise by 45 degrees.
#' 
#' @author Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de}, Guido
#'   Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
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
#' Infanger D and Schmidt-Trucksäss A (2019):
#' P value functions: An underused method to present research results
#' and to promote quantitative reasoning
#' \emph{Statistics in Medicine},
#' \bold{38}, 4189--97
#' 
#' @keywords hplot
#' 
#' @examples
#' data("lungcancer")
#' m1 <- metainc(d.smokers, py.smokers, d.nonsmokers, py.nonsmokers,
#'   data = lungcancer, studlab = study)
#' 
#' # Drapery plot
#' #
#' drapery(m1, xlim = c(0.5, 50))
#' 
#' \dontrun{
#' data(Fleiss1993bin)
#' m2 <- metabin(d.asp, n.asp, d.plac, n.plac,
#'   data = Fleiss1993bin, studlab = paste(study, year),
#'   sm = "OR", random = FALSE)
#'
#' # Produce drapery plot and print data frame with connection between
#' # IDs and study labels
#' #
#' (drapery(m2))
#' 
#' # For studies with a significant effect (p < 0.05), show
#' # study labels and print labels and lines in red
#' #
#' drapery(m2,
#'   labels = "studlab", subset.labels = pval < 0.05,
#'   srt.labels = 0, col.labels = "red",
#'   col.study = ifelse(pval < 0.05, "red", "darkgray"))
#' }
#' 
#' @rdname drapery
#' @export drapery


drapery <- function(x, type = "zvalue", layout = "grayscale",
                    ##
                    study.results = TRUE,
                    lty.study = 1, lwd.study = 1, col.study = "darkgray",
                    ##
                    labels, col.labels = "black", cex.labels = 0.7,
                    subset.labels, srt.labels,
                    ##
                    common = x$common, random = x$random,
                    lty.common = 1, lwd.common = max(3, lwd.study),
                    col.common = "blue",
                    lty.random = 1, lwd.random = lwd.common,
                    col.random = "red",
                    ##
                    sign = NULL,
                    lty.sign = 1, lwd.sign = 1, col.sign = "black",
                    prediction = random, col.predict = "lightblue",
                    ##
                    alpha = if (type == "zvalue") c(0.001, 0.01, 0.05, 0.1)
                            else c(0.01, 0.05, 0.1),
                    lty.alpha = 2, lwd.alpha = 1, col.alpha = "black",
                    cex.alpha = 0.7,
                    col.null.effect = "black",
                    ##
                    legend = TRUE, pos.legend = "topleft",
                    bg = "white", bty = "o",
                    ##
                    backtransf = x$backtransf,
                    xlab, ylab,
                    xlim, ylim,
                    lwd.max = 2.5,
                    lwd.study.weight = if (random) "random" else "common",
                    at = NULL,
                    n.grid = if (type == "zvalue") 10000 else 1000,
                    mar = c(5.1, 4.1, 4.1, 4.1),
                    plot = TRUE,
                    ##
                    warn.deprecated = gs("warn.deprecated"),
                    fixed, lwd.fixed, lty.fixed, col.fixed,
                    ##
                    ...) {
  
  
  pvf <- function(x, TE, seTE)
    2 * pnorm(abs(TE - x) / seTE, lower.tail = FALSE)
  
  
  if (plot) {
    oldpar <- par(mar = mar)
    on.exit(par(oldpar))
  }
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  chkclass(x, "meta")
  x <- updateversion(x)
  ##
  seq.TE <- seq(along = x$TE)
  k.all <- length(seq.TE)
  ##
  mf <- match.call()
  fun <- "drapery"
  ##
  types <- c("zvalue", "pvalue", "surprisal")
  type <- setchar(type, types)
  layout <- setchar(layout, c("grayscale", "linewidth", "equal"))
  ##
  chklogical(study.results)
  ##
  ## Argument lty.study
  ##
  if (!missing(lty.study)) {
    lty.study <- catchvar("lty.study", x, mf)
    chknumeric(lty.study, min = 0, zero = TRUE)
  }
  lty.study <- augment(lty.study, k.all, fun)
  ##
  ## Argument lwd.study
  ##
  missing.lwd.study <- missing(lwd.study)
  ##
  if (!missing.lwd.study) {
    lwd.study <- catchvar("lwd.study", x, mf)
    chknumeric(lwd.study, min = 0, zero = TRUE)
  }
  lwd.study <- augment(lwd.study, k.all, fun)
  ##
  ## Argument col.study
  ##
  missing.col.study <- missing(col.study)
  ##
  if (!missing.col.study)
    col.study <- catchvar("col.study", x, mf)
  col.study <- augment(col.study, k.all, fun)
  ##
  ## Argument labels
  ##
  if (missing(labels) || (is.logical(labels) && labels)) {
    if (max(nchar(x$studlab)) < 4)
      labels <- "studlab"
    else
      labels <- "id"
  }
  else {
    if (is.logical(labels) && !labels)
      labels <- ""
    else
      labels <- setchar(labels, c("", "id", "studlab"))
  }
  ##
  ## Argument col.labels
  ##
  if (!missing(col.labels))
    col.labels <- catchvar("col.labels", x, mf)
  col.labels <- augment(col.labels, k.all, fun)
  ##
  ## Argument cex.labels
  ##
  if (!missing(cex.labels)) {
    cex.labels <- catchvar("cex.labels", x, mf)
    chknumeric(cex.labels, min = 0, zero = TRUE)
  }
  cex.labels <- augment(cex.labels, k.all, fun)
  ##
  ## Argument subset.labels
  ##
  if (missing(subset.labels))
    subset.labels <- rep(TRUE, k.all)
  else {
    subset.labels <- catchvar("subset.labels", x, mf)
    chklength(subset.labels, k.all, fun)
    ##
    if (is.numeric(subset.labels)) {
      if (any(!(subset.labels %in% seq.TE)))
        stop("Numerical input for argument 'subset.labels' ",
             "must contain integers between 1 and ",
             k.all, ".",
             call. = FALSE)
      sel.labels <- rep(FALSE, k.all)
      sel.labels[subset.labels] <- TRUE
    }
    else if (!is.logical(subset.labels))
      stop("Argument 'subset.labels' must be a logical or ",
           "numeric vector.",
           call. = FALSE)
  }
  ##
  ## Argument srt.labels
  ##
  if (missing(srt.labels)) {
    if (labels == "studlab" & max(nchar(x$studlab)) > 3)
      srt.labels <- ifelse(x$TE < x$TE.common, 45, -45)
    else
      srt.labels <- rep(0, k.all)
  }
  else {
    srt.labels <- catchvar("srt.labels", x, mf)
    chknumeric(srt.labels, min = -90, max = 90)
  }
  srt.labels <- augment(srt.labels, k.all, fun)
  ##
  ## Other arguments
  ##
  common <-
    deprecated2(common, missing(common), fixed, missing(fixed),
                warn.deprecated)
  chklogical(common)
  chklogical(random)
  chklogical(prediction)
  ##
  lwd.common <-
    deprecated2(lwd.common, missing(lwd.common),
                lwd.fixed, missing(lwd.fixed),
                warn.deprecated)
  chknumeric(lwd.common, length = 1)
  ##
  lty.common <-
    deprecated2(lty.common, missing(lty.common),
                lty.fixed, missing(lty.fixed),
                warn.deprecated)
  chknumeric(lty.common, length = 1)
  ##
  col.common <-
    deprecated2(col.common, missing(col.common),
                col.fixed, missing(col.fixed),
                warn.deprecated)
  chklength(col.common, 1,
            text = "Argument 'col.common' must be of length 1.")
  chklength(col.random, 1,
            text = "Argument 'col.random' must be of length 1.")
  ##
  if (!is.null(sign)) {
    chklevel(sign, length = 0)
    if (type == "zvalue")
      sign <- qnorm(sign / 2)
    ##
    chknumeric(lty.sign, min = 0, zero = TRUE, length = 1)
    chknumeric(lwd.sign, min = 0, zero = TRUE, length = 1)
  }
  ##
  if (!all(is.na(alpha)))
    chklevel(alpha, length = 0)
  chknumeric(lty.alpha, min = 0, zero = TRUE, length = 1)
  chknumeric(lwd.alpha, min = 0, zero = TRUE, length = 1)
  chknumeric(cex.alpha, min = 0, zero = TRUE, length = 1)
  ##
  chklogical(legend)
  bty <- setchar(bty, c("o", "n"))
  chklogical(backtransf)
  ##
  if (missing(xlab)) {
    xlab <- ""
    xlab <- xlab_meta(x$sm, backtransf = backtransf)
  }
  else
    chkchar(xlab, length = 1)
  ##
  if (missing(ylab))
    ylab <- c("Test statistic",
              "P-value",
              "Binary S-value / surprisal")[charmatch(type, types)]
  else
    chkchar(ylab, length = 1)
  ##
  chknumeric(lwd.max, min = 0, zero = TRUE, length = 1)
  lwd.study.weight <- setchar(lwd.study.weight, c("random", "common"))
  ##
  chknumeric(n.grid, min = 0, zero = TRUE, length = 1)
  chklogical(plot)
  ##
  ## Additional check
  ##
  if (!any(c(study.results, common, random, prediction))) {
    warning("Nothing selected to show in drapery plot.")
    return(invisible(NULL))
  }
  
  
  ##
  ##
  ## (2) Set variables for plotting
  ##
  ##
  if (missing(ylim))
    ylim <-
      if (type == "pvalue") c(0, 1)
      else if (type == "zvalue") c(qnorm(0.0005), 0)
      else c(0, 12)
  else
    if (type == "pvalue") {
      warning("Argument 'ylim' ignored for p-value function plot." )
      ylim <- c(0, 1)
    }
  ##
  x.backtransf <- is_relative_effect(x$sm) & backtransf
  ##
  missing.xlim <- missing(xlim)
  ##
  if (missing.xlim) {
    mn <- min(x$lower, na.rm = TRUE)
    mx <- max(x$upper, na.rm = TRUE)
    ##
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
  grid <- c(seq(mn, mx, length.out = n.grid), x$TE, x$TE.common, x$TE.random)
  ##
  if (type == "zvalue")
    grid <- c(grid,
              x$TE + ylim[1] * x$seTE,
              x$TE - ylim[1] * x$seTE,              
              x$TE.common + ylim[1] * x$seTE.common,
              x$TE.common - ylim[1] * x$seTE.common,
              x$TE.random + ylim[1] * x$seTE.random,
              x$TE.random - ylim[1] * x$seTE.random)
  ##
  grid <- sort(grid)
  ##
  if (x.backtransf) {
    mn <- exp(mn)
    mx <- exp(mx)
    x.grid <- exp(grid)
    xlim <- exp(xlim)
    null.effect <- exp(x$null.effect)
  }
  else {
    x.grid <- grid
    null.effect <- x$null.effect
  }
  ##
  o <- order(x$seTE)
  ##
  TE <- x$TE[o]
  seTE <- x$seTE[o]
  lower <- x$lower[o]
  upper <- x$upper[o]
  pval <- x$pval[o]
  statistic <- x$statistic[o]
  w.common <- x$w.common[o]
  w.random <- x$w.random[o]
  seq.TE <- seq.TE[o]
  studlab <- x$studlab[o]
  ##
  lty.study <- lty.study[o]
  lwd.study <- lwd.study[o]
  col.study <- col.study[o]
  ##
  col.labels <- col.labels[o]
  cex.labels <- cex.labels[o]
  subset.labels <- subset.labels[o]
  srt.labels <- srt.labels[o]
  ##
  TE.common <- x$TE.common
  seTE.common <- x$seTE.common
  TE.random <- x$TE.random
  seTE.random <- x$seTE.random
  seTE.predict <- x$seTE.predict
  t.quantile <- (x$upper.predict - x$lower.predict) / seTE.predict / 2
  ##
  if (layout == "grayscale" & missing.col.study)
    col.study <- gray.colors(k.all)
  ##
  else if (layout == "linewidth" & missing.lwd.study) {
    if (lwd.study.weight == "common")
      lwd.study <- lwd.max * w.common / max(w.common)
    else
      lwd.study <- lwd.max * w.random / max(w.random)
  }
  ##  
  y.common <- pvf(grid, TE.common, seTE.common)
  y.random <- pvf(grid, TE.random, seTE.random)
  y.predict <- pvf(grid, TE.random, t.quantile / qnorm(0.975) * seTE.predict)
  y.alpha <- alpha
  ##
  if (type == "zvalue") {
    y.common <- qnorm(y.common / 2)
    y.random <- qnorm(y.random / 2)
    y.predict <- qnorm(y.predict / 2)
    y.alpha <- qnorm(y.alpha / 2)
  }
  else if (type == "surprisal") {
    y.common <- -log(y.common, base = 2)
    y.random <- -log(y.random, base = 2)
    y.predict <- -log(y.predict, base = 2)
    y.alpha <- -log(y.alpha, base = 2)
  }
  ##
  sel.common <- y.common >= min(ylim)
  sel.random <- y.random >= min(ylim)
  sel.predict <- y.predict >= min(ylim)
  ##
  if (!study.results & missing.xlim) {
    if (prediction) {
      xvals <- x.grid[sel.predict]
      if (common)
        xvals <- c(xvals, x.grid[sel.common])
      xlim <- range(xvals)
    }
    else {
      if (common & random)
        xvals <- c(x.grid[sel.common], x.grid[sel.random])
      else if (!common)
        xvals <- x.grid[sel.random]
      else
        xvals <- x.grid[sel.common]
      xlim <- range(xvals)
    }
  }
  
  
  ##
  ##
  ## (3) Generate drapery plot
  ##
  ##
  if (plot) {
    plot(x.grid, y.common,
         xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
         type = "n", las = 1, if (x.backtransf) log = "x" else log = "",
         axes = FALSE, ...)
    axis(1, at = at)
    axis(2)
    box()
    ##
    ## Add prediction region
    ##
    if (prediction) {
      polygon(x.grid[sel.predict], y.predict[sel.predict],
              col = col.predict, border = NA)
      ##
      if (type != "surprisal")
        polygon(c(min(x.grid[sel.predict]),
                  max(x.grid[sel.predict]),
                  max(x.grid[sel.predict])),
                c(min(y.predict[sel.predict]),
                  min(y.predict[sel.predict]),
                  tail(y.predict[sel.predict], n = 1)),
                col = col.predict, border = NA)
      else
        polygon(c(max(x.grid[sel.predict]),
                  min(x.grid[sel.predict]),
                  min(x.grid[sel.predict])),
                c(max(y.predict[sel.predict]),
                  max(y.predict[sel.predict]),
                  head(y.predict[sel.predict], n = 1)),
                col = col.predict, border = NA)
      ##
      if (all(y.predict > min(ylim)) & type != "surprisal") {
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
    ## Add studies
    ##
    if (study.results) {
      for (i in k.all:1) {
        y.i <- pvf(grid, TE[i], seTE[i])
        if (type == "zvalue")
          y.i <- qnorm(y.i / 2)
        else if (type == "surprisal")
          y.i <- -log(y.i, base = 2)
        sel.i <- y.i >= min(ylim)
        lines(x.grid[sel.i], y.i[sel.i],
              lty = lty.study[i], lwd = lwd.study[i], col = col.study[i])
        ##
        if (!is.null(sign)) {
          sel.sign.i <- sel.i & y.i < sign
          seq.i <- seq(along = y.i)[sel.sign.i]
          sel.seq.i.u <- seq.i[c(FALSE, diff(seq.i) != 1)]
          sel.seq.i.l <- seq.i[c(diff(seq.i) != 1, FALSE)]
          sel.sign.i.l <- sel.sign.i & seq(along = sel.sign.i) <= sel.seq.i.l
          sel.sign.i.u <- sel.sign.i & seq(along = sel.sign.i) >= sel.seq.i.u
          ##
          lines(x.grid[sel.sign.i.l], y.i[sel.sign.i.l],
                lty = lty.sign, lwd = lwd.sign, col = col.sign)
          lines(x.grid[sel.sign.i.u], y.i[sel.sign.i.u],
                lty = lty.sign, lwd = lwd.sign, col = col.sign)
        }
        ##
        if (subset.labels[i]) {
          if (type != "surprisal") {
            yval <- ylim[2] + 0.005 * abs(ylim[1] - ylim[2])
            y.adj <- 0
          }
          else {
            yval <- ylim[1] - 0.005 * abs(ylim[1] - ylim[2])
            y.adj <- 1
          }
          ##
          if (labels != "") {
            if (srt.labels[i] == 0)
              adj.i <- c(0.5, y.adj)
            else
              adj.i <- if (srt.labels[i] > 0) c(1, y.adj) else c(0, y.adj)
            ##
            text(if (x.backtransf) exp(TE[i]) else TE[i],
                 yval,
                 labels = if (labels == "id") seq.TE[i] else studlab[i],
                 col = col.labels[i], cex = cex.labels[i],
                 srt = srt.labels[i], adj = adj.i)
          }
        }
      }
    }
    ##
    ## Add common effect lines
    ##
    if (common)
      lines(x.grid[sel.common], y.common[sel.common],
            lty = lty.common, lwd = lwd.common, col = col.common)
    ##
    ## Add random effects lines
    ##
    if (random)
      lines(x.grid[sel.random], y.random[sel.random],
            lty = lty.random, lwd = lwd.random, col = col.random)
    ##
    ## Add null effect line
    ##
    abline(v = null.effect, col = col.null.effect)
    ##
    ## Add p-value lines
    ##
    if (type != "surprisal")
      for (alpha.i in y.alpha)
        abline(h = alpha.i, lty = lty.alpha, lwd = lwd.alpha, col = col.alpha)
    ##
    ## Add text for alpha levels
    ##
    if (type != "surprisal")
      for (i in seq(along = y.alpha))
        text(mn, y.alpha[i] + (ylim[2] - ylim[1]) / 200,
             paste("p =", alpha[i]),
             cex = cex.alpha, col = col.alpha, adj = c(0, 0))
    ##
    ## Add text for confidence levels
    ##
    if (type != "surprisal")
      for (i in seq(along = y.alpha))
        text(mx, y.alpha[i] + (ylim[2] - ylim[1]) / 200,
             paste0(100 * (1 - alpha[i]), "%-CI"),
             cex = cex.alpha, col = col.alpha, adj = c(1, 0))
    ##
    ## Add legend
    ##
    if (legend & any(c(common, random, prediction)))
      legend(pos.legend,
             col = c(if (common) col.common,
                     if (random) col.random,
                     if (prediction) col.predict),
             lwd = c(if (common) lwd.common,
                     if (random) lwd.random,
                     if (prediction) lwd.random),
             bty = bty, bg = bg,
             c(if (common) gs("text.common"),
               if (random) gs("text.random"),
               if (prediction) "Range of prediction"))
    ##
    ## Add y-axis on right side
    ##
    if (type == "pvalue") {
      par(new = TRUE)
      plot(x.grid, y.common,
           xlab = xlab, ylab = ylab,
           xlim = xlim, ylim = rev(ylim),
           type = "n", las = 1, if (x.backtransf) log = "x" else log = "",
           axes = FALSE,
           ...)
      axis(4, las = 1)
    }
    else if (type == "zvalue") {
      axis(4, las = 1,
           at = qnorm(c(0.001, 0.01, 0.05) / 2),
           labels = 1 - c(0.001, 0.01, 0.05))
      axis(4, las = 1,
           at = qnorm(seq(0.2, 1, by = 0.2) / 2),
           labels = format(seq(0.8, 0, by = -0.2)))
    }
    ##
    if (type %in% c("pvalue", "zvalue"))
      mtext("Confidence level", 4, line = 3)
  }
  
  
  if (study.results) {
    res <- data.frame(ID = seq.TE, studlab = studlab,
                      lty.study, lwd.study, col.study,
                      col.labels, cex.labels,
                      subset.labels, srt.labels,
                      TE, seTE, lower, upper, statistic, pval,
                      w.common, w.random,
                      stringsAsFactors = FALSE)
    ##
    res <- res[order(res$ID), , drop = FALSE]
    ##
    attr(res, "xlim") <- xlim
    attr(res, "ylim") <- ylim
  }
  else
    res <- NULL
  
  
  invisible(res)
}
