#' Bubble plot to display the result of a meta-regression
#' 
#' @description
#' Draw a bubble plot to display the result of a meta-regression.
#' 
#' @aliases bubble bubble.metareg
#' 
#' @param x An object of class \code{metareg}.
#' @param xlim The x limits (min,max) of the plot.
#' @param ylim The y limits (min,max) of the plot.
#' @param xlab A label for the x-axis.
#' @param ylab A label for the y-axis.
#' @param cex The magnification to be used for plotting symbols.
#' @param min.cex Minimal magnification for plotting symbols.
#' @param max.cex Maximal magnification for plotting symbols.
#' @param pch The plotting symbol used for individual studies.
#' @param col A vector with colour of plotting symbols.
#' @param bg A vector with background colour of plotting symbols (only
#'   used if \code{pch} in \code{21:25}).
#' @param lty The line type for the meta-regression line.
#' @param lwd The line width for the meta-regression line.
#' @param col.line Colour for the meta-regression line.
#' @param studlab A logical indicating whether study labels should be
#'   printed in the graph. A vector with study labels can also be
#'   provided (must be of same length as the numer of studies in the
#'   meta-analysis then).
#' @param cex.studlab The magnification for study labels.
#' @param pos A position specifier for study labels (see
#'   \code{\link{text}}).
#' @param offset Offset for study labels (see \code{\link{text}}).
#' @param regline A logical indicating whether a regression line
#'   should be added to the bubble plot.
#' @param axes A logical indicating whether axes should be printed.
#' @param box A logical indicating whether a box should be printed.
#' @param \dots Graphical arguments as in \code{par} may also be
#'   passed as arguments.
#' 
#' @details
#' A bubble plot can be used to display the result of a
#' meta-regression. It is a scatter plot with the treatment effect for
#' each study on the y-axis and the covariate used in the
#' meta-regression on the x-axis. Typically, the size of the plotting
#' symbol is inversely proportional to the variance of the estimated
#' treatment effect (Thompson & Higgins, 2002).
#' 
#' Argument \code{cex} specifies the plotting size for each individual
#' study.  If this argument is missing the weights from the
#' meta-regression model will be used (which typically is a random
#' effects model). Use \code{weight="fixed"} in order to utilise
#' weights from a fixed effect model to define the size of the plotted
#' symbols (even for a random effects meta-regression). If a vector
#' with individual study weights is provided, the length of this
#' vector must be of the same length as the number of studies.
#' 
#' Arguments \code{min.cex} and \code{max.cex} can be used to define
#' the size of the smallest and largest plotting symbol. The plotting
#' size of the most precise study is set to \code{max.cex} whereas the
#' plotting size of all studies with a plotting size smaller than
#' \code{min.cex} will be set to \code{min.cex}.
#' 
#' For a meta-regression with more than one covariate. Only a scatter plot of
#' the first covariate in the regression model is shown. In this case the
#' effect of the first covariate adjusted for other covariates in the
#' meta-regression model is shown.
#' 
#' For a factor or categorial covariate separate bubble plots for each
#' group compared to the baseline group are plotted.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{metagen}}, \code{\link{metainf}}
#' 
#' @references
#' Thompson SG, Higgins JP (2002):
#' How should meta-regression analyses be undertaken and interpreted?
#' \emph{Statistics in Medicine},
#' \bold{21}, 1559--73
#'
#' @keywords hplot
#' 
#' @examples
#' data(Fleiss93cont)
#' 
#' # Add some (fictitious) grouping variables:
#' Fleiss93cont$age <- c(55, 65, 52, 65, 58)
#' Fleiss93cont$region <- c("Europe", "Europe", "Asia", "Asia", "Europe")
#' 
#' m1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c,
#'                data = Fleiss93cont, sm = "MD")
#' 
#' mr1 <- metareg(m1, region)
#' mr1
#' 
#' bubble(mr1)
#' bubble(mr1, lwd = 2, col.line = "blue")
#' 
#' mr2 <- metareg(m1, age)
#' mr2
#' 
#' bubble(mr2, lwd = 2, col.line = "blue", xlim = c(50, 70))
#' bubble(mr2, lwd = 2, col.line = "blue", xlim = c(50, 70), cex = "fixed")
#' 
#' # Do not print regression line
#' #
#' bubble(mr2, lwd = 2, col.line = "blue", xlim = c(50, 70), regline = FALSE)
#' 
#' @rdname bubble
#' @export bubble


bubble <- function(x, ...) 
  UseMethod("bubble")





#' @rdname bubble
#' @method bubble metareg
#' @export
#' @export bubble.metareg


bubble.metareg <- function(x,
                           xlim, ylim,
                           xlab, ylab,
                           cex, min.cex = 0.5, max.cex = 5,
                           pch = 21, col = "black", bg = "darkgray",
                           lty = 1, lwd = 1, col.line = "black",
                           studlab = FALSE, cex.studlab = 0.8,
                           pos = 2, offset = 0.5,
                           regline = TRUE,
                           axes = TRUE, box = TRUE, ...) {


  ##
  ##
  ## (1) Check for meta object
  ##
  ##
  chkclass(x, "metareg")


  m0 <- x$.meta$x
  method.tau0 <- x$.meta$method.tau
  ##
  if (method.tau0 != "FE" & (method.tau0 != m0$method.tau))
    m1 <- update(m0, method.tau = method.tau0)
  else
    m1 <- m0
  ##
  TE <- m1$TE
  sm <- m1$sm


  if (is.logical(studlab)) {
    if (studlab)
      studlab <- m1$studlab
    else
      studlab <- rep("", length(TE))
  }
  else {
    studlab <- as.character(studlab)
    if (length(studlab) != length(TE))
      stop("Length of argument 'studlab' must be the same as number of studies in meta-analysis.")
  }


  charform <- as.character(x$.meta$formula)[2]
  splitform <- strsplit(charform, " ")[[1]]
  covar.name <- splitform[1]
  if (covar.name == "1" | covar.name == "-1")
    covar.name <- splitform[3]
  ##
  covar.names <- names(coef(x))
  covar.names.without.intrcpt <- covar.names[covar.names != "intrcpt"]
  if (length(covar.names.without.intrcpt) == 0) {
    warning("No covariate in meta-regression.")
    return(invisible(NULL))
  }
  ##
  nointrcpt <- ifelse("intrcpt" %in% covar.names, FALSE, TRUE)
  ##
  if (covar.name == ".byvar")
    covar.name <- x$.meta$x$bylab
  ##
  if (length(covar.names.without.intrcpt) > 1) {
    warning(paste("Only first covariate in meta-regression ",
                  "('", covar.name, "') considered in bubble plot. No regression line plotted.",
                  sep = "")
            )
    regline <- FALSE
    if (missing(xlab))
      xlab <- paste("Covariate ", covar.name,
                    " (meta-regression: ", charform, ")", sep = "")
  }
  else
    if (missing(xlab))
      xlab = paste("Covariate", covar.name)
  ##
  if (covar.name %in% names(x$.meta$x$data))
    covar <- x$.meta$x$data[[covar.name]]
  else if (".byvar" %in% names(x$.meta$x$data))
    covar <- x$.meta$x$data[[".byvar"]]
  else
    covar <- get(covar.name)
  ##
  if (!is.null(x$.meta$x$subset))
    covar <- covar[x$.meta$x$subset]
  ##
  if (is.character(covar))
    covar <- as.factor(covar)
  ##
  if (is.factor(covar)) {
    levs <- levels(covar)
    xs <- as.numeric(covar) - 1
    at.x <- sort(unique(xs))
    if (missing(xlim))
      xlim <- c(-0.5, 1.5)
    if (missing(xlab))
      xlab <- ""
  }
  else {
    xs <- covar
    if (missing(xlim))
      xlim <- range(xs, na.rm = TRUE)
  }
  ##
  alpha <- ifelse(nointrcpt, 0, coef(x)["intrcpt"])
  ##
  if (covar.name %in% names(coef(x)))
    beta <- coef(x)[covar.name]
  else
    beta <- coef(x)[".byvar"]


  ys <- TE
  ##
  if (missing(ylim))
    ylim <- range(ys)
  ##
  if (missing(ylab))
    ylab <- paste("Treatment effect (",
                  tolower(xlab(sm, backtransf = FALSE)),
                  ")", sep = "")


  missing.cex <- missing(cex)
  fixed.cex   <- !missing.cex && (length(cex) == 1 & all(cex == "fixed"))
  ##
  if (missing.cex)
    if (method.tau0 == "FE")
      cex <- m1$w.fixed
    else
      cex <- m1$w.random
  else if (fixed.cex)
    cex <- m1$w.fixed
  ##
  if (length(cex) != length(TE))
    stop("Length of argument 'cex' must be the same as number of studies in meta-analysis.")
  ##
  if (missing.cex | fixed.cex) {
    cexs <- max.cex*(cex / max(cex))
    cexs[cexs < min.cex] <- min.cex
  }
  else
    cexs <- cex


  ##
  ## Generate bubble plot
  ##
  if (is.factor(covar)) {
    for (i in 2:length(levs)) {
      sel <- xs %in% c(0, i - 1)
      xs.i <- xs[sel]
      xs.i[xs.i > 0] <- 1
      ys.i <- ys[sel]
      studlab.i <- studlab[sel]
      cexs.i <- cexs[sel]
      ##
      plot(xs.i, ys.i,
       pch = pch, cex = cexs.i,
       xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
       type = "n", axes = FALSE, ...)
      ##
      ## Add regression line
      ##
      if (regline)
        lines(c(0, 1),
              c(alpha, alpha + coef(x)[i]),
              lty = lty, lwd = lwd, col = col.line)
      ##
      points(xs.i, ys.i, cex = cexs.i, pch = pch, col = col, bg = bg)
      ##
      ## x-axis
      ##
      if (axes)
        axis(1, at = 0:1, labels = levs[c(1, i)], ...)
      ##
      ## y-axis
      ##
      if (axes)
        axis(2, ...)
      ##
      text(xs.i, ys.i, labels = studlab.i, cex = cex.studlab,
           pos = pos, offset = offset)
      ##
      if (box)
        box()
    }
  }
  else {
    plot(xs, ys,
         pch = pch, cex = cexs,
         xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
         type = "n", axes = FALSE, ...)
    ##
    ## Add regression line
    ##
    if (regline)
      abline(alpha, beta,
             lty = lty, lwd = lwd, col = col.line)
    ##
    points(xs, ys, cex = cexs, pch = pch, col = col, bg = bg)
    ##
    ## x-axis
    ##
    if (axes)
      axis(1, ...)
    ##
    ## y-axis
    ##
    if (axes)
      axis(2, ...)
    ##
    text(xs, ys, labels = studlab, cex = cex.studlab,
         pos = pos, offset = offset)
    ##
    if (box)
      box()
  }


  invisible(NULL)
}
