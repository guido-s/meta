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
#' @param pch The plotting symbol(s) used for individual studies.
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
#' @param pos.studlab Position of study labels, see argument
#'   \code{pos} in \code{\link{text}}.
#' @param offset Offset for study labels (see \code{\link{text}}).
#' @param regline A logical indicating whether a regression line
#'   should be added to the bubble plot.
#' @param col.line Colour for the meta-regression line.
#' @param backtransf A logical indicating whether results for relative
#'   summary measures (argument \code{sm} equal to \code{"OR"},
#'   \code{"RR"}, \code{"HR"}, or \code{"IRR"}) should be back
#'   transformed. If \code{backtransf=TRUE}, results for
#'   \code{sm="OR"} are printed as odds ratios rather than log odds
#'   ratios, for example.
#' @param ref A numerical giving the reference value to be plotted as
#'   a line in the bubble plot. No reference line is plotted if
#'   argument \code{ref} is equal to \code{NA}.
#' @param col.ref Colour of the reference line.
#' @param lty.ref The line type for the reference line.
#' @param lwd.ref The line width for the reference line.
#' @param pscale A numeric giving scaling factor for printing of
#'   probabilities.
#' @param irscale A numeric defining a scaling factor for printing of
#'   incidence rates.
#' @param axes Either a logical or a character string equal to \code{"x"},
#'   \code{"y"} or \code{"xy"} indicating whether x- and y-axis should be
#'   printed.
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
#' study. If this argument is missing the weights from the
#' meta-regression model will be used (which typically is a random
#' effects model). Use \code{cex="common"} in order to utilise weights
#' from a common effect model to define the size of the plotted
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
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
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
#' data(Fleiss1993cont)
#' 
#' # Add some (fictitious) grouping variables:
#' Fleiss1993cont$age <- c(55, 65, 52, 65, 58)
#' Fleiss1993cont$region <- c("Europe", "Europe", "Asia", "Asia", "Europe")
#' 
#' m1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, sm = "SMD")
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
#' bubble(mr2, lwd = 2, col.line = "blue", xlim = c(50, 70), cex = "common")
#' 
#' # Do not print regression line
#' #
#' bubble(mr2, lwd = 2, col.line = "blue", xlim = c(50, 70), regline = FALSE)
#' 
#' @rdname bubble.metareg
#' @method bubble metareg
#' @export


bubble.metareg <- function(x,
                           xlim, ylim,
                           xlab, ylab,
                           cex, min.cex = 0.5, max.cex = 5,
                           pch = 21, col = "black", bg = "darkgray",
                           lty = 1, lwd = 1, col.line = "black",
                           studlab = FALSE, cex.studlab = 0.8,
                           pos.studlab = 2, offset = 0.5,
                           regline = TRUE,
                           backtransf = x$.meta$x$backtransf,
                           ref,
                           col.ref = "lightgray", lty.ref = 1, lwd.ref = 1,
                           pscale = x$.meta$x$pscale,
                           irscale = x$.meta$x$irscale,
                           axes = TRUE, box = TRUE, ...) {
  
  
  #
  #
  # (1) Check for meta object
  #
  #
  chkclass(x, "metareg")
  
  
  #
  #
  # (2) Check other arguments
  #
  #
  chknumeric(min.cex)
  chknumeric(max.cex)
  chknumeric(lty)
  chknumeric(lwd)
  chknumeric(cex.studlab)
  pos.studlab <- as.numeric(setchar(pos.studlab, as.character(1:4)))
  chknumeric(offset)
  chklogical(regline)
  chklogical(backtransf)
  chknumeric(lty.ref)
  chknumeric(lwd.ref)
  #
  if (!is.null(pscale))
    chknumeric(pscale, length = 1)
  else
    pscale <- 1
  #
  if (!is.null(irscale))
    chknumeric(irscale, length = 1)
  else
    irscale <- 1
  #
  if (!is.character(axes)) {
    chklogical(axes)
    axis.x <- axis.y <- axes
  }
  else {
    axes <- setchar(axes, c("x", "y", "xy"))
    #
    axis.x <- grepl("x", axes)
    axis.y <- grepl("y", axes)
  }
  #
  chklogical(box)
  
  
  m0 <- update(x$.meta$x)
  method.tau0 <- x$.meta$method.tau
  #
  if (method.tau0 != "FE" & (method.tau0 != m0$method.tau))
    m1 <- update(m0, method.tau = method.tau0)
  else
    m1 <- m0
  #
  TE <- m1$TE
  sm <- m1$sm

  log.y <- backtransf &&
    (is_relative_effect(sm) || is_log_effect(sm) ||
     (!is.null(m1$func.backtransf) && m1$func.backtransf == "exp"))
  #
  if (log.y) {
    log <- "y"
    TE <- exp(TE)
  }
  else
    log <- ""
  #
  if (missing(ref)) {
    if (is_prop(sm) | is_rate(sm) | is_mean(sm))
      ref <- NA
    else if (log.y)
      ref <- 1
    else
      ref <- 0
  }
  #
  chknumeric(ref, length = 1)
  
  
  if (is.logical(studlab)) {
    if (studlab)
      studlab <- m1$studlab
    else
      studlab <- rep("", length(TE))
  }
  else {
    studlab <- as.character(studlab)
    if (length(studlab) != length(TE))
      stop("Length of argument 'studlab' must be the same as ",
           "number of studies in meta-analysis.")
  }
  
  
  charform <- as.character(x$.meta$formula)[2]
  splitform <- strsplit(charform, " ")[[1]]
  covar.name <- splitform[1]
  if (covar.name == "1" | covar.name == "-1")
    covar.name <- splitform[3]
  #
  covar.names <- names(coef(x))
  covar.names.without.intrcpt <- covar.names[covar.names != "intrcpt"]
  if (length(covar.names.without.intrcpt) == 0) {
    warning("No covariate in meta-regression.")
    return(invisible(NULL))
  }
  #
  nointrcpt <- ifelse("intrcpt" %in% covar.names, FALSE, TRUE)
  #
  if (covar.name == ".subgroup")
    covar.name <- x$.meta$x$subgroup.name
  #
  if (covar.name %in% names(x$.meta$x$data))
    covar <- x$.meta$x$data[[covar.name]]
  else if (covar.name %in% names(x$.meta$x))
    covar <- x$.meta$x[[covar.name]]
  else if (".subgroup" %in% names(x$.meta$x$data))
    covar <- x$.meta$x$data[[".subgroup"]]
  else
    covar <- get(covar.name)
  #
  if (!is.null(x$.meta$x$subset))
    covar <- covar[x$.meta$x$subset]
  #
  if (is.character(covar))
    covar <- as.factor(covar)
  #
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
  #
  alpha <- ifelse(nointrcpt, 0, coef(x)["intrcpt"])
  #
  if (covar.name %in% names(coef(x)))
    beta <- coef(x)[covar.name]
  else
    beta <- coef(x)[".subgroup"]
  #
  if (length(covar.names.without.intrcpt) > 1 & !is.factor(covar)) {
    warning(paste0("Only first covariate in meta-regression ",
                   "('", covar.name, "') considered in bubble plot. ",
                   "No regression line plotted.")
            )
    regline <- FALSE
    if (missing(xlab))
      xlab <- paste0("Covariate ", covar.name,
                     " (meta-regression: ", charform, ")")
  }
  else
    if (missing(xlab))
      xlab = paste("Covariate", covar.name)
  
  
  if (backtransf && (pscale != 1 || irscale != 1)) {
    if (pscale != 1 && irscale != 1)
      stop("Provide either arguments 'pscale' or 'irscale'",
           call. = FALSE)
    if (pscale != 1)
      scale <- pscale
    else
      scale <- irscale
  }
  else
    scale <- 1
  #
  ys <- TE
  #
  if (backtransf & !log.y) {
    if (!is.null(m1$func.backtransf))
      func.backtransf <- m1$func.backtransf
    else if (sm == "PLOGIT")
      func.backtransf <- logit2p
    else if (sm == "PAS")
      func.backtransf <- asin2p
    else if (sm == "IRS")
      func.backtransf <- function(x) x^2
    else if (sm == "ZCOR")
      func.backtransf <- z2cor
    else
      func.backtransf <- I
    #
    ys <- do.call(func.backtransf, list(ys))
  }
  #
  if (missing(ylim))
    ylim <- scale * range(ys)
  #
  if (missing(ylab)) {
    ylab <- xlab_meta(sm, backtransf, func.transf = m1$func.transf,
                 func.backtransf = m1$func.backtransf)
    #
    if (ylab == "") {
      if (sm == "PRAW" | (backtransf & sm %in% c("PLN", "PAS", "PLOGIT")))
        ylab <- "Proportion"
      else if (sm == "IR" | (backtransf & sm %in% c("IRLN", "IRS")))
        ylab <- "Incidence Rate"
      else if (sm == "MRAW" | (backtransf & sm %in% c("IRLN", "IRS")))
        ylab <- "Incidence Rate"
      else if (sm == "COR" | (backtransf & sm == "ZCOR"))
        ylab <- "Correlation"
    }
  }
  
  
  missing.cex <- missing(cex)
  #
  if (!missing.cex && is.character(cex)) {
    cex.type <- setchar(cex, c("common", "random", "fixed"),
                        "must be numeric or equal to \"common\" or \"random\"")
    cex.type[cex.type == "fixed"] <- "common"
    if (length(unique(cex.type)) != 1)
      stop("Argument 'cex' must be numeric or equal to ",
           "\"common\" or \"random\".")
    common.cex <- all(cex.type == "common")
    random.cex <- all(cex.type == "random")
  }
  else {
    common.cex <- FALSE
    random.cex <- FALSE
  }
  #
  if (missing.cex)
    if (method.tau0 == "FE")
      cex <- m1$w.common
    else
      cex <- m1$w.random
  else if (common.cex)
    cex <- m1$w.common
  else if (random.cex)
    cex <- m1$w.random
  else if (length(cex) == 1)
    cex <- rep_len(cex, length(TE))
  #
  if (all(is.na(cex))) {
    cex <- rep_len(1, length(TE))
    min.cex <- 1
    max.cex <- 1
  }
  #
  if (length(cex) != length(TE))
    stop("Length of argument 'cex' must be the same as ",
         "number of studies in meta-analysis.")
  #
  if (missing.cex | common.cex) {
    cexs <- max.cex * (cex / max(cex))
    cexs[cexs < min.cex] <- min.cex
  }
  else
    cexs <- cex
  #
  if (length(col) > 1 && length(col) != length(TE))
    stop("Length of argument 'col' must be 1 or the same as ",
         "number of studies in meta-analysis.")
  #
  if (length(bg) > 1 && length(bg) != length(TE))
    stop("Length of argument 'col' must be 1 or the same as ",
         "number of studies in meta-analysis.")
  #
  o <- rev(order(cex))
  xs <- xs[o]
  ys <- ys[o]
  studlab <- studlab[o]
  cexs <- cexs[o]
  #
  if (length(pch) > 1)
    pchs <- pch[o]
  else
    pchs <- rep(pch, length(o))
  #
  if (length(col) > 1)
    cols <- col[o]
  else
    cols <- rep(col, length(o))
  #
  if (length(bg) > 1)
    bgs <- bg[o]
  else
    bgs <- rep(bg, length(o))
  #
  if (length(cex.studlab) > 1)
    cex.studlab <- cex.studlab[o]
  #
  if (length(pos.studlab) > 1)
    pos.studlab <- pos.studlab[o]
  #
  if (length(offset) > 1)
    offset <- offset[o]
  
  
  #
  #
  # Generate bubble plot
  #
  #
  if (is.factor(covar)) {
    for (i in 2:length(levs)) {
      sel <- xs %in% c(0, i - 1)
      xs.i <- xs[sel]
      xs.i[xs.i > 0] <- 1
      ys.i <- ys[sel]
      studlab.i <- studlab[sel]
      cexs.i <- cexs[sel]
      pch.i <- pchs[sel]
      col.i <- cols[sel]
      bg.i <- bgs[sel]
      #
      plot(xs.i, scale * ys.i,
           pch = pch.i, cex = cexs.i,
           xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
           type = "n", axes = FALSE, log = log, ...)
      #
      # Add reference line
      #
      abline(h = scale * ref, col = col.ref, lty = lty.ref, lwd = lwd.ref)
      #
      # Add regression line
      #
      if (regline) {
        x.reg <- seq(0, 1, length.out = 500)
        y.reg <- seq(alpha, alpha + coef(x)[i], length.out = 500)
        #
        if (log.y)
          y.reg <- scale * exp(y.reg)
        else if (backtransf)
          y.reg <- scale * do.call(func.backtransf, list(y.reg))
        #
        lines(x.reg, y.reg, lty = lty, lwd = lwd, col = col.line)
      }
      #
      for (j in seq_along(xs.i))
        points(xs.i[j], scale * ys.i[j], cex = cexs.i[j], pch = pch.i[j],
               col = col.i[j], bg = bg.i[j])
      #
      # x-axis
      #
      if (axis.x)
        axis(1, at = 0:1, labels = levs[c(1, i)], ...)
      #
      # y-axis
      #
      if (axis.y)
        axis(2, ...)
      #
      text(xs.i, scale * ys.i, labels = studlab.i, cex = cex.studlab,
           pos = pos.studlab, offset = offset)
      #
      if (box)
        box()
    }
  }
  else {
    plot(xs, scale * ys,
         pch = pchs, cex = cexs,
         xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
         type = "n", axes = FALSE, log = log, ...)
    #
    # Add reference line
    #
    abline(h = scale * ref, col = col.ref, lty = lty.ref, lwd = lwd.ref)
    #
    # Add regression line
    #
    if (regline) {
      x.reg <- seq(xlim[1], xlim[2], length.out = 500)
      y.reg <- seq(alpha + beta * xlim[1], alpha + beta * xlim[2],
                   length.out = 500)
      #
      if (log.y)
        y.reg <- scale * exp(y.reg)
      else if (backtransf)
        y.reg <- scale * do.call(func.backtransf, list(y.reg))
      #
      lines(x.reg, y.reg, lty = lty, lwd = lwd, col = col.line)
    }
    #
    points(xs, scale * ys, cex = cexs, pch = pchs, col = cols, bg = bgs)
    #
    # x-axis
    #
    if (axis.x)
      axis(1, ...)
    #
    # y-axis
    #
    if (axis.y)
      axis(2, ...)
    #
    text(xs, scale * ys, labels = studlab, cex = cex.studlab,
         pos = pos.studlab, offset = offset)
    #
    if (box)
      box()
  }
  
  
  invisible(NULL)
}





#' @rdname bubble.metareg
#' @export bubble


bubble <- function(x, ...) 
  UseMethod("bubble")
