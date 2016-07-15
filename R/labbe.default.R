labbe.default <- function(x, y,
                          xlim, ylim,
                          xlab = NULL, ylab = NULL,
                          TE.fixed = NULL, TE.random = NULL,
                          comb.fixed = !is.null(TE.fixed), comb.random = !is.null(TE.random),
                          backtransf = TRUE,
                          axes = TRUE,
                          pch = 21, text = NULL, cex = 1,
                          col = "black", bg = "lightgray",
                          lwd = 1, lwd.fixed = lwd, lwd.random = lwd,
                          lty.fixed = 2, lty.random = 9,
                          col.fixed = col, col.random = col,
                          nulleffect = TRUE,
                          lwd.nulleffect = lwd, col.nulleffect = "lightgray",
                          sm = NULL, weight,
                          studlab = FALSE, cex.studlab = 0.8,
                          label.e = NULL, label.c = NULL,
                          ...) {
  
  xpos <- x
  ypos <- y
  
  if(length(xpos) != length(ypos))
    stop("arguments 'x' and 'y' must be of same length")
  
  sm <- setchar(sm, c("OR", "RD", "RR", "ASD"))
  
  if (!backtransf) {
    if (sm == "OR") {
      xpos <- log(xpos / (1 - xpos))
      ypos <- log(ypos / (1 - ypos))
    }
    else if (sm == "RR") {
      xpos <- log(xpos)
      ypos <- log(ypos)
    }
    else if (sm == "ASD") {
      xpos <- asin(sqrt(xpos))
      ypos <- asin(sqrt(ypos))
    }
  }
  
  if (!missing(weight))
    cex.i <- 4 * cex * sqrt(weight) / sqrt(max(weight))
  else
    cex.i <- rep(cex, length(xpos))
  ##
  if (min(cex.i) < 0.5)
    cex.i <- cex.i + (0.5 - min(cex.i))
  
  
  if (backtransf)
    minval <- 0
  else
    minval <- min(c(xpos, ypos), na.rm = TRUE)
  ##
  if (missing(xlim) & missing(ylim)) {
    xlim <- c(minval, max(c(xpos, ypos), na.rm = TRUE))
    ylim <- xlim
  }
  if (missing(xlim))
    xlim <- c(minval, max(c(xpos, ypos), na.rm = TRUE))
  if (missing(ylim))
    ylim <- xlim
  
  
  oldpar <- par(pty = "s")
  on.exit(par(oldpar))
  
  
  if (is.null(xlab)) {
    if (length(label.c) > 0) {
      xlab <- paste("Event rate (", label.c, ")", sep = "")
      if (!backtransf)
        if (sm == "OR")
          xlab <- paste("ln (odds) ", label.c, sep = "")
        else if (sm == "RR")
          xlab <- paste("ln (event rate) ", label.c, sep = "")
        else if (sm == "ASD")
          xlab <- paste("Arcsin-transformed event rate (", label.c, ")", sep = "")
    }
    else {
      xlab <- "Event rate (Control)"
      if (!backtransf)
        if (sm == "OR")
          xlab <- "ln (odds) Control"
        else if (sm == "RR")
          xlab <- "ln (event rate) Control"
        else if (sm == "ASD")
          xlab <- "Arcsin-transformed event rate (Control)"
    }
  }
  ##
  if (is.null(ylab)) {
    if (length(label.e) > 0) {
      ylab <- paste("Event rate (", label.e, ")", sep = "")
      if (!backtransf)
        if (sm == "OR")
          ylab <- paste("ln (odds) ", label.e, sep = "")
        else if (sm == "RR")
          ylab <- paste("ln (event rate) ", label.e, sep = "")
        else if (sm == "ASD")
          ylab <- paste("Arcsin-transformed event rate (", label.e, ")", sep = "")
    }
    else {
      ylab <- "Event rate (Experimental)"
      if (!backtransf)
        if (sm == "OR")
          ylab <- "ln (odds) Experimental"
        else if (sm == "RR")
          ylab <- "ln (event rate) Experimental"
        else if (sm == "ASD")
          ylab <- "Arcsin-transformed event rate (Experimental)"
    }
  }
  
  
  if (comb.fixed && length(lty.fixed) == 1 & length(TE.fixed) > 1)
    lty.fixed <- rep(lty.fixed, length(TE.fixed))
  ##
  if (comb.fixed && length(lwd.fixed) == 1 & length(TE.fixed) > 1)
    lwd.fixed <- rep(lwd.fixed, length(TE.fixed))
  ##
  if (comb.fixed && length(col.fixed) == 1 & length(TE.fixed) > 1)
    col.fixed <- rep(col.fixed, length(TE.fixed))
  
  if (comb.random && length(lty.random) == 1 & length(TE.random) > 1)
    lty.random <- rep(lty.random, length(TE.random))
  ##
  if (comb.random && length(lwd.random) == 1 & length(TE.random) > 1)
    lwd.random <- rep(lwd.random, length(TE.random))
  ##
  if (comb.random && length(col.random) == 1 & length(TE.random) > 1)
    col.random <- rep(col.random, length(TE.random))
  
  
  ##
  ## Generate L'Abbe plot
  ##
  plot(xpos, ypos, type = "n", 
       xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       axes = axes, ...)
  ##
  if (nulleffect)
    abline(0, 1, lwd = lwd.nulleffect, col = col.nulleffect)
  ##
  points(xpos, ypos, pch = pch, cex = cex.i, col = col, bg = bg, lwd = lwd)
  
  
  ##
  ## Auxillary function
  ##
  addlines <- function(x.line, y.line, ylim, lty, lwd, col) {
    sel <- min(ylim) <= y.line & y.line <= max(ylim)
    if (sum(sel) > 1)
      lines(x.line[sel], y.line[sel],
            lty = lty, lwd = lwd, col = col)
  }
  
  
  ##
  ## Add results for common effect model
  ##
  if (comb.fixed & length(TE.fixed) > 0) {
    x.line <- seq(min(xlim), max(xlim), len = 100)
    ##
    if (!backtransf)
      for (i in 1:length(TE.fixed))
        abline(TE.fixed[i], 1,
               lty = lty.fixed[i], lwd = lwd.fixed[i],
               col = col.fixed[i])
    else {
      if (sm == "RR") {
        for (i in 1:length(TE.fixed)) {
          y.line <- x.line * exp(TE.fixed[i])
          addlines(x.line, y.line, ylim,
                   lty.fixed[i], lwd.fixed[i], col.fixed[i])
        }
      }
      else if (sm == "RD") {
        for (i in 1:length(TE.fixed)) {
          y.line <- x.line + TE.fixed[i]
          addlines(x.line, y.line, ylim,
                   lty.fixed[i], lwd.fixed[i], col.fixed[i])
        }
      }
      else if (sm == "OR") {
        for (i in 1:length(TE.fixed)) {
          y.line <- exp(TE.fixed[i]) * (x.line / (1 - x.line)) /
            (1 + exp(TE.fixed[i]) * x.line / (1 - x.line))
          addlines(x.line, y.line, ylim,
                   lty.fixed[i], lwd.fixed[i], col.fixed[i])
        }
      }
      else if (sm == "ASD" & length(TE.fixed) > 0) {
        for (i in 1:length(TE.fixed)) {
          y.line <- sin(asin(sqrt(x.line)) + TE.fixed[i])^2
          addlines(x.line, y.line, ylim,
                   lty.fixed[i], lwd.fixed[i], col.fixed[i])
        }
      }
    }
  }
  
  
  ##
  ## Add results for random effects model
  ##
  if (comb.random & length(TE.random) > 0) {
    x.line <- seq(min(xlim), max(xlim), len = 100)
    ##
    if (!backtransf)
      for (i in 1:length(TE.random))
        abline(TE.random[i], 1,
               lty = lty.random[i], lwd = lwd.random[i],
               col = col.random[i])
    else {
      if (sm == "RR") {
        for (i in 1:length(TE.random)) {
          y.line <- x.line * exp(TE.random[i])
          addlines(x.line, y.line, ylim,
                   lty.random[i], lwd.random[i], col.random[i])
        }
      }
      else if (sm == "RD") {
        for (i in 1:length(TE.random)) {
          y.line <- x.line + TE.random[i]
          addlines(x.line, y.line, ylim,
                   lty.random[i], lwd.random[i], col.random[i])
        }
      }
      else if (sm == "OR") {
        for (i in 1:length(TE.random)) {
          y.line <- exp(TE.random[i]) * (x.line / (1 - x.line)) /
            (1 + exp(TE.random[i]) * x.line / (1 - x.line))
          addlines(x.line, y.line, ylim,
                   lty.random[i], lwd.random[i], col.random[i])
        }
      }
      else if (sm == "ASD" & length(TE.random) > 0) {
        for (i in 1:length(TE.random)) {
          y.line <- sin(asin(sqrt(x.line)) + TE.random[i])^2
          addlines(x.line, y.line, ylim,
                   lty.random[i], lwd.random[i], col.random[i])
        }
      }
    }
  }


  ##
  ## Add study labels
  ##
  if (!is.logical(studlab) && length(studlab) > 0)
    text(xpos, ypos, labels = studlab, pos = 2, cex = cex.studlab)
  
  
  invisible(NULL)
}
