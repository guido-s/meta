labbe.metabin <- function(x,
                          xlim, ylim,
                          xlab=NULL, ylab=NULL,
                          TE.fixed=x$TE.fixed,
                          TE.random=x$TE.random,
                          comb.fixed=x$comb.fixed,
                          comb.random=x$comb.random,
                          axes=TRUE,
                          pch=21, text=NULL, cex=1,
                          col="black", bg="lightgray",
                          lwd=1, lwd.fixed=lwd, lwd.random=lwd,
                          lty.fixed=2, lty.random=9,
                          sm=x$sm, weight,
                          studlab=FALSE, cex.studlab=0.8,
                          ...){
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(x, "metabin")
  x <- updateversion(x)
  
  
  pc <- x$event.c/x$n.c
  pe <- x$event.e/x$n.e
  
  if(length(pc) != length(pe))
    stop("event rates must be of same length")
  
  sm <- setchar(sm, c("OR", "RD", "RR", "ASD"))
  
  if (missing(weight))
    weight <- ifelse(comb.random & !comb.fixed, "random", "fixed")
  ##
  iweight <- charmatch(tolower(weight),
                       c("same", "fixed", "random"), nomatch = NA)
  ##
  if(is.na(iweight))
    stop("weight should be \"same\", \"fixed\", or \"random\"")
  ##
  weight <- c("same", "fixed", "random")[iweight]
  ##
  if (weight=="fixed")
    cex.i <- 4*cex*sqrt(x$w.fixed)/sqrt(max(x$w.fixed))
  else if (weight=="random")
    cex.i <- 4*cex*sqrt(x$w.random)/sqrt(max(x$w.random))
  else if (weight=="same")
    cex.i <- cex
  ##
  if (min(cex.i) < 0.5)
    cex.i <- cex.i + (0.5-min(cex.i))
  
  
  if (is.logical(studlab) && studlab)
    studlab <- x$studlab
  
  
  if (length(comb.fixed)==0)
    comb.fixed <- TRUE
  if (length(comb.random)==0)
    comb.random <- TRUE
  
  
  if (missing(xlim) & missing(ylim)){
    xlim <- c(0, max(c(pc, pe), na.rm=TRUE))
    ylim <- xlim
  }
  if (missing(xlim))
    xlim <- c(0, max(c(pc, pe), na.rm=TRUE))
  if (missing(ylim))
    ylim <- xlim
  
  
  oldpar <- par(pty="s")
  on.exit(par(oldpar))
  
  
  if (is.null(xlab))
    if (length(x$label.c)>0)
      xlab <- paste("Event rate (", x$label.c, ")", sep="")
    else
      xlab <- "Event rate (Control)"
  ##
  if (is.null(ylab))
    if (length(x$label.e)>0)
      ylab <- paste("Event rate (", x$label.e, ")", sep="")
    else
      ylab <- "Event rate (Experimental)"
  
  
  if (comb.fixed && length(lty.fixed)==1 & length(TE.fixed)>1)
    lty.fixed <- rep(lty.fixed, length(TE.fixed))
  ##
  if (comb.fixed && length(lwd.fixed)==1 & length(TE.fixed)>1)
    lwd.fixed <- rep(lwd.fixed, length(TE.fixed))
  
  if (comb.random && length(lty.random)==1 & length(TE.random)>1)
    lty.random <- rep(lty.random, length(TE.random))
  ##
  if (comb.random && length(lwd.random)==1 & length(TE.random)>1)
    lwd.random <- rep(lwd.random, length(TE.random))
  
  
  ##
  ## plot
  ##
  plot(pc, pe, type="n", 
       xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
       axes=axes, ...)
  
  points(pc, pe, pch=pch, cex=cex.i, col=col, bg=bg)
  
  if (comb.fixed){
    x.line <- seq(min(xlim), max(xlim), len=100)
    ##
    if (sm=="RR" & length(TE.fixed)>0){
      for (i in 1:length(TE.fixed)){
        y.line <- x.line*exp(TE.fixed[i])
        sel <- min(ylim) <= y.line & y.line <= max(ylim)
        if (sum(sel)>1)
          lines(x.line[sel], y.line[sel],
                lty=lty.fixed[i], lwd=lwd.fixed[i])
      }
    }
    ##
    if (sm=="RD" & length(TE.fixed)>0){
      for (i in 1:length(TE.fixed)){
        y.line <- x.line + TE.fixed[i]
        sel <- min(ylim) <= y.line & y.line <= max(ylim)
        if (sum(sel)>1)
          lines(x.line[sel], y.line[sel],
                lty=lty.fixed[i], lwd=lwd.fixed[i])
      }
    }
    ##
    if (sm=="OR" & length(TE.fixed)>0){
      for (i in 1:length(TE.fixed)){
        y.line <- exp(TE.fixed[i])*(x.line/(1-x.line)) /
          (1+exp(TE.fixed[i])*x.line/(1-x.line))
        sel <- min(ylim) <= y.line & y.line <= max(ylim)
        if (sum(sel, na.rm=TRUE)>1)
          lines(x.line[sel], y.line[sel],
                lty=lty.fixed[i], lwd=lwd.fixed[i])
      }
    }
    ##
    if (sm=="ASD" & length(TE.fixed)>0){
      for (i in 1:length(TE.fixed)){
        y.line <- sin(asin(sqrt(x.line)) + TE.fixed[i])^2
        sel <- min(ylim) <= y.line & y.line <= max(ylim)
        if (sum(sel)>1)
          lines(x.line[sel], y.line[sel],
                lty=lty.fixed[i], lwd=lwd.fixed[i])
      }
    }
  }
  
  if (comb.random){
    x.line <- seq(min(xlim), max(xlim), len=100)
    ##
    if (sm=="RR" & length(TE.random)>0){
      for (i in 1:length(TE.random)){
        y.line <- x.line*exp(TE.random[i])
        sel <- min(ylim) <= y.line & y.line <= max(ylim)
        if (sum(sel)>1)
          lines(x.line[sel], y.line[sel],
                lty=lty.random[i], lwd=lwd.random[i])
      }
    }
    ##
    if (sm=="RD" & length(TE.random)>0){
      for (i in 1:length(TE.random)){
        y.line <- x.line + TE.random[i]
        sel <- min(ylim) <= y.line & y.line <= max(ylim)
        if (sum(sel)>1)
          lines(x.line[sel], y.line[sel],
                lty=lty.random[i], lwd=lwd.random[i])
      }
    }
    ##
    if (sm=="OR" & length(TE.random)>0){
      for (i in 1:length(TE.random)){
        y.line <- exp(TE.random[i])*(x.line/(1-x.line)) /
          (1+exp(TE.random[i])*x.line/(1-x.line))
        sel <- min(ylim) <= y.line & y.line <= max(ylim)
        if (sum(sel, na.rm=TRUE)>1)
          lines(x.line[sel], y.line[sel],
                lty=lty.random[i], lwd=lwd.random[i])
      }
    }
    ##
    if (sm=="ASD" & length(TE.random)>0){
      for (i in 1:length(TE.random)){
        y.line <- sin(asin(sqrt(x.line)) + TE.random[i])^2
        sel <- min(ylim) <= y.line & y.line <= max(ylim)
        if (sum(sel)>1)
          lines(x.line[sel], y.line[sel],
                lty=lty.random[i], lwd=lwd.random[i])
      }
    }
  }
  
  if (!is.logical(studlab) && length(studlab)>0)
    text(pc, pe, labels=studlab, pos=2, cex=cex.studlab)
  
  invisible(NULL)
}
