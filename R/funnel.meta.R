funnel.meta <- function(x,
                        ##
                        xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
                        ##
                        comb.fixed=x$comb.fixed, comb.random=x$comb.random,
                        ##
                        axes=TRUE,
                        pch=if (!inherits(x, "trimfill")) 21 else ifelse(x$trimfill, 1, 21),
                        text=NULL, cex=1,
                        lty.fixed=2, lty.random=9,
                        lwd=1, lwd.fixed=lwd, lwd.random=lwd,
                        col="black", bg="darkgray",
                        col.fixed="black", col.random="black",
                        ##
                        log="", yaxis="se",
                        contour.levels=NULL, col.contour,
                        ##
                        ref=ifelse(backtransf & is.relative.effect(x$sm), 1, 0),
                        ##
                        level=x$level,
                        studlab=FALSE, cex.studlab=0.8,
                        ##
                        backtransf=x$backtransf,
                        ...){
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta object
  ##
  ##
  chkclass(x, "meta")
  ##
  if (inherits(x, "metacum"))
    stop("Funnel plot not meaningful for object of class \"metacum\"")
  ##
  if (inherits(x, "metainf"))
    stop("Funnel plot not meaningful for object of class \"metainf\"")
  ##
  x <- updateversion(x)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(axes)
  chknumeric(cex)
  chknumeric(lty.fixed)
  chknumeric(lty.random)
  chknumeric(lwd)
  chknumeric(lwd.fixed)
  chknumeric(lwd.random)
  yaxis <- setchar(yaxis, c("se", "size", "invvar", "invse"))
  if (!is.null(contour.levels))
    chklevel(contour.levels, single=FALSE, ci=FALSE)
  chknumeric(ref)
  if (!is.null(level))
    chklevel(level)
  chknumeric(cex.studlab)
  chklogical(backtransf)
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  TE <- x$TE
  k.All <- length(TE)
  seTE <- x$seTE
  seTE[is.infinite(seTE)] <- NA
  ##
  if (is.logical(studlab))
    if (studlab){
      slab <- TRUE
      studlab <- x$studlab
    }
    else
      slab <- FALSE
  else
    slab <- TRUE
  ##
  fun <- "funnel"
  chklength(seTE, k.All, fun)
  if (slab)
    chklength(studlab, k.All, fun)
  if (!is.null(text))
    chklength(text, k.All, fun)
  
  
  ##
  ##
  ## (4) Further assignments
  ##
  ##
  TE.fixed <- x$TE.fixed
  TE.random <- x$TE.random
  sm <- x$sm
  ##  
  if ( yaxis == "se" )
    seTE.min <- 0
  else
    seTE.min <- min(seTE, na.rm=TRUE)
  ##
  seTE.max <- max(seTE, na.rm=TRUE)
  ##
  if (!is.null(level)){
    ##
    seTE.seq <- seq(seTE.min, seTE.max, length.out=500)
    ##
    ciTE <- ci(TE.fixed, seTE.seq, level)
    ##
    TE.xlim <- c(min(c(TE, ciTE$lower), na.rm=TRUE)/1.025,
                 1.025*max(c(TE, ciTE$upper), na.rm=TRUE))
  }
  ##
  if (backtransf & is.relative.effect(sm)){
    TE <- exp(TE)
    TE.fixed <- exp(TE.fixed)
    TE.random <- exp(TE.random)
    if (!is.null(level)){
      ciTE$lower <- exp(ciTE$lower)
      ciTE$upper <- exp(ciTE$upper)
      TE.xlim <- exp(TE.xlim)
    }
    ##
    if (log=="") log <- "x"
  }
  ##
  ## y-value: weight
  ##
  if (yaxis=="invvar") weight <- 1/seTE^2
  if (yaxis=="invse")  weight <- 1/seTE
  if (yaxis=="se") weight <- seTE
  if (yaxis=="size")
    if (inherits(x, "metabin") || inherits(x, "metacont"))
      weight <- floor(x$n.e)+floor(x$n.c)
    else if (length(x$n.e)>0 & length(x$n.c)>0)
      weight <- floor(x$n.e)+floor(x$n.c)
    else if (inherits(x, "metaprop"))
      weight <- floor(x$n)
    else if (length(x$n)>0)
      weight <- floor(x$n)
    else
      stop("No information on sample size available in object '",
           deparse(substitute(x)), "'.")
  ##
  ## x-axis: labels / xlim
  ##
  if (is.null(xlab))
    if (is.relative.effect(sm))
      xlab <- xlab(sm, backtransf)
    else if (sm == "PRAW")
      xlab <- "Proportion"
    else
      xlab <- xlab(sm, FALSE)
  ##
  if (is.null(xlim) & !is.null(level) &
      (yaxis == "se" |
       yaxis == "invse" |
       yaxis == "invvar"))
    xlim <- TE.xlim
  else if (is.null(xlim))
    xlim <- range(TE, na.rm=TRUE)
  ##
  ## y-axis: labels / ylim
  ##
  if (yaxis=="se"     & is.null(ylab)) ylab <- "Standard error"
  if (yaxis=="size"   & is.null(ylab)) ylab <- "Study size"
  if (yaxis=="invvar" & is.null(ylab)) ylab <- "Inverse of variance"
  if (yaxis=="invse"  & is.null(ylab)) ylab <- "Inverse of standard error"
  ##
  if (is.null(ylim) & yaxis=="se") ylim <- c(max(weight, na.rm=TRUE), 0)
  if (is.null(ylim)              ) ylim <- range(weight, na.rm=TRUE)
  

  ##
  ##
  ## (5) Produce funnel plot
  ##
  ##
  plot(TE, weight, type="n",
       xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
       axes=axes, log=log, ...)
  ##
  ## Add contour shades (enhanced funnel plots)
  ##
  if (!is.null(contour.levels) & yaxis!="size"){
    ##
    if (missing(col.contour))
      if (length(contour.levels)<2)
        col.contour <- "gray50"
      else
        col.contour <- gray(seq(0.5, 0.9, len=length(contour.levels)))
    ##
    if (length(contour.levels)!=length(col.contour))
      stop("arguments 'contour.levels' and 'col.contour' must be of the same length")
    ##
    seTE.cont <- seq(seTE.max, seTE.min, length.out=500)
    ##
    if (is.relative.effect(sm))
      ref <- log(ref)
    ##
    j <- 0
    ##
    for (i in contour.levels){
      ##
      j <- j+1
      ##
      ciContour <- ci(ref, seTE.cont, i)
      ##
      if (is.relative.effect(sm)){
        ciContour$TE    <- exp(ciContour$TE)
        ciContour$lower <- exp(ciContour$lower)
        ciContour$upper <- exp(ciContour$upper)
      }
      sel.l <- ciContour$lower>min(xlim) & ciContour$lower<max(xlim)
      sel.u <- ciContour$upper>min(xlim) & ciContour$upper<max(xlim)
      ##
      min.l <- min(ciContour$lower)
      max.l <- max(ciContour$lower)
      min.u <- min(ciContour$upper)
      max.u <- max(ciContour$upper)
      ##
      if (yaxis=="se"){
        if (max.u<min(xlim) | min.l > max(xlim)){
          contour.u.x <- c(min(xlim), min(xlim), max(xlim), max(xlim))
          contour.u.y <- c(min(ylim), max(ylim), max(ylim), min(ylim))
          contour.l.x <- NA
          contour.l.y <- NA
        }
        else {
          contour.l.x <- c(min(xlim), min(xlim),
                           ciContour$lower[sel.l],
                           if (any(sel.l)) max(ciContour$lower[sel.l]) else NA)
          contour.l.y <- c(min(ylim), max(ylim), seTE.cont[sel.l],
                           min(ylim))
          ##
          contour.u.x <- c(max(xlim), max(xlim),
                           ciContour$upper[sel.u],
                           if (any(sel.u)) min(ciContour$upper[sel.u]) else NA)
          contour.u.y <- c(min(ylim), max(ylim), seTE.cont[sel.u],
                           min(ylim))
        }
      }
      ##
      if (yaxis=="invvar"){
        if (max.u<min(xlim) | min.l > max(xlim)){
          contour.u.x <- c(min(xlim), min(xlim), max(xlim), max(xlim))
          contour.u.y <- c(max(ylim), min(ylim), min(ylim), max(ylim))
          contour.l.x <- NA
          contour.l.y <- NA
        }
        else {
          contour.l.x <- c(min(xlim), min(xlim),
                           ciContour$lower[sel.l],
                           if (any(sel.l)) max(ciContour$lower[sel.l]) else NA)
          contour.l.y <- c(max(ylim), min(ylim), 1/seTE.cont[sel.l]^2,
                           max(ylim))
          ##
          contour.u.x <- c(max(xlim), max(xlim),
                           ciContour$upper[sel.u],
                           if (any(sel.u)) min(ciContour$upper[sel.u]) else NA)
          contour.u.y <- c(max(ylim), min(ylim), 1/seTE.cont[sel.u]^2,
                           max(ylim))
        }
      }
      ##
      if (yaxis=="invse"){
        if (max.u<min(xlim) | min.l > max(xlim)){
          contour.u.x <- c(min(xlim), min(xlim), max(xlim), max(xlim))
          contour.u.y <- c(max(ylim), min(ylim), min(ylim), max(ylim))
          contour.l.x <- NA
          contour.l.y <- NA
        }
        else {
          contour.l.x <- c(min(xlim), min(xlim),
                           ciContour$lower[sel.l],
                           if (any(sel.l)) max(ciContour$lower[sel.l]) else NA)
          contour.l.y <- c(max(ylim), min(ylim), 1/seTE.cont[sel.l],
                           max(ylim))
          ##
          contour.u.x <- c(max(xlim), max(xlim),
                           ciContour$upper[sel.u],
                           if (any(sel.u)) min(ciContour$upper[sel.u]) else NA)
          contour.u.y <- c(max(ylim), min(ylim), 1/seTE.cont[sel.u],
                           max(ylim))
        }
      }
      ##
      polygon(contour.l.x, contour.l.y,
              col=col.contour[j], border=FALSE)
      ##
      polygon(contour.u.x, contour.u.y,
              col=col.contour[j], border=FALSE)
    }
  }
  ##
  ## Add results for individual studies
  ##
  if (is.null(text))
    points(TE, weight, pch=pch, cex=cex, col=col, bg=bg)
  else
    text(TE, weight, labels=text, cex=cex, col=col)
  ##
  ## Add results for meta-analysis
  ##
  if (comb.fixed)
    lines(c(TE.fixed, TE.fixed), range(ylim), lty=lty.fixed, lwd=lwd.fixed, col=col.fixed)
  ##
  if (comb.random)
    lines(c(TE.random, TE.random), range(ylim), lty=lty.random, lwd=lwd.random, col=col.random)
  ##
  ## Add approximate confidence intervals
  ##
  if (!is.null(level))
    if (yaxis=="se"){
      points(ciTE$lower, seTE.seq, type="l", lty=lty.fixed, lwd=lwd.fixed)
      points(ciTE$upper, seTE.seq, type="l", lty=lty.fixed, lwd=lwd.fixed)
    }
    else if (yaxis=="invvar"){
      points(ciTE$lower, 1/seTE.seq^2, type="l", lty=lty.fixed, lwd=lwd.fixed)
      points(ciTE$upper, 1/seTE.seq^2, type="l", lty=lty.fixed, lwd=lwd.fixed)
    }
    else if (yaxis=="invse"){
      points(ciTE$lower, 1/seTE.seq, type="l", lty=lty.fixed, lwd=lwd.fixed)
      points(ciTE$upper, 1/seTE.seq, type="l", lty=lty.fixed, lwd=lwd.fixed)
    }
  ##
  ## Add study labels
  ##
  if (!is.logical(studlab) && length(studlab)>0)
    text(TE, weight, labels=studlab, pos=2, cex=cex.studlab)  
  
  
  ##
  ##
  ## (6) Return contour levels (if not NULL)
  ##
  ##
  if (!is.null(contour.levels) & yaxis!="size")
    res <- list(col.contour=col.contour)
  else
    res <- NULL
  
  
  invisible(res)
}
