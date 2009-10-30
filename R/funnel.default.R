funnel.default <- function(x, y,
                           xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
                           comb.fixed=FALSE, comb.random=FALSE,
                           axes=TRUE,
                           pch=21, text=NULL, cex=1,
                           lty.fixed=2, lty.random=9,
                           lwd=1, lwd.fixed=lwd, lwd.random=lwd,
                           col="black", bg="darkgray",
                           col.fixed="black", col.random="black",
                           log="", yaxis="se", sm=NULL,
                           contour.levels=NULL, col.contour,
                           ref=ifelse(sm %in% c("RR", "OR", "HR"), 1, 0),
                           level=NULL,
                           studlab=FALSE, cex.studlab=0.8,
                           ...){
  
  TE <- x
  seTE <- y
  TE.fixed  <- metagen(TE, seTE)$TE.fixed
  TE.random <- metagen(TE, seTE)$TE.random
  if (is.null(sm)) sm <- "other"
  if (missing(ref))
    ref <- ifelse(sm %in% c("RR", "OR", "HR"), 1, 0)
  if (is.logical(studlab) && studlab)
    studlab <- seq(along=TE)
  
  if (length(level)==0){
    level <- NULL
  }
  
  
  if(length(TE) != length(seTE))
    stop("length of argument TE and seTE must be equal")
  
  if (!is.null(level) && (level<=0|level>=1))
    stop("no valid level for confidence interval")
  
  if (!is.null(contour.levels) &&
      any(contour.levels<=0|contour.levels>=1))
    stop("contour levels must be between 0 and 1")
  
  iyaxis <- charmatch(yaxis,
                     c("se", "size", "invvar", "invse"),
                     nomatch = NA)
  if(is.na(iyaxis) | iyaxis==0)
    stop("yaxis should be \"se\", \"size\", \"invvar\", or \"invse\"")
  ##
  yaxis <- c("se", "size", "invvar", "invse")[iyaxis]
  
  seTE[is.infinite(seTE)] <- NA
  
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
    TE.xlim <- 1.025*c(min(c(TE, ciTE$lower), na.rm=TRUE),
                       max(c(TE, ciTE$upper), na.rm=TRUE))
  }
  ##
  if (match(sm, c("OR", "RR", "HR"), nomatch=0)>0){
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
    else if (inherits(x, "meta") && (length(x$n.e)>0 & length(x$n.c)>0))
      weight <- floor(x$n.e)+floor(x$n.c)
    else if (inherits(x, "meta") && length(x$n)>0)
      weight <- floor(x$n)
    else if (inherits(x, "meta"))
      stop("no information on sample size available in object '",
           deparse(substitute(x)), "'")
    else
      weight <- y
  
  
  ##
  ## x-axis: labels / xlim
  ##
  if (is.null(xlab)){
    if      (sm=="OR" ) xlab <- "Odds Ratio"
    else if (sm=="RD" ) xlab <- "Risk Difference"
    else if (sm=="RR" ) xlab <- "Relative Risk"
    else if (sm=="SMD") xlab <- "Standardised mean difference"
    else if (sm=="WMD"|sm=="MD") xlab <- "Mean difference"
    else if (sm=="HR" ) xlab <- "Hazard Ratio"
    else if (sm=="AS" ) xlab <- "Arcus Sinus Transformation"
    else if (sm=="proportion" ) xlab <- "Proportion"
    else xlab <- sm
  }
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
  if (yaxis=="se"   & is.null(ylab)) ylab <- "Standard error"
  if (yaxis=="size" & is.null(ylab)) ylab <- "Study size"
  if (yaxis=="invvar" & is.null(ylab)) ylab <- "Inverse of variance"
  if (yaxis=="invse" & is.null(ylab)) ylab <- "Inverse of standard error"
  ##
  if (is.null(ylim) & yaxis=="se") ylim <- c(max(weight, na.rm=TRUE), 0)
  if (is.null(ylim)              ) ylim <- range(weight, na.rm=TRUE)
  
  ##
  ## plot
  ##
  plot(TE, weight, type="n",
       xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
       axes=axes, log=log, ...)
  
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
    if (match(sm, c("OR", "RR", "HR"), nomatch=0)>0)
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
      if (match(sm, c("OR", "RR", "HR"), nomatch=0)>0){
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
  
  if (is.null(text))
    points(TE, weight, pch=pch, cex=cex, col=col, bg=bg)
  else
    text(TE, weight, labels=text, cex=cex, col=col)
  ##  if (is.null(text)){
  ##    if (is.null(col))
  ##      points(TE, weight, pch=pch, cex=cex)
  ##    else
  ##      if (all(pch==21))
  ##        points(TE, weight, pch=pch, cex=cex, bg=col)
  ##      else
  ##        points(TE, weight, pch=pch, cex=cex, col=col)
  ##  }
  ##  else{
  ##    if (is.null(col))
  ##      text(TE, weight, labels=text, cex=cex)
  ##    else
  ##      text(TE, weight, labels=text, cex=cex, col=col)
  ##  }
  
  if (comb.fixed)
    lines(c(TE.fixed, TE.fixed), range(ylim), lty=lty.fixed, lwd=lwd.fixed)
  ##abline(v=TE.fixed, lty=lty.fixed, lwd=lwd.fixed)
  ##abline(v=TE.fixed, lty=2)
  
  if (comb.random)
    lines(c(TE.random, TE.random), range(ylim), lty=lty.random, lwd=lwd.random)
    ##abline(v=TE.random, lty=lty.random, lwd=lwd.random)
  
  if (!is.null(level) & yaxis=="se"){
    points(ciTE$lower, seTE.seq, type="l", lty=lty.fixed, lwd=lwd.fixed)
    points(ciTE$upper, seTE.seq, type="l", lty=lty.fixed, lwd=lwd.fixed)
  }
  
  if (!is.null(level) & yaxis=="invvar"){
    points(ciTE$lower, 1/seTE.seq^2, type="l", lty=lty.fixed, lwd=lwd.fixed)
    points(ciTE$upper, 1/seTE.seq^2, type="l", lty=lty.fixed, lwd=lwd.fixed)
  }
  
  if (!is.null(level) & yaxis=="invse"){
    points(ciTE$lower, 1/seTE.seq, type="l", lty=lty.fixed, lwd=lwd.fixed)
    points(ciTE$upper, 1/seTE.seq, type="l", lty=lty.fixed, lwd=lwd.fixed)
  }

  if (!is.null(contour.levels) & yaxis!="size")
    res <- list(col.contour=col.contour)
  else
    res <- NULL

  if (!is.logical(studlab) && length(studlab)>0)
    text(TE, weight, labels=studlab, pos=2, cex=cex.studlab)  
  
  invisible(res)
}
