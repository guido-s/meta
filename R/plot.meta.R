plot.meta <- function(x,
                      byvar=x$byvar,
                      bylab=x$bylab,
                      print.byvar=x$print.byvar,
                      sortvar,
                      studlab=TRUE,
                      level=x$level,
                      level.comb=x$level.comb,
                      ##
                      comb.fixed=x$comb.fixed, comb.random=x$comb.random,
                      overall=TRUE,
                      text.fixed="Fixed effect model",
                      text.random="Random effects model",
                      lty.fixed=2, lty.random=3,
                      ##
                      xlab=NULL,
                      xlim, ylim, lwd=1, cex=1,
                      cex.comb=1.2*cex,
                      cex.axis=cex,
                      cex.lab=cex,
		      log=ifelse(x$sm %in% c("RR", "OR", "HR"), "x", ""),
                      axes=TRUE,
                      allstudies=TRUE,
                      weight=ifelse(comb.random, "random", "fixed"),
                      scale.diamond=1,
                      scale.square=1,
                      col.i="black",
                      clim=xlim,
                      arrow.length=0.1,
		      ref=ifelse(x$sm %in% c("RR", "OR", "HR"), 1, 0),
                      ...){
  
  
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")
  
  
  byvar.name <- deparse(substitute(byvar))
  ##
  if (!missing(byvar) & length(x$byvar)==0){
    if (!is.null(x[[byvar.name]]))
      byvar <- x[[byvar.name]]
  }
  
  
  if (length(comb.fixed)==0){
    comb.fixed <- FALSE
  }
  ##
  if (length(comb.random)==0){
    comb.random <- FALSE
  }
  ##
  if (length(print.byvar)==0){
    print.byvar <- TRUE
  }

  
  oldpar <- par(las=1, mar=c(5.1, 2.1, 4.1, 2.1))
  on.exit(par(oldpar))
  
  
  if (inherits(x, "metainf")|inherits(x, "metacum")){
    ##
    x$TE.fixed    <- rev(x$TE)[1]
    x$seTE.fixed  <- rev(x$seTE)[1]
    x$TE.random   <- rev(x$TE)[1]
    x$seTE.random <- rev(x$seTE)[1]
    ##
    x$TE <- rev(rev(x$TE)[-(1:2)])
    x$seTE <- rev(rev(x$seTE)[-(1:2)])
    x$studlab <- rev(rev(x$studlab)[-(1:2)])
    ##
    if (x$pooled=="fixed")
      text.fixed <- "Fixed effect model"
    else
      text.fixed <- "Random effects model"
    ##
    comb.fixed <- TRUE
  }
  
  sm <- x$sm

  
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
  
  
  iweight <- charmatch(tolower(weight),
                       c("same", "fixed", "random"), nomatch = NA)
  ##
  if(is.na(iweight))
    stop("weight should be \"same\", \"fixed\", or \"random\"")
  ##
  weight <- c("same", "fixed", "random")[iweight]
  
  
  ##
  ## total number of trials to plot (*not* number of trials combined)
  ##
  k.all <- length(x$TE)


  if (allstudies) n.stud <- k.all
  else n.stud <- x$k # number of trials combined in meta-analysis

  
  by <- length(byvar)>0
  sort <- !missing(sortvar)
  
  
  if (!by) byvar <- rep(1, k.all)
  if (!sort) sortvar <- rep(1, k.all)

  if (sort & length(sortvar) != k.all)
    stop("'x' and 'sortvar' have different length")
  if (by & length(byvar) != k.all)
    stop("'x' and 'byvar' have different length")
  if (by & any(is.na(byvar)))
    stop("Missing values in 'byvar'")


  if (length(studlab) == 1 & is.logical(studlab))
    if (studlab == FALSE) studlab <- rep("", k.all)
    else studlab <- x$studlab
  ##
  if (length(studlab) != k.all)
    stop("'x' and 'studlab' have different length")
  
  if (allstudies)
    sel <- 1:length(x$TE)
  else
    sel <- !is.na(x$TE)
  ##
  if (n.stud != sum(sel>0)) warning("n.stud != sum(sel)")
  ##
  ##
  ##
  x$n <- x$n[sel]
  x$n.e <- x$n.e[sel]
  x$n.c <- x$n.c[sel]
  ##
  x$event <- x$event[sel]
  x$event.e <- x$event.e[sel]
  x$event.c <- x$event.c[sel]
  ##
  x$mean.e <- x$mean.e[sel]
  x$mean.c <- x$mean.c[sel]
  ##
  x$sd.e <- x$sd.e[sel]
  x$sd.c <- x$sd.c[sel]
  ##
  x$TE <- x$TE[sel]
  x$seTE <- x$seTE[sel]
  ##
  x$w.fixed <- x$w.fixed[sel]
  x$w.random <- x$w.random[sel]
  x$studlab <- studlab[sel]
  ##
  byvar <- byvar[sel]
  sortvar <- sortvar[sel]

  
  if (sort | by){
    ##
    o <- order(byvar, sortvar)
    ##
    x$n <- x$n[o]
    x$n.e <- x$n.e[o]
    x$n.c <- x$n.c[o]
    ##
    x$event <- x$event[o]
    x$event.e <- x$event.e[o]
    x$event.c <- x$event.c[o]
    ##
    x$mean.e <- x$mean.e[o]
    x$mean.c <- x$mean.c[o]
    ##
    x$sd.e <- x$sd.e[o]
    x$sd.c <- x$sd.c[o]
    ##
    x$TE <- x$TE[o]
    x$seTE <- x$seTE[o]
    ##
    x$w.fixed <- x$w.fixed[o]
    x$w.random <- x$w.random[o]
    x$studlab <- x$studlab[o]
    ##
    byvar <- byvar[o]
    sortvar <- sortvar[o]
  }

  
  by.levs <- unique(byvar)
  
  
  if (length(level)==0){
    warning("level set to 0.95")
    level <- 0.95
  }
  ##
  if (length(level.comb)==0){
    if (comb.fixed | comb.random)
      warning("level.comb set to 0.95")
    level.comb <- 0.95
  }
  
  sum.meta <- summary(x, level=level, level.comb=level.comb, warn=FALSE)
  ##
  TE <- sum.meta$study$TE
  lowTE <- sum.meta$study$lower
  uppTE <- sum.meta$study$upper
  ##
  TE.fixed <- sum.meta$fixed$TE
  lowTE.fixed <- sum.meta$fixed$lower
  uppTE.fixed <- sum.meta$fixed$upper
  ##
  TE.random <- sum.meta$random$TE
  lowTE.random <- sum.meta$random$lower
  uppTE.random <- sum.meta$random$upper
  
  
  if (by){

    res.w <- matrix(NA, ncol=7, nrow=length(by.levs))
    j <- 0
    ##
    for ( i in by.levs){
      j <- j+1
      sel <- byvar == i
      ##
      if (inherits(x, "metabin")){
        meta1 <- metabin(x$event.e[sel], x$n.e[sel],
                         x$event.c[sel], x$n.c[sel],
                         method=x$method, sm=x$sm,
                         incr=x$incr, allincr=x$allincr,
                         addincr=x$addincr, allstudies=x$allstudies,
                         MH.exact=x$MH.exact, RR.cochrane=x$RR.cochrane,
                         warn=x$warn)
      }
      ##
      if (inherits(x, "metacont")){
        meta1 <- metacont(x$n.e[sel], x$mean.e[sel], x$sd.e[sel],
                          x$n.c[sel], x$mean.c[sel], x$sd.c[sel],
                          sm=x$sm)
      }
      ##
      if (inherits(x, "metagen")){
        meta1 <- metagen(x$TE[sel], x$seTE[sel], sm=x$sm)
      }
      ##
      if (inherits(x, "metaprop")){
        meta1 <- metaprop(event=x$event[sel], n=x$n[sel],
                          freeman.tukey=x$freeman.tukey)
      }
      ##
      res.w[j,] <- c(meta1$TE.fixed, meta1$seTE.fixed,
                     meta1$TE.random, meta1$seTE.random,
                     meta1$Q, meta1$k, length(meta1$TE))
      ##
    }
    ##
    k.w     <- res.w[,6]
    k.all.w <- res.w[,7]
    ##
    if (allstudies) sel <- 1:length(k.w)
    else sel <- k.w>0
    ##
    TE.fixed.w    <- res.w[sel,1]
    seTE.fixed.w  <- res.w[sel,2]
    TE.random.w   <- res.w[sel,3]
    seTE.random.w <- res.w[sel,4]
    Q.w           <- res.w[sel,5]
    ##
    by.levs <- by.levs[sel]
    k.all.w <- k.all.w[sel]
    ##
    if (comb.fixed | !(comb.fixed|comb.random)){
      if (comb.random)
        warning("Estimate from fixed effect model used in groups defined by 'byvar'")
      ##
      meta.w <- ci(TE.fixed.w, seTE.fixed.w, level)
      text.by <- text.fixed
    }
    if (!comb.fixed & comb.random){
      meta.w <- ci(TE.random.w, seTE.random.w, level)
      text.by <- text.random
    }
    ##
    TE.w <- meta.w$TE
    lowTE.w <- meta.w$lower
    uppTE.w <- meta.w$upper
  }
  
  
  if (by) n.by <- length(by.levs)
  else n.by <- 0

  
  if (comb.fixed & comb.random){
    dy.comb <- 3.0
    yTE.fixed <- -0.5
    yTE.random <- -2.0
  }
  if (comb.fixed & !comb.random){
    dy.comb <- 1.5
    yTE.fixed <- -0.5
  }
  if (!comb.fixed & comb.random){
    dy.comb <- 1.5
    yTE.random <- -0.5
  }
  if (!comb.fixed & !comb.random){
    dy.comb <- 0
  }
  if (!overall){
    dy.comb <- 0
  }
  
  
  if (!comb.fixed) text.fixed <- ""
  if (!comb.random) text.random <- ""
  
  
  ##
  ## Plot symbol
  ##
  if (weight=="same"|inherits(x, "metainf")|inherits(x, "metacum"))
    cex.i <- rep(cex, n.stud)
  else
    if (weight=="fixed")
      cex.i <- 3*x$w.fixed/max(x$w.fixed)*scale.square
    else
      cex.i <- 3*x$w.random/max(x$w.random)*scale.square


  ##
  ## y-axis:
  ##
  rows.w <- 3
  N <- n.stud + rows.w*n.by
  ##
  if (!by){
    yTE <- N:1
    yTE.w <- NA
  }
  else{
    ##
    ## Don't ask me why I did it this way, but it works.
    ##
    sel <- rep(rep(c(FALSE,TRUE), n.by),
               as.vector(rbind(rep(rows.w,n.by),
                               rev(k.all.w))))
    ##
    yTE <- rev((1:N)[sel])
    yTE.w <- rev((1:N)[!sel][2+cumsum(c(0, rep(rows.w, n.by-1)))])
    ##
    ##
    if ( (!comb.fixed & !comb.random) | !overall ){
      yTE <- yTE - 1
      yTE.w <- yTE.w - 1
      N <- N-1
    }
  }
  ##
  if (missing(ylim))
    ylim <- c(0.5-dy.comb, N+0.5)

  dy <- 0.25*scale.diamond

  
  if (sm %in% c("RR", "OR", "HR")){
    TE <- exp(TE)
    lowTE <- exp(lowTE)
    uppTE <- exp(uppTE)
    ##
    TE.fixed <- exp(TE.fixed)
    lowTE.fixed <- exp(lowTE.fixed)
    uppTE.fixed <- exp(uppTE.fixed)
    ##
    TE.random <- exp(TE.random)
    lowTE.random <- exp(lowTE.random)
    uppTE.random <- exp(uppTE.random)
    ##
    if (by){
      TE.w <- exp(TE.w)
      lowTE.w <- exp(lowTE.w)
      uppTE.w <- exp(uppTE.w)
    }
  }
  
  
  if (inherits(x, "metaprop")){
    ref <- NA
    ##
    denum <- 1 + x$freeman.tukey
    ##
    TE.fixed <- sin(TE.fixed/denum)^2
    lowTE.fixed <- sin(lowTE.fixed/denum)^2
    uppTE.fixed <- sin(uppTE.fixed/denum)^2
    ##
    TE.random <- sin(TE.random/denum)^2
    lowTE.random <- sin(lowTE.random/denum)^2
    uppTE.random <- sin(uppTE.random/denum)^2
    ##
    if (by){
      TE.w <- sin(TE.w/denum)^2
      lowTE.w <- sin(lowTE.w/denum)^2
      uppTE.w <- sin(uppTE.w/denum)^2
    }
  }
  
  ##
  ## x-axis:
  ##
  if (missing(xlim))
    xlim <- c(min(lowTE, na.rm=TRUE), max(uppTE, na.rm=TRUE))

  
  ##
  ## Change code to create bylab
  ## (sc, 4.6.2008)
  ##
  ##if (missing(bylab))
  ##  bylab <- paste("byvar", by.levs, sep=" = ")
  ##else
  ##  bylab <- paste(bylab, by.levs, sep=" = ")
  ##
  if (print.byvar){
    if (missing(bylab))
      bylab <- paste(byvar.name,
                     " = ",
                     format(by.levs), sep="")
    else
      bylab <- paste(bylab,
                     " = ",
                     format(by.levs), sep="")
  }
  else
    bylab <- format(by.levs)

  
  ##
  ## Forest plot
  ##
  plot(NA, NA,
       xlim=xlim, ylim=ylim, log=log,
       xlab="", ylab="", type="n", axes=FALSE,
       ...)
  
  ##
  ## Add axis
  ##
  if (axes){
    axis(1, cex.axis=cex.axis)
    mtext(xlab, 1, line=2, cex=cex.lab)
  }

  ##
  ## Add no effect line
  ##
  ##abline(v=ref, lty=2, lwd=2, col="dimgray")
  abline(v=ref, lty=1, lwd=lwd)
  
  ##
  ## Add line at pooled effect
  ##
  if (comb.fixed)
    segments(TE.fixed, yTE.fixed+dy, TE.fixed, N,
             lty=lty.fixed, lwd=lwd, col="dimgray")
  if (comb.random)
    segments(TE.random, yTE.random+dy, TE.random, N,
             lty=lty.random, lwd=lwd, col="dimgray")
  ##  if (inherits(x, "metainf")|inherits(x, "metacum"))
  ##    abline(v=TE.fixed, lty=3, lwd=2, col="dimgray")
  
  ##
  ## Add single trial results
  ##
  if (length(col.i) == 1)
    col.i <- rep(col.i, length(TE))
  ##
  for ( i in seq(along=TE)){
    if (cex.i[i] < 0.3)
      points(TE[i], yTE[i], pch=3, col=col.i[i], lwd=lwd)
    else
      points(TE[i], yTE[i], cex=as.numeric(cex.i[i]), pch=15, col=col.i[i])
  }
  
  ##
  ## Plot confidence intervals
  ##
  ## Cut confidence limits at value of clim
  ## and add arrows at that position
  ## (add by sc, 4.6.2008):
  ##
  ##segments(lowTE, yTE, uppTE, yTE, lwd=lwd, col=col.i)
  ##
  for (i in 1:length(lowTE)){
    if (!is.na(lowTE[i])){
      if (lowTE[i] < clim[1])
        arrows(TE[i], yTE[i], clim[1], yTE[i],
               lwd=lwd, col=col.i, length=arrow.length)
      else
        segments(lowTE[i], yTE[i], TE[i], yTE[i],
                 lwd=lwd, col=col.i)
      ##
      if (uppTE[i] > clim[2])
        arrows(TE[i], yTE[i], clim[2], yTE[i],
               lwd=lwd, col=col.i, length=arrow.length)
      else
        segments(TE[i], yTE[i], uppTE[i], yTE[i],
                 lwd=lwd, col=col.i)
    }
  }
  
  ##
  ## Add trial labels
  ##
  mtext(x$studlab, side=2, line=-0.2, outer=FALSE,
        at=yTE+0.25, adj=0, cex=0.7*cex)
  
  ##
  ## Add group results
  ##
  if (by){
    ##
    if (comb.fixed|comb.random){
      for ( i in 1:length(lowTE.w)){
        polygon(c(lowTE.w[i], TE.w[i], uppTE.w[i], TE.w[i]),
                c(yTE.w[i], yTE.w[i]-dy, yTE.w[i], yTE.w[i]+dy),
                col="lightgray")
      }
    }
    ##
    mtext(bylab,
          side=2, line=-0.2, outer=FALSE,
          at=yTE.w+0.25,
          adj=0, cex=0.7*cex)
  }

  ##
  ## Add pooled estimates
  ##
  if (comb.fixed & overall){
    ##
    ## FE estimate
    ##
    polygon(c(lowTE.fixed, TE.fixed, uppTE.fixed, TE.fixed),
            c(yTE.fixed, yTE.fixed-dy, yTE.fixed, yTE.fixed+dy),
            col="black")
    
    mtext(text.fixed,
          side=2, line=-0.2, outer=FALSE,
          at=yTE.fixed+0.25, adj=0, cex=cex.comb)
  }
  ##
  if (comb.random & overall){
    ##
    ## RE estimate
    ##
    polygon(c(lowTE.random, TE.random, uppTE.random, TE.random),
            c(yTE.random, yTE.random-dy, yTE.random, yTE.random+dy),
            col="darkgray")
    ##
    mtext(text.random, 
          side=2, line=-0.2, outer=FALSE,
          at=yTE.random+0.25, adj=0, cex=cex.comb)
  }
  
  invisible(NULL)
}
