bubble.metareg <- function(x,
                           xlim, ylim,
                           xlab, ylab,
                           cex, min.cex=0.5, max.cex=5,
                           pch=21, col="black", bg="darkgray",
                           lty=1, lwd=1, col.line="black",
                           studlab=FALSE, cex.studlab=0.8,
                           pos=2, offset=0.5,
                           regline=TRUE,
                           axes=TRUE, box=TRUE, ...){
  
  
  ##
  ##
  ## (1) Check for meta object
  ##
  ##
  chkclass(x, "metareg")
  
  
  m0 <- x$.meta$x
  method.tau0 <- x$.meta$method.tau
  ##
  if (method.tau0!="FE" & (method.tau0 != m0$method.tau))
    m1 <- update(m0, method.tau=method.tau0)
  else
    m1 <- m0
  ##
  TE <- m1$TE
  sm <- m1$sm
  
  
  if (is.logical(studlab)){
    if (studlab)
      studlab <- m1$studlab
    else
      studlab <- rep("", length(TE))
  }
  else{
    studlab <- as.character(studlab)
    if (length(studlab) != length(TE))
      stop("Length of argument 'studlab' must be the same as number of studies in meta-analysis.")
  }
  
  
  charform <- as.character(x$.meta$formula)[2]
  splitform <- strsplit(charform, " ")[[1]]
  covar.name <- splitform[1]
  if (covar.name=="1" | covar.name=="-1")
    covar.name <- splitform[3]
  ##
  covar.names <- names(coef(x))
  covar.names.without.intrcpt <- covar.names[covar.names!="intrcpt"]
  if (length(covar.names.without.intrcpt)==0){
    warning("No covariate in meta-regression.")
    return(invisible(NULL))
  }
  ##
  nointrcpt <- ifelse("intrcpt" %in% covar.names, FALSE, TRUE)
  ##
  if (covar.name==".byvar")
    covar.name <- x$.meta$x$bylab
  ##
  if (length(covar.names.without.intrcpt) > 1){
    warning(paste("Only first covariate in meta-regression ",
                  "('", covar.name, "') considered in bubble plot. No regression line plotted.",
                  sep="")
            )
    regline <- FALSE
    if (missing(xlab))
      xlab <- paste("Covariate ", covar.name,
                    " (meta-regression: ", charform, ")", sep="")
  }
  else
    if (missing(xlab))
      xlab=paste("Covariate", covar.name)
  ##
  covar <- x$.meta$x$data[,covar.name]
  ##
  if (!is.null(x$.meta$x$subset))
    covar <- covar[x$.meta$x$subset]
  ##
  if (is.character(covar))
    covar <- as.factor(covar)
  ##
  if (is.factor(covar)){
    levs <- levels(covar)
    xs <- as.numeric(covar) - 1
    at.x <- sort(unique(xs))
    if (missing(xlim))
      xlim <- c(-0.5, 1.5)
    if (missing(xlab))
      xlab <- ""
  }
  else{
    xs <- covar
    if (missing(xlim))
      xlim <- range(xs)
  }
  ##
  alpha <- ifelse(nointrcpt, 0, coef(x)["intrcpt"])
  beta <- coef(x)[covar.name]
  
  
  ys <- TE
  ##
  if (missing(ylim))
    ylim <- range(ys)
  ##
  if (missing(ylab))
    ylab <- paste("Treatment effect (",
                  tolower(xlab(sm, backtransf=FALSE)),
                  ")", sep="")
  
  
  missing.cex <- missing(cex)
  fixed.cex   <- !missing.cex && (length(cex)==1 & all(cex=="fixed"))
  ##
  if (missing.cex)
    if (method.tau0=="FE")
      cex <- m1$w.fixed
    else
      cex <- m1$w.random
  else if (fixed.cex)
    cex <- m1$w.fixed
  ##
  if (length(cex) != length(TE))
    stop("Length of argument 'cex' must be the same as number of studies in meta-analysis.")
  ##
  if (missing.cex | fixed.cex){
    cexs <- max.cex*(cex/max(cex))
    cexs[cexs<min.cex] <- min.cex
  }
  else
    cexs <- cex
  
  
  ##
  ## Generate bubble plot
  ##
  if (is.factor(covar)){
    for (i in 2:length(levs)){
      sel <- xs %in% c(0, i-1)
      xs.i <- xs[sel]
      xs.i[xs.i>0] <- 1
      ys.i <- ys[sel]
      studlab.i <- studlab[sel]
      cexs.i <- cexs[sel]
      ##
      plot(xs.i, ys.i,
       pch=pch, cex=cexs.i,
       xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
       type="n", axes=FALSE, ...)
      ##
      ## Add regression line
      ##
      if (regline)
        lines(c(0, 1),
              c(alpha, alpha + coef(x)[i]),
              lty=lty, lwd=lwd, col=col.line)
      ##
      points(xs.i, ys.i, cex=cexs.i, pch=pch, col=col, bg=bg)
      ##
      ## x-axis
      ##
      if (axes)
        axis(1, at=0:1, labels=levs[c(1,i)], ...)
      ##
      ## y-axis
      ##
      if (axes)
        axis(2, ...)
      ##  
      text(xs.i, ys.i, labels=studlab.i, cex=cex.studlab,
           pos=pos, offset=offset)
      ##
      if (box)
        box()
    }
  }
  else{
    plot(xs, ys,
         pch=pch, cex=cexs,
         xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
         type="n", axes=FALSE, ...)
    ##
    ## Add regression line
    ##
    if (regline)
      abline(alpha, beta,
             lty=lty, lwd=lwd, col=col.line)
    ##
    points(xs, ys, cex=cexs, pch=pch, col=col, bg=bg)
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
    text(xs, ys, labels=studlab, cex=cex.studlab,
         pos=pos, offset=offset)
    ##
    if (box)
      box()
  }
  
  
  invisible(NULL)
}
