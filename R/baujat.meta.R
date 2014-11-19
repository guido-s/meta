baujat.meta <- function(x,
                        yscale=1,
                        xlim, ylim,
                        xlab="Contribution to overall heterogeneity",
                        ylab="Influence on overall result",
                        pch=21, cex=1, col="black", bg="darkgray",
                        studlab=TRUE, cex.studlab=0.8,
                        xmin=0, ymin=0, pos=2, offset=0.5,
                        grid=TRUE, col.grid="lightgray", lty.grid="dotted", lwd.grid=par("lwd"),
                        pty="s",
                        ...){
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(x, "meta")
  ##
  if (inherits(x, "metacum"))
    stop("Baujat plot not meaningful for object of class \"metacum\"")
  if (inherits(x, "metainf"))
    stop("Baujat plot not meaningful for object of class \"metainf\"")
  if (inherits(x, "trimfill"))
    stop("Baujat plot not meaningful for object of class \"trimfill\"")
  ##
  x <- updateversion(x)
  
  
  oldpar <- par(pty=pty)
  on.exit(par(oldpar))
  
  
  TE <- x$TE
  seTE <- x$seTE
  TE.fixed <- metagen(TE, seTE)$TE.fixed
  k <- x$k

  if (is.logical(studlab)){
    if (studlab)
      studlab <- x$studlab
    else
      studlab <- rep("", length(TE))
  }
  else{
    studlab <- as.character(studlab)
    if (length(studlab) != length(TE))
      stop("Length of argument 'studlab' must be the same as number of studies in meta-analysis.")
  }
  
  
  m.inf <- metainf(x, pooled="fixed")
  TE.inf <- m.inf$TE[1:m.inf$k]
  seTE.inf <- m.inf$seTE[1:m.inf$k]
  ##
  ys <- (TE.inf - TE.fixed)^2 / seTE.inf^2
  ys <- ys*yscale
  ##  
  xs <- (TE - TE.fixed)^2 / seTE^2
  
  
  if (missing(xlim))
    xlim <- c(0, max(xs, na.rm=TRUE))
  ##
  if (missing(ylim))
    ylim <- c(0, max(ys, na.rm=TRUE))
  
  
  ##
  ## Do not print labels for studies with x and/or y values below
  ## limits
  ##
  if (!missing(xmin) & !missing(ymin))
    studlab[xs<xmin & ys<ymin] <- ""
  else if (!missing(xmin) & missing(ymin))
    studlab[xs<xmin] <- ""
  else if (missing(xmin) & !missing(ymin))
    studlab[ys<ymin] <- ""
  
  
  plot(xs, ys,
       xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
       type="n")
  ##
  if (grid)
    grid(col=col.grid, lty=lty.grid, lwd=lwd.grid)
  ##
  points(xs, ys, pch=pch, cex=cex, col=col, bg=bg)
  ##  
  text(xs, ys, labels=studlab, cex=cex.studlab, pos=pos, offset=offset)
  
  
  res <- data.frame(x=xs, y=ys)
  
  
  invisible(res)
}
