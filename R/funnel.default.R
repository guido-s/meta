funnel.default <- function(x, y,
                           ##
                           xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
                           ##
                           comb.fixed=FALSE, comb.random=FALSE,
                           ##
                           axes=TRUE,
                           pch=21, text=NULL, cex=1,
                           ##
                           lty.fixed=2, lty.random=9,
                           lwd=1, lwd.fixed=lwd, lwd.random=lwd,
                           col="black", bg="darkgray",
                           col.fixed="black", col.random="black",
                           ##
                           log="", yaxis="se", sm="",
                           contour.levels=NULL, col.contour,
                           ##
                           ref=ifelse(backtransf & is.relative.effect(sm), 1, 0),
                           ##
                           level=NULL,
                           studlab=FALSE, cex.studlab=0.8,
                           ##
                           backtransf=TRUE,
                           ...){
  
  
  ##
  ##
  ## (1) Check essential arguments
  ##
  ##
  TE <- x
  k.All <- length(TE)
  seTE <- y
  ##
  chknumeric(TE)
  chknumeric(seTE)
  ##
  fun <- "funnel"
  chklength(seTE, k.All, fun)
  
  
  ##
  ##
  ## (2) Do meta-analysis
  ##
  ##
  m <- metagen(TE, seTE, sm=sm)
  
  
  ##
  ##
  ## (3) Produce funnel plot
  ##
  ##
  res <- funnel(m,
                xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
                ##
                comb.fixed=comb.fixed, comb.random=comb.random,
                ##
                axes=axes,
                pch=pch, text=text, cex=cex,
                ##
                lty.fixed=lty.fixed, lty.random=lty.random,
                lwd=lwd, lwd.fixed=lwd.fixed, lwd.random=lwd.random,
                col=col, bg=bg,
                col.fixed=col.fixed, col.random=col.random,
                ##
                log=log, yaxis=yaxis,
                contour.levels=contour.levels, if (!missing(col.contour)) col.contour=col.contour,
                ##
                ref=ref,
                ##
                level=level,
                studlab=studlab, cex.studlab=cex.studlab,
                ##
                backtransf=backtransf,
                ...)
  
  
  invisible(res)
}
