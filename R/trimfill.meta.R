trimfill.meta <- function(x, left=NULL, ma.fixed=TRUE,
                          type="L", n.iter.max=50,
                          level=x$level, level.comb=x$level.comb,
                          comb.fixed=FALSE, comb.random=TRUE,
                          hakn=x$hakn,
                          method.tau=x$method.tau,
                          prediction=x$prediction, level.predict=x$level.predict,
                          backtransf=x$backtransf,
                          silent=TRUE, ...){
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(x, "meta")
  if (inherits(x, "metacum"))
    stop("This function is not usable for an object of class \"metacum\"")
  if (inherits(x, "metainf"))
    stop("This function is not usable for an object of class \"metainf\"")
  x <- updateversion(x)
  
  
  ##
  ## Check arguments
  ##
  type <- setchar(type, c("L", "R"))
  ##
  chklevel(level)
  chklevel(level.comb)
  chklevel(level.predict)
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  
  
  TE <- x$TE
  seTE <- x$seTE
  sm <- x$sm
  studlab <- x$studlab
  n.e <- x$n.e
  event.e <- x$event.e
  n.c <- x$n.c
  event.c <- x$event.c
  n <- x$n
  event <- x$event
  cor <- x$cor
  mean.e <- x$mean.e
  sd.e <- x$sd.e
  mean.c <- x$mean.c
  sd.c <- x$sd.c
  ##
  data.name <- deparse(substitute(x))
  
  
  if(length(TE) != length(seTE))
    stop("length of argument TE and seTE must be equal")
  ##
  if(length(TE) != length(studlab))
    stop("length of argument TE and studlab must be equal")
  ##
  sel <- !is.na(TE) & !is.na(seTE)
  if (length(TE) != sum(sel))
    warning(paste(length(TE) - sum(sel),
                  "observation(s) dropped due to missing values"))
  ##
  TE <- TE[sel]
  seTE <- seTE[sel]
  studlab <- studlab[sel]
  if (!is.null(n.e))
    n.e <- n.e[sel]
  if (!is.null(n.c))
    n.c <- n.c[sel]
  if (!is.null(n))
    n <- n[sel]
  if (!is.null(event.e))
    event.e <- event.e[sel]
  if (!is.null(event.c))
    event.c <- event.c[sel]
  if (!is.null(event))
    event <- event[sel]
  if (!is.null(cor))
    cor <- cor[sel]
  if (!is.null(mean.e))
    mean.e <- mean.e[sel]
  if (!is.null(mean.c))
    mean.c <- mean.c[sel]
  if (!is.null(sd.e))
    sd.e <- sd.e[sel]
  if (!is.null(sd.c))
    sd.c <- sd.c[sel]
  ##
  k <- length(TE)
  ##
  if (k<=2){
    warning("Minimal number of three studies for trim-and-fill method")
    return(invisible(NULL))
  }
  
  
  if (is.null(left))
    left <- as.logical(sign(metabias(TE, seTE, method="linreg", k.min=3)$estimate[1])==1)
  ##
  if (!left) TE <- -TE
  ##
  ord <- order(TE)
  TE <- TE[ord]
  seTE <- seTE[ord]
  studlab <- studlab[ord]
  if (!is.null(n.e))
    n.e <- n.e[ord]
  if (!is.null(n.c))
    n.c <- n.c[ord]
  if (!is.null(n))
    n <- n[ord]
  if (!is.null(event.e))
    event.e <- event.e[ord]
  if (!is.null(event.c))
    event.c <- event.c[ord]
  if (!is.null(event))
    event <- event[ord]
  if (!is.null(cor))
    cor <- cor[ord]
  if (!is.null(mean.e))
    mean.e <- mean.e[ord]
  if (!is.null(mean.c))
    mean.c <- mean.c[ord]
  if (!is.null(sd.e))
    sd.e <- sd.e[ord]
  if (!is.null(sd.c))
    sd.c <- sd.c[ord]
  
  
  if (ma.fixed)
    TE.sum <- metagen(TE, seTE)$TE.fixed
  else
    TE.sum <- metagen(TE, seTE, method.tau=method.tau)$TE.random
  
  
  if (k==1){
    n.iter <- 0
    k0 <- -9
  }
  else{
    n.iter  <-  0
    k0.last <- -1
    k0      <-  0
    ##
    while (k0.last != k0 & k0 <= (k-1) & n.iter < n.iter.max){
      ##
      n.iter <- n.iter + 1
      ##
      k0.last <- k0
      ##
      sel <- 1:(k-k0)
      ##
      if (ma.fixed)
        TE.sum <- metagen(TE[sel], seTE[sel])$TE.fixed
      else
        TE.sum <- metagen(TE[sel], seTE[sel],
                          method.tau=method.tau)$TE.random
      ##
      trim1 <- estimate.missing(TE, TE.sum, type)
      ##
      if (!silent){
        cat("n.iter = ", n.iter, "\n", sep="")
        if (type=="L")
          cat("L0 = ", round(trim1$res0, 2), "\n\n", sep="")
        if (type=="R")
          cat("R0 = ", round(trim1$res0+0.5, 2), "\n\n", sep="")
      }
      ##
      k0 <- trim1$res0.plus
    }
  }
  
  
  if (k0 > (k-1)) k0 <- k-1
  ##
  if (k0 > 0){
    TE.star   <- 2 * TE.sum - TE[(k-k0+1):k]
    seTE.star <- seTE[(k-k0+1):k]
    ##
    trimfill  <- c(rep(FALSE, length(TE)),
                   rep(TRUE, length(TE.star)))
    ##
    TE      <- c(TE[order(ord)], TE.star)
    seTE    <- c(seTE[order(ord)], seTE.star)
    studlab <- c(studlab[order(ord)],
                 paste("Filled:", studlab[(k-k0+1):k]))
    if (!is.null(n.e))
      n.e <- c(n.e[order(ord)], n.e[(k-k0+1):k])
    if (!is.null(n.c))
      n.c <- c(n.c[order(ord)], n.c[(k-k0+1):k])
    if (!is.null(n))
      n <- c(n[order(ord)], n[(k-k0+1):k])
    if (!is.null(event.e))
      event.e <- c(event.e[order(ord)], event.e[(k-k0+1):k])
    if (!is.null(event.c))
      event.c <- c(event.c[order(ord)], event.c[(k-k0+1):k])
    if (!is.null(event))
      event <- c(event[order(ord)], event[(k-k0+1):k])
    if (!is.null(cor))
      cor <- c(cor[order(ord)], cor[(k-k0+1):k])
    if (!is.null(mean.e))
      mean.e <- c(mean.e[order(ord)], mean.e[(k-k0+1):k])
    if (!is.null(mean.c))
      mean.c <- c(mean.c[order(ord)], mean.c[(k-k0+1):k])
    if (!is.null(sd.e))
      sd.e <- c(sd.e[order(ord)], sd.e[(k-k0+1):k])
    if (!is.null(sd.c))
      sd.c <- c(sd.c[order(ord)], sd.c[(k-k0+1):k])
  }
  else{
    TE.star   <- NA
    seTE.star <- NA
    trimfill  <- rep(FALSE, length(TE))
    TE        <- TE[order(ord)]
    seTE      <- seTE[order(ord)]
    studlab   <- studlab[order(ord)]
    if (!is.null(n.e))
      n.e <- n.e[order(ord)]
    if (!is.null(n.c))
      n.c <- n.c[order(ord)]
    if (!is.null(n))
      n <- n[order(ord)]
    if (!is.null(event.e))
      event.e <- event.e[order(ord)]
    if (!is.null(event.c))
      event.c <- event.c[order(ord)]
    if (!is.null(event))
      event <- event[order(ord)]
    if (!is.null(cor))
      cor <- cor[order(ord)]
    if (!is.null(mean.e))
      mean.e <- mean.e[order(ord)]
    if (!is.null(mean.c))
      mean.c <- mean.c[order(ord)]
    if (!is.null(sd.e))
      sd.e <- sd.e[order(ord)]
    if (!is.null(sd.c))
      sd.c <- sd.c[order(ord)]
  }
  
  
  if (!left)
    m <- metagen(-TE, seTE, studlab=studlab,
                 level=level, level.comb=level.comb,
                 hakn=hakn, method.tau=method.tau,
                 prediction=prediction, level.predict=level.predict)
  else
    m <- metagen(TE, seTE, studlab=studlab,
                 level=level, level.comb=level.comb,
                 hakn=hakn, method.tau=method.tau,
                 prediction=prediction, level.predict=level.predict)
  
  
  ##
  ## Calculate H and I-Squared
  ##
  Hres  <- calcH(m$Q, m$df.Q, level.comb)
  I2res <- isquared(m$Q, m$df.Q, level.comb)
  
  
  res <- list(studlab=m$studlab,
              TE=m$TE, seTE=m$seTE,
              lower=m$lower, upper=m$upper,
              zval=m$zval, pval=m$pval,
              w.fixed=m$w.fixed, w.random=m$w.random,
              TE.fixed=m$TE.fixed, seTE.fixed=m$seTE.fixed,
              lower.fixed=m$lower.fixed, upper.fixed=m$upper.fixed,
              zval.fixed=m$zval.fixed, pval.fixed=m$pval.fixed,
              ##
              TE.random=m$TE.random, seTE.random=m$seTE.random,
              lower.random=m$lower.random, upper.random=m$upper.random,
              zval.random=m$zval.random, pval.random=m$pval.random,
              ##
              seTE.predict=m$seTE.predict,
              lower.predict=m$lower.predict,
              upper.predict=m$upper.predict,
              level.predict=level.predict,
              ##
              k=m$k, Q=m$Q, df.Q=m$df.Q, tau=m$tau,
              ##
              H=Hres$TE,
              lower.H=Hres$lower,
              upper.H=Hres$upper,
              ##
              I2=I2res$TE,
              lower.I2=I2res$lower,
              upper.I2=I2res$upper,
              ##
              sm=sm,
              method=m$method,
              ##
              call=match.call(),
              left=left,
              ma.fixed=ma.fixed,
              type=type,
              n.iter.max=n.iter.max,
              n.iter=n.iter,
              trimfill=trimfill,
              hakn=m$hakn,
              df.hakn=m$df.hakn,
              method.tau=m$method.tau,
              prediction=prediction,
              title=x$title,
              complab=x$complab,
              outclab=x$outclab,
              label.e=x$label.e,
              label.c=x$label.c,
              label.left=x$label.left,
              label.right=x$label.right,
              k0=sum(trimfill),
              level=level, level.comb=level.comb,
              comb.fixed=comb.fixed,
              comb.random=comb.random,
              n.e=n.e,
              n.c=n.c,
              event.e=event.e,
              event.c=event.c,
              mean.e=mean.e,
              mean.c=mean.c,
              sd.e=sd.e,
              sd.c=sd.c,
              n=n,
              event=event,
              cor=cor,
              class.x=class(x)[1]
              )
  
  res$backtransf <- backtransf
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metagen", "meta", "trimfill")
  ##
  res
}
