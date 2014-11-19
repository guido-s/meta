trimfill.default <- function(x, seTE, left=NULL, ma.fixed=TRUE,
                             type="L", n.iter.max=50,
                             sm=NULL, studlab=NULL,
                             level=0.95, level.comb=level,
                             comb.fixed=FALSE, comb.random=TRUE,
                             hakn=FALSE,
                             method.tau="DL",
                             prediction=FALSE, level.predict=level,
                             backtransf=TRUE,
                             silent=TRUE, ...){
  
  
  TE <- x
  ##
  if (is.null(sm)) sm <- ""
  if (is.null(studlab)) studlab <- seq(along=x)
  data.name <- paste(deparse(substitute(x)),
                     deparse(substitute(seTE)),
                     sep=", ")
  
  
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
  ##
  k <- length(TE)
  ##
  if (k<=2){
    warning("Minimal number of three studies for trim-and-fill method.")
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
    TE        <- c(TE[order(ord)], TE.star)
    seTE      <- c(seTE[order(ord)], seTE.star)
    studlab   <- c(studlab[order(ord)],
                   paste("Filled:", studlab[(k-k0+1):k]))
  }
  else{
    TE.star   <- NA
    seTE.star <- NA
    trimfill  <- rep(FALSE, length(TE))
    TE        <- TE[order(ord)]
    seTE      <- seTE[order(ord)]
    studlab   <- studlab[order(ord)]
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
              k0=sum(trimfill),
              level=level, level.comb=level.comb,
              comb.fixed=comb.fixed, comb.random=comb.random)
  
  res$backtransf <- backtransf
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metagen", "meta", "trimfill")
  ##
  res
}
