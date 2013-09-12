trimfill.default <- function(x, seTE, left=NULL, ma.fixed=TRUE,
                             type="L", n.iter.max=50,
                             sm=NULL, studlab=NULL,
                             level=0.95, level.comb=level,
                             comb.fixed=FALSE, comb.random=TRUE,
                             hakn=FALSE,
                             method.tau="DL",
                             prediction=FALSE, level.predict=level,
                             silent=TRUE, ...){
  
  estimate.missing <- function(TE, TE.sum, type){
    ##
    ## 1. Centre around mean
    ##
    TE.c <- TE - TE.sum
    n <- length(TE.c)
    ##
    ## 2. Rank absolute values of centred values
    ##
    r.star <- rank(abs(TE.c))*sign(TE.c)
    ##
    if (type=="L"){
      ##
      ## 3. Sum the positive ranks only
      ##
      S.rank <- sum(r.star[r.star>0])
      ##
      ## 4. Estimate for L0
      ##
      res0 <- (4*S.rank - n*(n+1))/(2*n-1)
      res0.plus <- max(0, res0 + 0.5) %/% 1
    }
    if (type=="R"){
      ##
      ## 5. Estimate for R0
      ##
      res0 <- n - abs(min(r.star)) - 1.5
      res0.plus <- max(0, res0 + 0.5) %/% 1
    }
    ##
    res <- list(res0=res0, res0.plus=res0.plus)
    res
  }
  
  
  TE <- x
  ##
  if (is.null(sm)) sm <- ""
  if (is.null(studlab)) studlab <- seq(along=x)
  data.name <- paste(deparse(substitute(x)),
                     deparse(substitute(seTE)),
                     sep=", ")
  
  
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
  
  
  if (match(type, c("L", "R"), nomatch=0) == 0)
    stop("type must be either 'L' or 'R'")
  
  
  ##
  ## Check for levels of confidence interval
  ##
  if (!is.numeric(level) | length(level)!=1)
    stop("parameter 'level' must be a numeric of length 1")
  if (level <= 0 | level >= 1)
    stop("parameter 'level': no valid level for confidence interval")
  ##
  if (!is.numeric(level.comb) | length(level.comb)!=1)
    stop("parameter 'level.comb' must be a numeric of length 1")
  if (level.comb <= 0 | level.comb >= 1)
    stop("parameter 'level.comb': no valid level for confidence interval")
  ##
  if (!is.numeric(level.predict) | length(level.predict)!=1)
    stop("parameter 'level.predict' must be a numeric of length 1")
  if (level.predict <= 0 | level.predict >= 1)
    stop("parameter 'level.predict': no valid level for confidence interval")
  
  
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
  res <- list(studlab=m$studlab,
              TE=m$TE, seTE=m$seTE,
              w.fixed=m$w.fixed, w.random=m$w.random,
              TE.fixed=m$TE.fixed, seTE.fixed=m$seTE.fixed,
              TE.random=m$TE.random, seTE.random=m$seTE.random,
              ##
              seTE.predict=m$seTE.predict,
              lower.predict=m$lower.predict,
              upper.predict=m$upper.predict,
              level.predict=level.predict,
              ##
              k=m$k, Q=m$Q, tau=m$tau,
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
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metagen", "meta", "trimfill")
  ##
  res
}
