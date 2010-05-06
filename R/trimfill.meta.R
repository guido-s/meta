trimfill.meta <- function(x, left=NULL, ma.fixed=TRUE,
                          type="L", n.iter.max=50,
                          sm=NULL, studlab=NULL,
                          level=x$level, level.comb=x$level.comb,
                          comb.fixed=x$comb.fixed, comb.random=x$comb.random,
                          silent=TRUE, ...){
  
  
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")
  if (inherits(x, "metacum"))
    stop("This function is not usable for an object of class \"metacum\"")
  if (inherits(x, "metainf"))
    stop("This function is not usable for an object of class \"metainf\"")
  
  
  if (length(comb.fixed)==0){
    comb.fixed <- TRUE
  }
  ##
  if (length(comb.random)==0){
    comb.random <- TRUE
  }
  
  
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
  
  
  TE <- x$TE
  seTE <- x$seTE
  sm <- x$sm
  studlab <- x$studlab
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
  ##
  k <- length(TE)
  
  
  if (match(type, c("L", "R"), nomatch=0) == 0)
    stop("type must be either 'L' or 'R'")
  
  
  if (is.null(left))
    left <- as.logical(sign(metabias(TE, seTE, meth="linreg")$estimate[1])==1)
  ##
  if (!left) TE <- -TE
  ##
  ord <- order(TE)
  ##print(data.frame(TE, studlab)[ord,])
  TE <- TE[ord]
  seTE <- seTE[ord]
  studlab <- studlab[ord]
  
  if (ma.fixed)
    TE.sum <- metagen(TE, seTE)$TE.fixed
  else
    TE.sum <- metagen(TE, seTE)$TE.random
  
  
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
        TE.sum <- metagen(TE[sel], seTE[sel])$TE.random
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
                 level=level, level.comb=level.comb)
  else
    m <- metagen(TE, seTE, studlab=studlab,
                 level=level, level.comb=level.comb)
  
  ##
  res <- list(studlab=m$studlab,
              TE=m$TE, seTE=m$seTE,
              w.fixed=m$w.fixed, w.random=m$w.random,
              TE.fixed=m$TE.fixed, seTE.fixed=m$seTE.fixed,
              TE.random=m$TE.random, seTE.random=m$seTE.random,
              k=m$k, Q=m$Q, tau=m$tau,
              sm=sm,
              method=m$method,
              ##paste("Inverse variance method (Trim and fill -",
              ##      ifelse(ma.fixed, "FE model)", "RE model)")),
              call=match.call(),
              left=left,
              ma.fixed=ma.fixed,
              type=type,
              n.iter.max=n.iter.max,
              n.iter=n.iter,
              trimfill=trimfill,
              k0=sum(trimfill),
              level=level, level.comb=level.comb,
              comb.fixed=comb.fixed, comb.random=comb.random)
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metagen", "meta", "trimfill")
  ##
  res
}
