metainf <- function(x, pooled, sortvar, level=x$level, level.comb=x$level.comb){
  
  
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")
  
  
  if (missing(pooled)){
    if (length(x$comb.fixed)==0 & length(x$comb.random)==0)
      pooled <- "fixed"
    if (length(x$comb.fixed)>0 & length(x$comb.random)==0)
      if (x$comb.fixed) pooled <- "fixed"
      else pooled <- "NoNe"
    if (length(x$comb.fixed)==0 & length(x$comb.random)>0)
      if (x$comb.random) pooled <- "random"
      else pooled <- "NoNe"
    if (length(x$comb.fixed)>0 & length(x$comb.random)>0){
      if (x$comb.fixed)
        pooled <- "fixed"
      if (!x$comb.fixed & x$comb.random)
        pooled <- "random"
      if (!x$comb.fixed & !x$comb.random)
        pooled <- "NoNe"
    }
  }
  ##
  if (pooled=="NoNe")
    stop("Parameters \"comb.fixed\" and \"comb.random\" in object '",
        deparse(substitute(x)),
         "' are either 'FALSE' or 'NULL'. ",
         "Please use argument \"pooled=fixed\" or \"pooled=random\" ",
         "to select meta-analytical model.")
  
  
  imeth <- charmatch(tolower(pooled), c("fixed", "random"), nomatch = NA)
  ##
  if (is.na(imeth)) 
        stop("'pooled' should be \"fixed\" or \"random\"")
  ##
  pooled <- c("fixed", "random")[imeth]
  
  
  if (length(level)==0){
    warning("level set to 0.95")
    level <- 0.95
  }
  ##
  if (length(level.comb)==0){
    warning("level.comb set to 0.95")
    level.comb <- 0.95
  }
  
  
  k.all <- length(x$TE)
  sort <- !missing(sortvar)
  
  sm <- x$sm
  
  if (!sort) sortvar <- rep(1, k.all)

  if (sort & length(sortvar) != k.all)
    stop("'x' and 'sortvar' have different length")

  
  n.e <- x$n.e
  n.c <- x$n.c
  n   <- x$n
  ##
  event.e <- x$event.e
  event.c <- x$event.c
  event   <- x$event
  ##
  mean.e <- x$mean.e
  mean.c <- x$mean.c
  ##
  sd.e <- x$sd.e
  sd.c <- x$sd.c
  ##
  cor <- x$cor
  ##
  TE <- x$TE
  seTE <- x$seTE
  ##
  studlab <- x$studlab
  sortvar <- sortvar
  
  
  if (sort){
    ##
    o <- order(sortvar)
    ##
    n.e <- n.e[o]
    n.c <- n.c[o]
    n   <- n[o]
    ##
    event.e <- event.e[o]
    event.c <- event.c[o]
    event   <- event[o]
    ##
    mean.e <- mean.e[o]
    mean.c <- mean.c[o]
    ##
    sd.e <- sd.e[o]
    sd.c <- sd.c[o]
    ##
    cor <- cor[o]
    ##
    TE <- TE[o]
    seTE <- seTE[o]
    ##
    studlab <- studlab[o]
    sortvar <- sortvar[o]
  }

  
  res.i <- matrix(NA, ncol=6, nrow=k.all)
  ##
  for (i in 1:k.all){
    sel <- -i
    ##
    if (inherits(x, "metabin")){
      ##
      ## To get rid of warning message
      ## "For a single trial, inverse variance method used
      ##  instead of Mantel Haenszel method."
      ##
      oldopt <- options(warn=-1)
      m <- metabin(event.e[sel], n.e[sel], event.c[sel], n.c[sel],
                   method=x$method, sm=sm,
                   incr=x$incr, allincr=x$allincr, addincr=x$addincr,
                   allstudies=x$allstudies, MH.exact=x$MH.exact,
                   RR.cochrane=x$RR.cochrane,
                   level=level, level.comb=level.comb,
                   warn=x$warn)
      options(oldopt)
    }
    ##
    if (inherits(x, "metacont")){
      m <- metacont(n.e[sel], mean.e[sel], sd.e[sel],
                    n.c[sel], mean.c[sel], sd.c[sel], sm=sm,
                    level=level, level.comb=level.comb)
    }
    ##
    if (inherits(x, "metagen")){
      m <- metagen(TE[sel], seTE[sel], sm=sm,
                   level=level, level.comb=level.comb)
    }
    ##
    if (inherits(x, "metaprop")){
      m <- metaprop(event[sel], n[sel],
                    studlab=studlab[sel],
                    sm=sm,
                    level=level, level.comb=level.comb,
                    warn=x$warn)
    }
    ##
    if (inherits(x, "metacor")){
      m <- metacor(cor[sel], n[sel],
                   studlab=studlab[sel],
                   sm=sm,
                   level=level, level.comb=level.comb)
    }
    ##
    tsum.i <- summary(m)
    ##
    if (pooled == "fixed"){
      res.i[i,] <- c(m$TE.fixed, m$seTE.fixed,
                     tsum.i$fixed$p, tsum.i$I2$TE,
                     tsum.i$tau, sum(m$w.fixed, na.rm=TRUE))
    }
    ##
    else if (pooled == "random"){
      res.i[i,] <- c(m$TE.random, m$seTE.random,
                     tsum.i$random$p, tsum.i$I2$TE,
                     tsum.i$tau, sum(m$w.random, na.rm=TRUE))
    }
  }
  
  
  tsum <- summary(x, level=level, level.comb=level.comb)
  ##
  if (pooled == "fixed"){
    TE.sum <- x$TE.fixed
    seTE.sum <- x$seTE.fixed
    pval.sum <- tsum$fixed$p
    w.sum <- sum(x$w.fixed, na.rm=TRUE)
  }
  ##
  else if (pooled == "random"){
    TE.sum <- x$TE.random
    seTE.sum <- x$seTE.random
    pval.sum <- tsum$random$p
    w.sum <- sum(x$w.random, na.rm=TRUE)
  }
  
  
  slab <- c(paste("Omitting", studlab), "Pooled estimate")
  
  res <- list(TE=c(res.i[,1], NA, TE.sum),
              seTE=c(res.i[,2], NA, seTE.sum),
              studlab=c(rev(rev(slab)[-1]), " ", rev(slab)[1]),
              p.value=c(res.i[,3], NA, pval.sum),
              w=c(res.i[,6], NA, w.sum),
              I2=c(res.i[,4], NA, tsum$I2$TE),
              tau=c(res.i[,5], NA, tsum$tau),
              sm=sm, method=x$method, k=x$k,
              pooled=pooled,
              TE.fixed=NA, seTE.fixed=NA,
              TE.random=NA, seTE.random=NA,
              Q=NA, tau=NA,
              level=level, level.comb=level.comb)
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metainf", "meta")
  
  res
}
