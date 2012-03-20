metacum <- function(x, pooled, sortvar, level.comb=x$level.comb){
  
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
  
  
  if (length(level.comb)==0){
    warning("level.comb set to 0.95")
    level.comb <- 0.95
  }
  
  
  k.all <- length(x$TE)
  ##
  if (k.all==1){
    warning("Nothing calculated (minimum number of studies: 2)")
    return(invisible(NULL))
  }
  
  sort <- !missing(sortvar)
  
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

  
  if (pooled == "fixed" | (pooled == "random" & !x$hakn))
    res.i <- matrix(NA, ncol=8, nrow=k.all)
  ##
  else if (pooled == "random" & x$hakn)
    res.i <- matrix(NA, ncol=9, nrow=k.all)
  ##
  for (i in 1:k.all){
    sel <- 1:i
    ##
    if (inherits(x, "metabin"))
      m <- metabin(event.e[sel], n.e[sel], event.c[sel], n.c[sel],
                   method=x$method, sm=x$sm,
                   incr=x$incr, allincr=x$allincr, addincr=x$addincr,
                   allstudies=x$allstudies, MH.exact=x$MH.exact,
                   RR.cochrane=x$RR.cochrane,
                   level=level.comb, level.comb=level.comb,
                   hakn=x$hakn,
                   method.tau=x$method.tau,
                   tau.preset=x$tau.preset, TE.tau=x$TE.tau,
                   warn=FALSE)
    ##
    if (inherits(x, "metacont"))
      m <- metacont(n.e[sel], mean.e[sel], sd.e[sel],
                    n.c[sel], mean.c[sel], sd.c[sel], sm=x$sm,
                    level=level.comb, level.comb=level.comb,
                    hakn=x$hakn,
                    method.tau=x$method.tau,
                    tau.preset=x$tau.preset, TE.tau=x$TE.tau)
    ##
    if (inherits(x, "metagen"))
      m <- metagen(TE[sel], seTE[sel], sm=x$sm,
                   level=level.comb, level.comb=level.comb,
                   hakn=x$hakn,
                   method.tau=x$method.tau,
                   tau.preset=x$tau.preset, TE.tau=x$TE.tau)
    ##
    if (inherits(x, "metaprop"))
      m <- metaprop(event[sel], n[sel],
                    studlab=studlab[sel],
                    sm=x$sm,
                    level=level.comb, level.comb=level.comb,
                    hakn=x$hakn,
                    method.tau=x$method.tau,
                    tau.preset=x$tau.preset, TE.tau=x$TE.tau,
                    warn=FALSE)
    ##
    if (inherits(x, "metacor"))
      m <- metacor(cor[sel], n[sel],
                   studlab=studlab[sel],
                   sm=x$sm,
                   level=level.comb, level.comb=level.comb,
                   hakn=x$hakn,
                   method.tau=x$method.tau,
                   tau.preset=x$tau.preset, TE.tau=x$TE.tau)
    ##
    s.i <- summary(m, level=level.comb, level.comb=level.comb)
    ##
    if (pooled == "fixed"){
      res.i[i,] <- c(m$TE.fixed, m$seTE.fixed,
                     s.i$fixed$lower, s.i$fixed$upper,
                     m$pval.fixed, s.i$I2$TE,
                     m$tau, sum(m$w.fixed, na.rm=TRUE))
    }
    ##
    else if (pooled == "random" & !x$hakn){
      res.i[i,] <- c(m$TE.random, m$seTE.random,
                     s.i$random$lower, s.i$random$upper,
                     m$pval.random, s.i$I2$TE,
                     m$tau, sum(m$w.random, na.rm=TRUE))
    }
    ##
    else if (pooled == "random" & x$hakn){
      res.i[i,] <- c(m$TE.random, m$seTE.random,
                     s.i$random$lower, s.i$random$upper,
                     m$pval.random, s.i$I2$TE,
                     m$tau, sum(m$w.random, na.rm=TRUE),
                     m$df.hakn)
    }
  }
  ##
  TE.i <- res.i[,1]
  seTE.i <- res.i[,2]
  lower.i <- res.i[,3]
  upper.i <- res.i[,4]
  pval.i <- res.i[,5]
  I2.i <- res.i[,6]
  tau.i <- res.i[,7]
  weight.i <- res.i[,8]
  if (pooled == "random" & x$hakn)
    df.hakn.i <- res.i[,9]
  
  
  sm1 <- summary(x, level=level.comb, level.comb=level.comb)
  ##
  if (pooled == "fixed"){
    TE.s <- sm1$fixed$TE
    seTE.s <- sm1$fixed$seTE
    TE.s.lower <- sm1$fixed$lower
    TE.s.upper <- sm1$fixed$upper
    pval.s <- sm1$fixed$p
    w.s <- sum(x$w.fixed, na.rm=TRUE)
  }
  ##
  else if (pooled == "random"){
    TE.s <- sm1$random$TE
    seTE.s <- sm1$random$seTE
    TE.s.lower <- sm1$random$lower
    TE.s.upper <- sm1$random$upper
    pval.s <- sm1$random$p
    w.s <- sum(x$w.random, na.rm=TRUE)
  }
  
  
  slab <- character(k.all)
  for (i in 1:k.all)
    slab[i] <- paste("Adding ", studlab[i],
                     " (k=", i, ")", sep="")
  ##
  slab <- c(slab, "Pooled estimate")
  
  res <- list(TE=c(TE.i, NA, TE.s),
              seTE=c(seTE.i, NA, seTE.s),
              lower=c(lower.i, NA, TE.s.lower),
              upper=c(upper.i, NA, TE.s.upper),
              studlab=c(rev(rev(slab)[-1]), " ", rev(slab)[1]),
              p.value=c(pval.i, NA, pval.s),
              w=c(weight.i, NA, w.s),
              I2=c(I2.i, NA, sm1$I2$TE),
              tau=c(tau.i, NA, sm1$tau),
              df.hakn=if (pooled=="random" & x$hakn) c(df.hakn.i, NA, x$df.hakn) else NULL,
              sm=x$sm, method=x$method, k=x$k,
              pooled=pooled,
              TE.fixed=NA, seTE.fixed=NA,
              TE.random=NA, seTE.random=NA,
              Q=NA,
              level.comb=level.comb,
              hakn=x$hakn,
              method.tau=x$method.tau,
              tau.preset=x$tau.preset,
              TE.tau=x$TE.tau)
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metacum", "meta")
  
  res
}
