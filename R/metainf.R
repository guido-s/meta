metainf <- function(x, pooled, sortvar){
  
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")
  
  
  ## Upgrade meta objects created with older versions of meta
  ##
  if (!(!is.null(x$version) &&
        as.numeric(unlist(strsplit(x$version, "-"))[1]) >= 3.8))
    x <- update(x, warn=FALSE)
  
  
  level.comb <- x$level.comb
  
  
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
    warning("Nothing calculated (minimum number of studies: 2).")
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
  time.e <- x$time.e
  time.c <- x$time.c
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
    time.e <- time.e[o]
    time.c <- time.c[o]
    ##
    cor <- cor[o]
    ##
    TE <- TE[o]
    seTE <- seTE[o]
    ##
    studlab <- studlab[o]
    sortvar <- sortvar[o]
  }
  
  
  res.i <- matrix(NA, ncol=10, nrow=k.all)
  ##
  for (i in 1:k.all){
    sel <- -i
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
                    n.c[sel], mean.c[sel], sd.c[sel],
                    sm=x$sm, pooledvar=x$pooledvar,
                    level=level.comb, level.comb=level.comb,
                    hakn=x$hakn,
                    method.tau=x$method.tau,
                    tau.preset=x$tau.preset, TE.tau=x$TE.tau)
    ##
    if (inherits(x, "metagen"))
      if (!is.null(x$tau.preset))
        m <- metagen(TE[sel], seTE[sel], sm=x$sm,
                     level=level.comb, level.comb=level.comb,
                     hakn=x$hakn,
                     method.tau=x$method.tau,
                     tau.preset=x$tau.preset, TE.tau=x$TE.tau)
      else
        m <- metagen(TE[sel], seTE[sel], sm=x$sm,
                     level=level.comb, level.comb=level.comb,
                     hakn=x$hakn,
                     method.tau=x$method.tau,
                     TE.tau=x$TE.tau)
    ##
    if (inherits(x, "metaprop"))
      m <- metaprop(event[sel], n[sel],
                    studlab=studlab[sel],
                    sm=x$sm,
                    incr=x$incr, allincr=x$allincr, addincr=x$addincr,
                    method.ci=x$method.ci,
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
    if (inherits(x,"metainc"))
      m <- metainc(event.e[sel], time.e[sel],
                   event.c[sel], time.c[sel],
                   studlab=studlab[sel],
                   method=x$method,
                   sm=x$sm,
                   incr=x$incr, allincr=x$allincr, addincr=x$addincr,
                   level=x$level, level.comb=x$level.comb,
                   hakn=x$hakn, method.tau=x$method.tau,
                   tau.preset=x$tau.preset, TE.tau=x$TE.tau)
    ##
    sel.pas <- inherits(x, "metaprop") & m$sm=="PAS"
    sel.pft <- inherits(x, "metaprop") & m$sm=="PFT"
    ##
    if (pooled == "fixed"){
      res.i[i,] <- c(m$TE.fixed, m$seTE.fixed,
                     m$lower.fixed, m$upper.fixed,
                     m$pval.fixed, m$I2,
                     m$tau, sum(m$w.fixed, na.rm=TRUE),
                     if (sel.pft) 1/mean(1/n[sel]) else NA,
                     NA)
    }
    ##
    else if (pooled == "random"){
      res.i[i,] <- c(m$TE.random, m$seTE.random,
                     m$lower.random, m$upper.random,
                     m$pval.random, m$I2,
                     m$tau, sum(m$w.random, na.rm=TRUE),
                     if (sel.pft) 1/mean(1/n[sel]) else NA,
                     if (x$hakn) m$df.hakn else NA)
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
  n.harmonic.mean.i <- res.i[,9]
  if (pooled == "random" & x$hakn)
    df.hakn.i <- res.i[,10]
  
  
  if (pooled == "fixed"){
    TE.s <- x$TE.fixed
    seTE.s <- x$seTE.fixed
    TE.s.lower <- x$lower.fixed
    TE.s.upper <- x$upper.fixed
    pval.s <- x$pval.fixed
    w.s <- sum(x$w.fixed, na.rm=TRUE)
  }
  ##
  else if (pooled == "random"){
    TE.s <- x$TE.random
    seTE.s <- x$seTE.random
    TE.s.lower <- x$lower.random
    TE.s.upper <- x$upper.random
    pval.s <- x$pval.random
    w.s <- sum(x$w.random, na.rm=TRUE)
  }
  
  
  slab <- c(paste("Omitting", studlab), "Pooled estimate")
  
  
  res <- list(TE=c(TE.i, NA, TE.s),
              seTE=c(seTE.i, NA, seTE.s),
              lower=c(lower.i, NA, TE.s.lower),
              upper=c(upper.i, NA, TE.s.upper),
              studlab=c(rev(rev(slab)[-1]), " ", rev(slab)[1]),
              p.value=c(pval.i, NA, pval.s),
              w=c(weight.i, NA, w.s),
              I2=c(I2.i, NA, x$I2),
              tau=c(tau.i, NA, x$tau),
              df.hakn=if (pooled=="random" & x$hakn) c(df.hakn.i, NA, x$df.hakn) else NULL,
              sm=x$sm, method=x$method, k=x$k,
              pooled=pooled,
              comb.fixed=ifelse(pooled=="fixed", TRUE, FALSE),
              comb.random=ifelse(pooled=="random", TRUE, FALSE),
              TE.fixed=NA, seTE.fixed=NA,
              TE.random=NA, seTE.random=NA,
              Q=NA,
              level.comb=level.comb,
              hakn=x$hakn,
              method.tau=x$method.tau,
              tau.preset=x$tau.preset,
              TE.tau=x$TE.tau,
              n.harmonic.mean=c(n.harmonic.mean.i, NA, 1/mean(1/n)),
              call=match.call())
  
  res$backtransf <- x$backtransf
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metainf", "meta")
  if (inherits(x, "trimfill"))
    class(res) <- c(class(res), "trimfill")
  
  res
}
