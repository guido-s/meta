metacum <- function(x, pooled, sortvar){
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(x, "meta")
  x <- updateversion(x)
  ##
  k.all <- length(x$TE)
  if (k.all < 2){
    warning("Nothing calculated (minimum number of studies: 2).")
    return(invisible(NULL))
  }
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  if (!missing(pooled))
    pooled <- setchar(pooled, c("fixed", "random"))
  else
    if (!x$comb.fixed & x$comb.random)
      pooled <- "random"
    else
      pooled <- "fixed"
  ##
  mf <- match.call()
  error <- try(sortvar <- eval(mf[[match("sortvar", names(mf))]],
                               as.data.frame(x, stringsAsFactors=FALSE),
                               enclos = sys.frame(sys.parent())),
               silent=TRUE)
  if (class(error)=="try-error"){
    xd <- x$data
    sortvar <- eval(mf[[match("sortvar", names(mf))]],
                    xd, enclos = NULL)
    if (!is.null(x$data$.subset))
      sortvar <- sortvar[x$data$.subset]
  }
  sort <- !is.null(sortvar)
  if (sort && (length(sortvar) != k.all))
    stop("Number of studies in object 'x' and argument 'sortvar' have different length.")
  if (!sort)
    sortvar <- 1:k.all
  
  
  ##
  ##
  ## (3) Sort variables
  ##
  ##
  o <- order(sortvar)
  ##
  n.e <- x$n.e[o]
  n.c <- x$n.c[o]
  n   <- x$n[o]
  ##
  event.e <- x$event.e[o]
  event.c <- x$event.c[o]
  event   <- x$event[o]
  ##
  mean.e <- x$mean.e[o]
  mean.c <- x$mean.c[o]
  ##
  sd.e <- x$sd.e[o]
  sd.c <- x$sd.c[o]
  ##
  time.e <- x$time.e[o]
  time.c <- x$time.c[o]
  ##
  cor <- x$cor[o]
  ##
  TE <- x$TE[o]
  seTE <- x$seTE[o]
  ##
  studlab <- x$studlab[o]
  slab <- character(k.all)
  for (i in 1:k.all)
    slab[i] <- paste("Adding ", studlab[i],
                     " (k=", i, ")", sep="")
  slab <- c(slab, "Pooled estimate")
  studlab <- c(rev(rev(slab)[-1]), " ", rev(slab)[1])
  
  
  ##
  ##
  ## (4) Do sensitivity analysis
  ##
  ##
  res.i <- matrix(NA, ncol=10, nrow=k.all)
  ##
  for (i in 1:k.all){
    sel <- 1:i
    ##
    if (inherits(x, "metabin"))
      m <- metabin(event.e[sel], n.e[sel], event.c[sel], n.c[sel],
                   ##
                   method=x$method, sm=x$sm,
                   incr=x$incr, allincr=x$allincr, addincr=x$addincr,
                   allstudies=x$allstudies, MH.exact=x$MH.exact,
                   RR.cochrane=x$RR.cochrane,
                   ##
                   level.comb=x$level.comb,
                   ##
                   hakn=x$hakn,
                   method.tau=x$method.tau,
                   tau.preset=x$tau.preset, TE.tau=x$TE.tau,
                   ##
                   keepdata=FALSE,
                   warn=FALSE)
    ##
    if (inherits(x, "metacont"))
      m <- metacont(n.e[sel], mean.e[sel], sd.e[sel],
                    n.c[sel], mean.c[sel], sd.c[sel],
                    ##
                    sm=x$sm, pooledvar=x$pooledvar,
                    ##
                    level.comb=x$level.comb,
                    ##
                    hakn=x$hakn,
                    method.tau=x$method.tau,
                    tau.preset=x$tau.preset, TE.tau=x$TE.tau,
                    ##
                    keepdata=FALSE,
                    warn=FALSE)
    ##
    if (inherits(x, "metacor"))
      m <- metacor(cor[sel], n[sel],
                   ##
                   sm=x$sm,
                   ##
                   level.comb=x$level.comb,
                   ##
                   hakn=x$hakn,
                   method.tau=x$method.tau,
                   tau.preset=x$tau.preset, TE.tau=x$TE.tau,
                   ##
                   keepdata=FALSE)
    ##
    if (inherits(x, "metagen"))
      m <- metagen(TE[sel], seTE[sel],
                   ##
                   sm=x$sm,
                   ##
                   level.comb=x$level.comb,
                   ##
                   hakn=x$hakn,
                   method.tau=x$method.tau,
                   tau.preset=x$tau.preset, TE.tau=x$TE.tau,
                   ##
                   keepdata=FALSE,
                   warn=FALSE)
    ##
    if (inherits(x,"metainc"))
      m <- metainc(event.e[sel], time.e[sel],
                   event.c[sel], time.c[sel],
                   ##
                   method=x$method,
                   sm=x$sm,
                   incr=x$incr, allincr=x$allincr, addincr=x$addincr,
                   ##
                   level.comb=x$level.comb,
                   ##
                   hakn=x$hakn, method.tau=x$method.tau,
                   tau.preset=x$tau.preset, TE.tau=x$TE.tau,
                   ##
                   keepdata=FALSE,
                   warn=FALSE)
    ##
    if (inherits(x, "metaprop"))
      m <- metaprop(event[sel], n[sel],
                    ##
                    sm=x$sm,
                    incr=x$incr, allincr=x$allincr, addincr=x$addincr,
                    method.ci=x$method.ci,
                    ##
                    level.comb=x$level.comb,
                    ##
                    hakn=x$hakn,
                    method.tau=x$method.tau,
                    tau.preset=x$tau.preset, TE.tau=x$TE.tau,
                    ##
                    keepdata=FALSE,
                    warn=FALSE)
    ##
    sel.pft <- inherits(x, "metaprop") & x$sm=="PFT"
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
  ##  
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
  
  
  ##
  ##
  ## (5) Generate R object
  ##
  ##
  res <- list(TE=c(TE.i, NA, TE.s),
              seTE=c(seTE.i, NA, seTE.s),
              lower=c(lower.i, NA, TE.s.lower),
              upper=c(upper.i, NA, TE.s.upper),
              studlab=studlab,
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
              level.comb=x$level.comb,
              hakn=x$hakn,
              method.tau=x$method.tau,
              tau.preset=x$tau.preset,
              TE.tau=x$TE.tau,
              n.harmonic.mean=c(n.harmonic.mean.i, NA, 1/mean(1/n)),
              prediction=FALSE,
              ##
              backtransf=x$backtransf,
              title=x$title, complab=x$complab,
              outclab=x$outclab,
              ##
              call=match.call())
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metacum", "meta")
  if (inherits(x, "trimfill"))
    class(res) <- c(class(res), "trimfill")
  
  res
}
