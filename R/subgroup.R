subgroup <- function(x, tau.preset=NULL){
  
  
  byvar <- x$byvar
  ##
  bylevs <- bylevs(byvar)
  
  
  if (!(length(byvar)>0)){
    warning("Argument 'byvar' is missing.")
    return(NULL)
  }
  
  
  bin.cont <- inherits(x, "metabin") | inherits(x, "metacont")
  bin.inc <- inherits(x, "metabin") | inherits(x, "metainc")
  prop <- inherits(x, "metaprop")
  cor.prop <- inherits(x, "metacor") | inherits(x, "metaprop")
  
  
  sumNA <- function(x)
    if (all(is.na(x)))
      NA
    else sum(x, na.rm=TRUE)
  
  
  res.w <- matrix(NA, ncol=24, nrow=length(bylevs))
    j <- 0
  ##
  for (i in bylevs){
    j <- j+1
    sel <- byvar == i
    ##
    if (all(is.na(x$studlab[sel])))
      stop("No data available for byvar = ", i)
    ##
    ##
    if (inherits(x, "metabin"))
      meta1 <- metabin(x$event.e[sel], x$n.e[sel],
                       x$event.c[sel], x$n.c[sel],
                       studlab=x$studlab[sel],
                       method=x$method,
                       sm=x$sm,
                       incr=x$incr,
                       allincr=x$allincr,
                       addincr=x$addincr,
                       allstudies=x$allstudies,
                       MH.exact=x$MH.exact,
                       RR.cochrane=x$RR.cochrane,
                       level=x$level, level.comb=x$level.comb,
                       hakn=x$hakn,
                       method.tau=x$method.tau,
                       tau.preset=tau.preset,
                       TE.tau=x$TE.tau,
                       warn=x$warn)
    ##
    if (inherits(x, "metacont"))
      meta1 <- metacont(x$n.e[sel], x$mean.e[sel],
                        x$sd.e[sel],
                        x$n.c[sel], x$mean.c[sel],
                        x$sd.c[sel],
                        sm=x$sm, pooledvar=x$pooledvar,
                        studlab=x$studlab[sel],
                        level=x$level, level.comb=x$level.comb,
                        hakn=x$hakn,
                        method.tau=x$method.tau,
                        tau.preset=tau.preset, TE.tau=x$TE.tau,
                        warn=x$warn)
    ##
    if (inherits(x, "metagen"))
      meta1 <- metagen(x$TE[sel], x$seTE[sel],
                       sm=x$sm,
                       studlab=x$studlab[sel],
                       level=x$level, level.comb=x$level.comb,
                       hakn=x$hakn,
                       method.tau=x$method.tau,
                       tau.preset=tau.preset, TE.tau=x$TE.tau,
                       warn=x$warn)
    ##
    if (inherits(x, "metaprop"))
      meta1 <- metaprop(x$event[sel], x$n[sel],
                        sm=x$sm,
                        studlab=x$studlab[sel],
                        level=x$level, level.comb=x$level.comb,
                        incr=x$incr,
                        allincr=x$allincr,
                        addincr=x$addincr,
                        hakn=x$hakn,
                        method.tau=x$method.tau,
                        tau.preset=tau.preset, TE.tau=x$TE.tau,
                        warn=x$warn)
    ##
    if (inherits(x, "metacor"))
      meta1 <- metacor(x$cor[sel], x$n[sel],
                       sm=x$sm,
                       studlab=x$studlab[sel],
                       level=x$level, level.comb=x$level.comb,
                       hakn=x$hakn,
                       method.tau=x$method.tau,
                       tau.preset=tau.preset, TE.tau=x$TE.tau)
    ##
    if (inherits(x, "metainc"))
      meta1 <- metainc(x$event.e[sel], x$time.e[sel],
                       x$event.c[sel], x$time.c[sel],
                       studlab=x$studlab[sel],
                       method=x$method,
                       sm=x$sm,
                       incr=x$incr,
                       allincr=x$allincr,
                       addincr=x$addincr,
                       level=x$level, level.comb=x$level.comb,
                       hakn=x$hakn,
                       method.tau=x$method.tau,
                       tau.preset=tau.preset,
                       TE.tau=x$TE.tau,
                       warn=x$warn)
    ##
    ##
    res.w[j,] <- c(meta1$TE.fixed,                                      #  1
                   meta1$seTE.fixed,                                    #  2
                   meta1$Q,                                             #  3
                   meta1$k,                                             #  4
                   length(meta1$TE),                                    #  5
                   meta1$TE.random,                                     #  6
                   meta1$seTE.random,                                   #  7
                   meta1$H,                                             #  8
                   meta1$lower.H,                                       #  9
                   meta1$upper.H,                                       # 10
                   meta1$I2,                                            # 11
                   meta1$lower.I2,                                      # 12
                   meta1$upper.I2,                                      # 13
                   meta1$tau,                                           # 14
                   meta1$C,                                             # 15
                   mean(1/x$n[sel]),                                    # 16
                   sum(x$w.fixed[sel]),                                 # 17
                   sum(x$w.random[sel]),                                # 18
                   if (bin.inc) sumNA(meta1$event.e) else NA, # 19
                   if (bin.cont) sumNA(meta1$n.e) else NA,    # 20
                   if (bin.inc) sumNA(meta1$event.c) else NA, # 21
                   if (bin.cont) sumNA(meta1$n.c) else NA,    # 22
                   if (prop) sumNA(meta1$event) else NA,      # 23
                   if (cor.prop) sumNA(meta1$n) else NA       # 24
                   )
  }
  ##
  TE.fixed.w    <- res.w[,1]
  seTE.fixed.w  <- res.w[,2]
  Q.w           <- res.w[,3]
  k.w           <- res.w[,4]
  k.all.w       <- res.w[,5]
  TE.random.w   <- res.w[,6]
  seTE.random.w <- res.w[,7]
  ##
  H.w     <- res.w[,8]
  H.w.low <- res.w[,9]
  H.w.upp <- res.w[,10]
  ##
  I2.w     <- res.w[,11]
  I2.w.low <- res.w[,12]
  I2.w.upp <- res.w[,13]
  ##
  tau.w <- res.w[,14]
  ##
  C.w   <- res.w[,15]
  ##
  n.harmonic.mean.w <- res.w[,16]
  ##
  w.fixed.w <- res.w[,17]
  w.random.w <- res.w[,18]
  ##
  event.e.w <- res.w[,19]
  n.e.w <- res.w[,20]
  event.c.w <- res.w[,21]
  n.c.w <- res.w[,22]
  event.w <- res.w[,23]
  n.w <- res.w[,24]
  ##
  ci.fixed.w  <- ci(TE.fixed.w, seTE.fixed.w, x$level.comb)
  ##
  if (!is.null(x$hakn) && x$hakn)
    ci.random.w <- ci(TE.random.w, seTE.random.w, x$level.comb, df=k.w-1)
  else
    ci.random.w <- ci(TE.random.w, seTE.random.w, x$level.comb)
  ##
  Q.b.fixed  <- metagen(TE.fixed.w, seTE.fixed.w)$Q
  Q.b.random <- metagen(TE.random.w, seTE.random.w)$Q
  
  
  res <- list(bylevs=bylevs,
              TE.fixed.w=ci.fixed.w$TE,
              seTE.fixed.w=ci.fixed.w$seTE,
              lower.fixed.w=ci.fixed.w$lower,
              upper.fixed.w=ci.fixed.w$upper,
              zval.fixed.w=ci.fixed.w$z,
              pval.fixed.w=ci.fixed.w$p,
              w.fixed.w=w.fixed.w,
              ##
              TE.random.w=ci.random.w$TE,
              seTE.random.w=ci.random.w$seTE,
              lower.random.w=ci.random.w$lower,
              upper.random.w=ci.random.w$upper,
              zval.random.w=ci.random.w$z,
              pval.random.w=ci.random.w$p,
              df.hakn.w=ci.random.w$df,
              w.random.w=w.random.w,
              ##
              n.harmonic.mean.w=n.harmonic.mean.w,
              ##
              event.e.w=event.e.w,
              n.e.w=n.e.w,
              event.c.w=event.c.w,
              n.c.w=n.c.w,
              n.w=n.w,
              event.w=event.w,
              ##
              k.w=k.w,
              k.all.w=k.all.w,
              Q.w=Q.w,
              Q.w.fixed=sum(Q.w, na.rm=TRUE),
              Q.w.random=NA,
              df.Q.w=sum((k.w-1)[!is.na(Q.w)]),
              Q.b.fixed=Q.b.fixed,
              Q.b.random=Q.b.random,
              df.Q.b=x$k-1 - sum((k.w-1)[!is.na(Q.w)]),
              tau.w=tau.w,
              C.w=C.w,
              ##
              H.w=H.w,
              lower.H.w=H.w.low,
              upper.H.w=H.w.upp,
              ##
              I2.w=I2.w,
              lower.I2.w=I2.w.low,
              upper.I2.w=I2.w.upp
              )

  
  res
}
