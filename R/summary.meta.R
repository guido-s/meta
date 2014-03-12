summary.meta <- function(object,
                         comb.fixed=object$comb.fixed,
                         comb.random=object$comb.random,
                         prediction=object$prediction,
                         bylab=object$bylab,
                         print.byvar=object$print.byvar,
                         bystud=FALSE,
                         print.CMH=object$print.CMH,
                         warn=object$warn,
                         ...){
  
  if (!inherits(object, "meta"))
    stop("Argument 'object' must be an object of class \"meta\"")
  
  if (length(warn)==0){
    warn <- TRUE
  }
  
  if (warn){
    if (inherits(object, "metacum"))
      warning("Summary method not defined for objects of class \"metacum\".")
    ##
    if (inherits(object, "metainf"))
      warning("Summary method not defined for objects of class \"metainf\".")
  }
  
  
  k <- object$k
  Q <- object$Q
  ##
  if (is.null(object$df.Q))
    df.Q <- k-1
  else
    df.Q <- object$df.Q
  
  
  cl <- class(object)[1]
  addargs <- names(list(...))
  ##
  fun <- "summary.meta"
  ##
  warnarg("byvar", addargs, fun, cl)
  warnarg("level", addargs, fun, cl)
  warnarg("level.comb", addargs, fun, cl)
  warnarg("level.predict", addargs, fun, cl)
  ##
  byvar <- object$byvar
  level <- object$level
  level.comb <- object$level.comb
  level.predict <- object$level.predict
  
  
  if (length(comb.fixed)==0)
    comb.fixed <- TRUE
  ##
  if (length(comb.random)==0)
    comb.random <- TRUE
  ##
  if (length(prediction)==0)
    prediction <- FALSE
  ##
  if (length(print.byvar)==0)
    print.byvar <- TRUE
  ##
  if (length(print.CMH)==0)
    print.CMH <- FALSE
  ##
  if (length(object$tau.common)==0)
    object$tau.common <- FALSE
  
  if (length(level)==0){
    if (warn)
      warning("level set to 0.95")
    level <- 0.95
  }
  ##
  if (length(level.comb)==0){
    if ((comb.fixed | comb.random) & warn)
      warning("level.comb set to 0.95")
    level.comb <- 0.95
  }
  ##
  if (length(level.predict)==0){
    if (prediction & comb.random & warn)
      warning("level.predict set to 0.95")
    level.predict <- 0.95
  }
  
  
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
  
  
  ##
  ## Check for older version of R package meta
  ##
  oldmeta <- !(!is.null(object$version) &&
               as.numeric(unlist(strsplit(object$version, "-"))[1]) >= 3.2)
  ##
  if (oldmeta){
    ci.H <- calcH(Q, df.Q, level.comb)
    ci.I2 <- isquared(Q, df.Q, level.comb)
    }
  else{
    ci.H <- list(TE=object$H, lower=object$lower.H, upper=object$upper.H)
    ci.I2 <- list(TE=object$I2, lower=object$lower.I2, upper=object$upper.I2)
  }
  
  
  ci.lab <- paste(round(100*level.comb, 1), "%-CI", sep="")
  ##
  ci.study <- ci(object$TE, object$seTE, level)
  ##
  ## Check for very old versions of R package meta
  ##
  ancientmeta <- !(!is.null(object$version) &&
                   as.numeric(unlist(strsplit(object$version, "\\."))[1]) >= 2)
  ##
  if (ancientmeta){
    ci.f <- ci(object$TE.fixed , object$seTE.fixed , level.comb)
    ci.r <- ci(object$TE.random, object$seTE.random, level.comb, df=object$df.hakn)
  }
  else{
    ## Use available values
    ci.f <- list(TE=object$TE.fixed,
                 seTE=object$seTE.fixed,
                 lower=object$lower.fixed,
                 upper=object$upper.fixed,
                 z=object$zval.fixed,
                 p=object$pval.fixed,
                 level=object$level.comb)
    ##
    ci.r <- list(TE=object$TE.random,
                 seTE=object$seTE.random,
                 lower=object$lower.random,
                 upper=object$upper.random,
                 z=object$zval.random,
                 p=object$pval.random,
                 level=object$level.comb,
                 df=if (!is.null(object$df.hakn)) object$df.hakn else NA)
  }
  ##
  ## Calculate exact confidence intervals for individual studies
  ##
  if (!(inherits(object, "metainf")|inherits(object, "metacum")) &
      inherits(object, "metaprop")){
    for ( i in 1:length(ci.study$TE)){
      cint <- binom.test(object$event[i], object$n[i], conf.level=level)
      ci.study$TE[i]    <- cint$estimate
      ci.study$lower[i] <- cint$conf.int[[1]]
      ci.study$upper[i] <- cint$conf.int[[2]]
      ci.study$seTE[i]  <- NA
      ci.study$z[i]     <- NA
      ci.study$p[i]     <- NA
    }
    ci.f$z    <- NA
    ci.f$p    <- NA
    ci.f$harmonic.mean <- mean(1/object$n)
    ##
    ci.r$z    <- NA
    ci.r$p    <- NA
    ci.r$harmonic.mean <- mean(1/object$n)
  }
  
  
  if (length(byvar)>0){
    
    if (any(is.na(byvar))) stop("Missing values in 'byvar'")
    
    if (is.factor(byvar))
      by.levs <- levels(byvar)
    else
      by.levs <- unique(byvar)
    
    if (length(bylab)==0) bylab <- ""
    
    if (object$tau.common){
      if (!is.null(object$tau.preset) & warn)
        warning("Value for argument 'tau.preset' not considered as argument 'tau.common=TRUE'.")
      object$tau.preset <- object$tau
    }
    
    res.w <- matrix(NA, ncol=15, nrow=length(by.levs))
    j <- 0
    ##
    for (i in by.levs){
      j <- j+1
      sel <- byvar == i
      ##
      if (all(is.na(object$studlab[sel])))
        stop("No data available for byvar = ", i)
      ##
      ##
      if (inherits(object, "metabin")){
        meta1 <- metabin(object$event.e[sel], object$n.e[sel],
                         object$event.c[sel], object$n.c[sel],
                         studlab=object$studlab[sel],
                         method=object$method,
                         sm=object$sm,
                         incr=object$incr,
                         allincr=object$allincr,
                         addincr=object$addincr,
                         allstudies=object$allstudies,
                         MH.exact=object$MH.exact,
                         RR.cochrane=object$RR.cochrane,
                         level=level, level.comb=level.comb,
                         comb.fixed=comb.fixed,
                         comb.random=comb.random,
                         hakn=object$hakn,
                         method.tau=object$method.tau,
                         tau.preset=object$tau.preset,
                         TE.tau=object$TE.tau,
                         warn=warn)
      }
      ##
      if (inherits(object, "metacont")){
        meta1 <- metacont(object$n.e[sel], object$mean.e[sel],
                          object$sd.e[sel],
                          object$n.c[sel], object$mean.c[sel],
                          object$sd.c[sel],
                          sm=object$sm,
                          studlab=object$studlab[sel],
                          level=level, level.comb=level.comb,
                          comb.fixed=comb.fixed,
                          comb.random=comb.random,
                          hakn=object$hakn,
                          method.tau=object$method.tau,
                          tau.preset=object$tau.preset, TE.tau=object$TE.tau,
                          warn=warn)
      }
      ##
      if (inherits(object, "metagen")){
        if (!is.null(object$tau.preset))
          meta1 <- metagen(object$TE[sel], object$seTE[sel],
                           sm=object$sm,
                           studlab=object$studlab[sel],
                           level=level, level.comb=level.comb,
                           comb.fixed=comb.fixed,
                           comb.random=comb.random,
                           hakn=object$hakn,
                           method.tau=object$method.tau,
                           tau.preset=object$tau.preset, TE.tau=object$TE.tau,
                           warn=warn)
        else
          meta1 <- metagen(object$TE[sel], object$seTE[sel],
                           sm=object$sm,
                           studlab=object$studlab[sel],
                           level=level, level.comb=level.comb,
                           comb.fixed=comb.fixed,
                           comb.random=comb.random,
                           hakn=object$hakn,
                           method.tau=object$method.tau,
                           TE.tau=object$TE.tau,
                           warn=warn)
      }
      ##
      if (inherits(object, "metaprop")){
        meta1 <- metaprop(object$event[sel], object$n[sel],
                          sm=object$sm,
                          studlab=object$studlab[sel],
                          level=level, level.comb=level.comb,
                          comb.fixed=comb.fixed,
                          comb.random=comb.random,
                          incr=object$incr,
                          allincr=object$allincr,
                          addincr=object$addincr,
                          hakn=object$hakn,
                          method.tau=object$method.tau,
                          tau.preset=object$tau.preset, TE.tau=object$TE.tau,
                          warn=warn)
      }
      ##
      if (inherits(object, "metacor")){
        meta1 <- metacor(object$cor[sel], object$n[sel],
                         sm=object$sm,
                         studlab=object$studlab[sel],
                         level=level, level.comb=level.comb,
                         comb.fixed=comb.fixed,
                         comb.random=comb.random,
                         hakn=object$hakn,
                         method.tau=object$method.tau,
                         tau.preset=object$tau.preset, TE.tau=object$TE.tau)
      }
      ##
      if (inherits(object, "metainc")){
        meta1 <- metainc(object$event.e[sel], object$time.e[sel],
                         object$event.c[sel], object$time.c[sel],
                         studlab=object$studlab[sel],
                         method=object$method,
                         sm=object$sm,
                         incr=object$incr,
                         allincr=object$allincr,
                         addincr=object$addincr,
                         level=level, level.comb=level.comb,
                         comb.fixed=comb.fixed,
                         comb.random=comb.random,
                         hakn=object$hakn,
                         method.tau=object$method.tau,
                         tau.preset=object$tau.preset,
                         TE.tau=object$TE.tau,
                         warn=warn)
      }
      ##
      if (bystud){
        if (print.byvar & bylab!="")
          bylab2 <- paste(bylab, " = ", i, sep="")
        else
          bylab2 <- i
        ##
        lab <- paste(rep("-", nchar(bylab2)), collapse="")
        ##
        cat(lab, "\n", sep="")
        cat(bylab2)
        cat("\n", lab, "\n", sep="")
        ##
        print(meta1, details=FALSE, ma=FALSE)
      }
      ##
      sm1 <- summary(meta1)
      res.w[j,] <- c(meta1$TE.fixed, meta1$seTE.fixed,
                     meta1$Q, meta1$k,
                     meta1$TE.random, meta1$seTE.random,
                     unlist(sm1$H),
                     unlist(sm1$I2),
                     sm1$tau,
                     meta1$C,
                     mean(1/object$n[sel]))
    }
    ##
    TE.fixed.w    <- res.w[,1]
    seTE.fixed.w  <- res.w[,2]
    Q.w           <- res.w[,3]
    k.w           <- res.w[,4]
    TE.random.w   <- res.w[,5]
    seTE.random.w <- res.w[,6]
    ##
    H.w     <- res.w[,7]
    H.w.low <- res.w[,8]
    H.w.upp <- res.w[,9]
    ##
    I2.w     <- res.w[,10]
    I2.w.low <- res.w[,11]
    I2.w.upp <- res.w[,12]
    ##
    tau.w <- res.w[,13]
    ##
    C.w   <- res.w[,14]
    ##
    harmonic.mean.w <- res.w[,15]
    ##
    ci.fixed.w  <- ci(TE.fixed.w, seTE.fixed.w, level.comb)
    ##
    if (!is.null(object$hakn) && object$hakn)
      ci.random.w <- ci(TE.random.w, seTE.random.w, level.comb, df=k.w-1)
    else
      ci.random.w <- ci(TE.random.w, seTE.random.w, level.comb)
    ##
    ci.fixed.w$harmonic.mean <- harmonic.mean.w
    ci.random.w$harmonic.mean <- harmonic.mean.w
    ##
    Q.b <- summary(metagen(TE.fixed.w, seTE.fixed.w))$Q
    ##
    if (bystud) cat("\n")
  }


  ## Calculate prediction interval
  ##
  if (k>=3){
    ci.p <- ci(object$TE.random,
               sqrt(object$seTE.random^2+object$tau^2),
               level.predict, object$k-2)
    ci.p$TE <- NA
  }
  else
    ci.p <- list(TE=NA, seTE=NA,
                 lower=NA, upper=NA, z=NA, p=NA)
  
  
  res <- list(study=ci.study,
              fixed=ci.f, random=ci.r,
              predict=ci.p,
              k=k, Q=Q, df.Q=df.Q,
              tau=object$tau, H=ci.H, I2=ci.I2,
              tau.preset=object$tau.preset,
              k.all=length(object$TE),
              Q.CMH=object$Q.CMH,
              sm=object$sm, method=object$method,
              call=match.call(),
              ci.lab=ci.lab,
              comb.fixed=comb.fixed,
              comb.random=comb.random,
              prediction=prediction)
  
  
  res$se.tau2    <- object$se.tau2
  res$hakn       <- object$hakn
  res$df.hakn    <- object$df.hakn
  res$method.tau <- object$method.tau
  res$TE.tau     <- object$TE.tau
  res$C          <- object$C
  
  
  if (length(byvar)>0){
    res$within.fixed  <- ci.fixed.w
    res$within.random <- ci.random.w
    res$k.w           <- k.w
    res$Q.w           <- Q.w
    res$Q.b.fixed     <- Q.b
    res$Q.b.random    <- summary(metagen(TE.random.w, seTE.random.w))$Q
    res$tau.w         <- tau.w
    res$C.w           <- C.w 
    res$H.w           <- list(TE=H.w, lower=H.w.low, upper=H.w.upp)
    res$I2.w          <- list(TE=I2.w, lower=I2.w.low, upper=I2.w.upp)
    res$bylab         <- bylab
    res$tau.common    <- object$tau.common
    res$by.levs       <- by.levs
    res$within        <- "Returned list 'within' replaced by lists 'within.fixed' and 'within.random'."
  }
  
  
  class(res) <- "summary.meta"
  ##
  if (inherits(object, "metaprop")){
    res$event  <- object$event
    res$n      <- object$n
    res$sparse <- object$sparse
    res$incr <- object$incr
    res$allincr <- object$allincr
    res$addincr <- object$addincr
    ##res$freeman.tukey <- object$freeman.tukey
    ##
    class(res) <- c(class(res), "metaprop")
  }
  ##
  if (inherits(object, "metacor")){
    res$cor       <- object$cor
    res$n         <- object$n
    ##
    class(res) <- c(class(res), "metacor")
  }
  ##
  if (inherits(object, "metabin")){
    class(res) <- c(class(res), "metabin")
    res$sparse <- object$sparse
    res$incr <- object$incr
    res$allincr <- object$allincr
    res$addincr <- object$addincr
    res$MH.exact <- object$MH.exact
  }
  ##
  if (inherits(object, "metainc")){
    class(res) <- c(class(res), "metainc")
    res$sparse <- object$sparse
    res$incr <- object$incr
    res$allincr <- object$allincr
    res$addincr <- object$addincr
  }
  ##
  if (inherits(object, "trimfill")){
    res$object <- object
    res$k0 <- object$k0
    ##
    class(res) <- c(class(res), "trimfill")
  }
  ##
  if (inherits(object, "metacum"))
    class(res) <- c(class(res), "metacum")
  ##
  if (inherits(object, "metainf"))
    class(res) <- c(class(res), "metainf")
  
  res$complab <- object$complab
  res$outclab <- object$outclab
  res$title   <- object$title
  ##
  res$print.byvar <- print.byvar
  
  res$print.CMH <- print.CMH

  res$data <- object$data
  res$subset <- object$subset
  
  res$version <- packageDescription("meta")$Version
  
  res
}
