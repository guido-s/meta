summary.meta <- function(object,
                         comb.fixed=object$comb.fixed,
                         comb.random=object$comb.random,
                         prediction=object$prediction,
                         backtransf=object$backtransf,
                         bylab=object$bylab,
                         print.byvar=object$print.byvar,
                         bystud=FALSE,
                         print.CMH=object$print.CMH,
                         warn=object$warn,
                         ...){
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(object, "meta")
  ##
  if (inherits(object, "metacum")){
    warning("Summary method not defined for objects of class \"metacum\".")
    return(object)
  }
  ##
  if (inherits(object, "metainf")){
    warning("Summary method not defined for objects of class \"metainf\".")
    return(object)
  }
  ##
  if (length(warn)==0)
    warn <- .settings$warn
  object <- updateversion(object)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(prediction)
  chklogical(backtransf)
  if (!is.null(print.byvar))
    chklogical(print.byvar)
  chklogical(bystud)
  if (!is.null(print.CMH))
    chklogical(print.CMH)
  chklogical(warn)
  ##  
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
  ##
  ## (3) Results for individual studies
  ##
  ##
  ci.study <- list(TE=object$TE,
                   seTE=object$seTE,
                   lower=object$lower,
                   upper=object$upper,
                   z=object$zval,
                   p=object$pval,
                   level=object$level,
                   df=NA)
  ##
  if (inherits(object, "metaprop")){
    ci.study$event <- object$event
    ci.study$n <- object$n
  }
  
  
  ##
  ##
  ## (4) Results for meta-analysis
  ##
  ##
  ci.f <- list(TE=object$TE.fixed,
               seTE=object$seTE.fixed,
               lower=object$lower.fixed,
               upper=object$upper.fixed,
               z=object$zval.fixed,
               p=object$pval.fixed,
               level=object$level.comb)
  if (inherits(object, "metaprop"))
    ci.f$harmonic.mean <- mean(1/object$n)
  ##
  ci.r <- list(TE=object$TE.random,
               seTE=object$seTE.random,
               lower=object$lower.random,
               upper=object$upper.random,
               z=object$zval.random,
               p=object$pval.random,
               level=object$level.comb,
               df=if (!is.null(object$df.hakn)) object$df.hakn else NA)
  if (inherits(object, "metaprop"))
    ci.r$harmonic.mean <- mean(1/object$n)
  ##
  ci.H <- list(TE=object$H, lower=object$lower.H, upper=object$upper.H)
  ##
  ci.I2 <- list(TE=object$I2, lower=object$lower.I2, upper=object$upper.I2)
  ##
  ci.p <- list(TE=NA,
               seTE=object$seTE.predict,
               lower=object$lower.predict,
               upper=object$upper.predict,
               z=NA,
               p=NA,
               level=object$level.predict,
               df=object$k-2)
  ##  
  ci.lab <- paste(round(100*object$level.comb, 1), "%-CI", sep="")
  
  
  ##
  ##
  ## (5) Generate R object
  ##
  ##
  res <- list(study=ci.study,
              fixed=ci.f, random=ci.r,
              predict=ci.p,
              k=object$k, Q=object$Q, df.Q=object$df.Q,
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
  ##  
  res$se.tau2    <- object$se.tau2
  res$hakn       <- object$hakn
  res$df.hakn    <- object$df.hakn
  res$method.tau <- object$method.tau
  res$TE.tau     <- object$TE.tau
  res$C          <- object$C
  ##
  ## Add results from subgroup analysis
  ##
  if (length(object$byvar)>0){
    ##
    ci.fixed.w <- list(TE=object$TE.fixed.w,
                       seTE=object$seTE.fixed.w,
                       lower=object$lower.fixed.w,
                       upper=object$upper.fixed.w,
                       z=object$zval.fixed.w,
                       p=object$pval.fixed.w,
                       level=object$level.comb,
                       harmonic.mean=object$n.harmonic.mean.w)
    ##
    ci.random.w <- list(TE=object$TE.random.w,
                       seTE=object$seTE.random.w,
                       lower=object$lower.random.w,
                       upper=object$upper.random.w,
                       z=object$zval.random.w,
                       p=object$pval.random.w,
                       level=object$level.comb,
                       df=object$df.hakn.w,
                       harmonic.mean=object$n.harmonic.mean.w)
    ##
    ci.H <- list(TE=object$H.w, lower=object$lower.H.w, upper=object$upper.H.w)
    ci.I2 <- list(TE=object$I2.w, lower=object$lower.I2.w, upper=object$upper.I2.w)
    ## 
    res$within.fixed  <- ci.fixed.w
    res$within.random <- ci.random.w
    res$k.w           <- object$k.w
    res$Q.w           <- object$Q.w
    res$Q.w.fixed     <- object$Q.w.fixed
    res$Q.w.random    <- object$Q.w.random
    res$df.Q.w        <- object$df.Q.w
    res$Q.b.fixed     <- object$Q.b.fixed
    res$Q.b.random    <- object$Q.b.random
    res$df.Q.b        <- object$df.Q.b
    res$tau.w         <- object$tau.w
    res$C.w           <- object$C.w
    res$H.w           <- ci.H
    res$I2.w          <- ci.I2
    res$bylab         <- object$bylab
    res$tau.common    <- object$tau.common
    res$bylevs        <- object$bylevs
    res$within        <- "Returned list 'within' replaced by lists 'within.fixed' and 'within.random'."
  }
  ##  
  class(res) <- "summary.meta"
  ##
  if (inherits(object, "metabin")){
    res$sparse   <- object$sparse
    res$incr     <- object$incr
    res$allincr  <- object$allincr
    res$addincr  <- object$addincr
    res$MH.exact <- object$MH.exact
    ##
    class(res) <- c(class(res), "metabin")
  }
  ##
  if (inherits(object, "metacont")){
    res$pooledvar <- object$pooledvar
    res$method.smd <- object$method.smd
    res$sd.glass <- object$sd.glass
    res$exact.smd <- object$exact.smd
    ##
    class(res) <- c(class(res), "metacont")
  }
  ##
  if (inherits(object, "metacor")){
    res$cor <- object$cor
    res$n   <- object$n
    ##
    class(res) <- c(class(res), "metacor")
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
  if (inherits(object, "metaprop")){
    res$event  <- object$event
    res$n      <- object$n
    res$sparse <- object$sparse
    res$incr <- object$incr
    res$allincr <- object$allincr
    res$addincr <- object$addincr
    res$method.ci <- object$method.ci
    ##
    class(res) <- c(class(res), "metaprop")
  }
  ##
  if (inherits(object, "trimfill")){
    res$object <- object
    res$k0 <- object$k0
    ##
    class(res) <- c(class(res), "trimfill")
  }
  ##
  res$complab <- object$complab
  res$outclab <- object$outclab
  res$title   <- object$title
  ##
  res$print.byvar <- print.byvar
  res$print.CMH <- print.CMH
  ##
  res$data <- object$data
  res$subset <- object$subset
  ##
  res$backtransf <- backtransf
  ##
  res$version <- packageDescription("meta")$Version
  
  
  res
}
