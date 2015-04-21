print.summary.meta <- function(x,
                               digits=max(3, .Options$digits - 3),
                               comb.fixed=x$comb.fixed,
                               comb.random=x$comb.random,
                               prediction=x$prediction,
                               print.byvar=x$print.byvar,
                               print.CMH=x$print.CMH,
                               header=TRUE,
                               backtransf=x$backtransf,
                               bylab.nchar=35,
                               ...){
  
  
  ##
  ##
  ## (1) Check for summary.meta object
  ##
  ##
  chkclass(x, "summary.meta")
  ##
  if (inherits(x, "metacum") | inherits(x, "metainf"))
    return(invisible(NULL))
  ##
  by <- !is.null(x$bylab)
  
  
  ##
  ##
  ## (2) Check other arguments
  ##
  ##
  chknumeric(digits)
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(prediction)
  if (by)
    chklogical(print.byvar)
  if (!is.null(print.CMH))
    chklogical(print.CMH)
  chklogical(header)
  chklogical(backtransf)
  chknumeric(bylab.nchar)
  ##
  ## Additional arguments / checks for metacont objects
  ##
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "print.summary.meta"
  ##
  warnarg("logscale", addargs, fun, otherarg="backtransf")
  
  
  ##
  ##
  ## (3) Some additional settings
  ##
  ##
  k <- x$k
  sm <- x$sm
  ##
  bip <- inherits(x, c("metabin", "metainc", "metaprop"))
  ##
  prediction <- prediction & k>=3
  ##
  if (is.null(x$df.Q))
    df.Q <- k-1
  else
    df.Q <- x$df.Q
  ##
  if (by){
    k.w <- x$k.w
    if (is.null(x$df.Q.w))
      df.Q.w <- sum((k.w-1)[!is.na(x$Q.w)])
    else
      df.Q.w <- x$df.Q.w
    ##
    if (is.null(x$df.Q.b))
      df.Q.b <- (k-1) - sum((k.w-1)[!is.na(x$Q.w)])
    else
      df.Q.b <- x$df.Q.b
  }
  ##  
  sm.lab <- sm
  ##
  if (backtransf){
    if (sm=="ZCOR")
      sm.lab <- "COR"
    if (sm %in% c("PFT", "PAS", "PRAW", "PLOGIT", "PLN"))
      sm.lab <- "proportion"
  }
  else 
    if (is.relative.effect(sm))
      sm.lab <- paste("log", sm, sep="")
  ##
  if (length(x$tau.common)==0)
    x$tau.common <- FALSE
  ##
  if (length(x$tau.common)==0)
    x$tau.common <- FALSE
  ##
  if (by)
    bylevs <- ifelse(nchar(x$bylevs) > bylab.nchar,
                     paste(substring(x$bylevs, 1, bylab.nchar-4), " ...", sep=""),#
                     x$bylevs)
  
  
  ##
  ##
  ## (4) Set and backtransform results of meta-analysis
  ##
  ##
  TE.fixed    <- x$fixed$TE
  lowTE.fixed <- x$fixed$lower
  uppTE.fixed <- x$fixed$upper
  ##
  TE.random    <- x$random$TE
  lowTE.random <- x$random$lower
  uppTE.random <- x$random$upper
  ##
  lowTE.predict <- x$predict$lower
  uppTE.predict <- x$predict$upper
  ##
  if (by){
    TE.fixed.w     <- x$within.fixed$TE
    lowTE.fixed.w  <- x$within.fixed$lower
    uppTE.fixed.w  <- x$within.fixed$upper
    pval.fixed.w   <- x$within.fixed$p
    n.harmonic.mean.w <- x$within.fixed$harmonic.mean
    TE.random.w    <- x$within.random$TE
    lowTE.random.w <- x$within.random$lower
    uppTE.random.w <- x$within.random$upper
    pval.random.w   <- x$within.random$p
  }
  ##  
  if (backtransf){
    npft.ma <- 1/mean(1/x$n)
    ##
    TE.fixed    <- backtransf(TE.fixed, sm, "mean",
                              npft.ma, warn=comb.fixed)
    lowTE.fixed <- backtransf(lowTE.fixed, sm, "lower",
                              npft.ma, warn=comb.fixed)
    uppTE.fixed <- backtransf(uppTE.fixed, sm, "upper",
                              npft.ma, warn=comb.fixed)
    ##
    TE.random <- backtransf(TE.random, sm, "mean",
                            npft.ma, warn=comb.random)
    lowTE.random <- backtransf(lowTE.random, sm, "lower",
                               npft.ma, warn=comb.random)
    uppTE.random <- backtransf(uppTE.random, sm, "upper",
                               npft.ma, warn=comb.random)
    ##
    lowTE.predict <- backtransf(lowTE.predict, sm, "lower",
                                npft.ma, warn=prediction)
    uppTE.predict <- backtransf(uppTE.predict, sm, "upper",
                                npft.ma, warn=prediction)
    ##
    if (by){
      npft.w <- n.harmonic.mean.w
      ##
      TE.fixed.w     <- backtransf(TE.fixed.w, sm, "mean",
                                   npft.w, warn=comb.fixed)
      lowTE.fixed.w  <- backtransf(lowTE.fixed.w, sm, "lower",
                                   npft.w, warn=comb.fixed)
      uppTE.fixed.w  <- backtransf(uppTE.fixed.w, sm, "upper",
                                   npft.w, warn=comb.fixed)
      ##
      TE.random.w    <- backtransf(TE.random.w, sm, "mean",
                                   npft.w, warn=comb.random)
      lowTE.random.w <- backtransf(lowTE.random.w, sm, "lower",
                                   npft.w, warn=comb.random)
      uppTE.random.w <- backtransf(uppTE.random.w, sm, "upper",
                                   npft.w, warn=comb.random)
    }
  }
  ##  
  TE.fixed    <- round(TE.fixed, digits)
  lowTE.fixed <- round(lowTE.fixed, digits)
  uppTE.fixed <- round(uppTE.fixed, digits)
  pTE.fixed <- x$fixed$p
  zTE.fixed <- round(x$fixed$z, digits)
  ##
  TE.random    <- round(TE.random, digits)
  lowTE.random <- round(lowTE.random, digits)
  uppTE.random <- round(uppTE.random, digits)
  pTE.random <- x$random$p
  zTE.random <- round(x$random$z, digits)
  ##
  lowTE.predict <- round(lowTE.predict, digits)
  uppTE.predict <- round(uppTE.predict, digits)
  ##
  if (by){
    TE.fixed.w     <- round(TE.fixed.w, digits)
    lowTE.fixed.w  <- round(lowTE.fixed.w, digits)
    uppTE.fixed.w  <- round(uppTE.fixed.w, digits)
    TE.random.w    <- round(TE.random.w, digits)
    lowTE.random.w <- round(lowTE.random.w, digits)
    uppTE.random.w <- round(uppTE.random.w, digits)
  }
  ##
  H <- x$H$TE
  lowH <- x$H$lower
  uppH <- x$H$upper
  ##
  I2 <- x$I2$TE
  lowI2 <- x$I2$lower
  uppI2 <- x$I2$upper
  
  
  ##
  ##
  ## (5) Print result for meta-analysis
  ##
  ##
  if (header)
    crtitle(x)
  ##
  if (x$k.all == 1){
    ##
    ## Print results for a single study
    ##
    res <- cbind(TE.fixed,
                 p.ci(format(lowTE.fixed), format(uppTE.fixed)),
                 format(round(zTE.fixed,4)),
                 format.p(pTE.fixed))
    dimnames(res) <- list("", c(sm.lab, x$ci.lab, "z", "p-value"))
    prmatrix(res, quote=FALSE, right=TRUE, ...)
    ## Print information on summary method:
    catmeth(method=x$method,
            sm=sm,
            k.all=x$k.all,
            metaprop=inherits(x, "metaprop"),
            metabin=inherits(x, "metabin"),
            metainc=inherits(x, "metainc"),
            sparse=ifelse(bip, x$sparse, FALSE),
            incr=ifelse(bip, x$incr, FALSE),
            allincr=ifelse(bip, x$allincr, FALSE),
            addincr=ifelse(bip, x$addincr, FALSE),
            method.ci=x$method.ci,
            metacont=inherits(x, "metacont"),
            pooledvar=x$pooledvar,
            method.smd=x$method.smd,
            sd.glass=x$sd.glass,
            exact.smd=x$exact.smd)
  }
  else{
    ##
    ## Print results for meta-analysis with more than one study
    ##
    if (comb.fixed|comb.random|prediction){
      if (!inherits(x, "trimfill"))
        cat(paste("Number of studies combined: k=", k, "\n\n", sep=""))
      else
        cat(paste("Number of studies combined: k=", k,
                  " (with ", x$k0, " added studies)\n\n", sep=""))
      res <- cbind(format(c(if (comb.fixed) TE.fixed,
                            if (comb.random) TE.random,
                            if (prediction) NA)),
                   p.ci(format(c(if (comb.fixed) lowTE.fixed,
                                 if (comb.random) lowTE.random,
                                 if (prediction) lowTE.predict)),
                        format(c(if (comb.fixed) uppTE.fixed,
                                 if (comb.random) uppTE.random,
                                 if (prediction) uppTE.predict))),
                   format(round(c(if (comb.fixed) zTE.fixed,
                                  if (comb.random) zTE.random,
                                  if (prediction) NA),4)),
                   format.p(c(if (comb.fixed) pTE.fixed,
                              if (comb.random) pTE.random,
                              if (prediction) NA)))
      if (prediction)
        res[dim(res)[1], c(1,3:4)] <- ""
      if (!is.null(x$hakn) && x$hakn){
        if (comb.fixed & comb.random)
          zlab <- "z|t"
        else if (comb.fixed & !comb.random)
          zlab <- "z"
        else if (!comb.fixed & comb.random)
          zlab <- "t"
      }
      else
        zlab <- "z"
      ##
      dimnames(res) <- list(c(if (comb.fixed) "Fixed effect model",
                              if (comb.random) "Random effects model",
                              if (prediction) "Prediction interval"),  
                            c(sm.lab, x$ci.lab, zlab, "p-value"))
      prmatrix(res, quote=FALSE, right=TRUE, ...)
      ##
      if (inherits(x, "metabin") && print.CMH){
        Qdata <- cbind(round(x$Q.CMH, 2), 1,
                       format.p(1-pchisq(x$Q.CMH, df=1)))
        dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
        ##
        cat("\nCochran-Mantel-Haenszel (CMH) test for overall effect: \n")
        prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
      }
    }
    else
      cat(paste("Number of studies: k=", k, "\n", sep=""))
    ##
    ## Print information on heterogeneity
    ##
    if (!is.na(x$tau))
      cat(paste("\nQuantifying heterogeneity:\n",
                if (x$tau^2 > 0 & x$tau^2 < 0.0001)
                paste("tau^2", format.tau(x$tau^2))
                else
                paste("tau^2 = ",
                      ifelse(x$tau==0,
                             "0",
                             format(round(x$tau^2, 4), 4, nsmall=4, scientific=FALSE)),
                      sep=""),
                paste("; H = ", round(H, 2),
                      ifelse(k>2,
                             paste(" ", p.ci(round(lowH, 2), round(uppH, 2)), sep=""),
                             ""),
                      "; ",
                      "I^2 = ", round(100*I2, 1), "%",
                      ifelse(k>2,
                             paste(" ",
                                   p.ci(paste(round(100*lowI2, 1), "%", sep=""),
                                        paste(round(100*uppI2, 1), "%", sep="")),
                                   sep=""),
                             ""),
                      sep=""),
                "\n", sep=""))
    ##    
    if (k > 1 & (comb.fixed|comb.random)){
      Qdata <- cbind(round(x$Q, 2), df.Q,
                     format.p(1-pchisq(x$Q, df=df.Q)))
      dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
      ##
      cat("\nTest of heterogeneity:\n")
      prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
      ##
      if (by){
        ##
        ## Print information for subgroup analysis
        ##
        if (comb.fixed){
          ##
          ## Subgroup analysis based on fixed effect model
          ##
          Q.w <- ifelse(is.na(x$Q.w),
                        "--",
                        format(round(x$Q.w, 2)))
          I2.w <- ifelse(is.na(x$I2.w$TE),
                         "--",
                         paste(round(100*x$I2.w$TE, 1), "%", sep=""))
          tau2.w <- ifelse(k.w==1, "--", format.tau(x$tau.w^2))
          ##
          Tdata <- cbind(format(k.w),
                         format(TE.fixed.w),
                         c(p.ci(format(lowTE.fixed.w), format(uppTE.fixed.w))),
                         Q.w, tau2.w, I2.w
                         ) #, format.p(pval.fixed.w))
          if (print.byvar)
            bylab <- paste(x$bylab,
                           " = ", 
                           format(bylevs), sep="")
          else
            bylab <- format(bylevs)
          dimnames(Tdata) <- list(bylab,
                                  c("  k", sm.lab, x$ci.lab,
                                    "Q", "tau^2", "I^2")
                                  ) #, "p-value"))
          cat("\nResults for subgroups (fixed effect model):\n")
          prmatrix(Tdata, quote=FALSE, right=TRUE, ...)
          ##
          cat("\nTest for subgroup differences (fixed effect model):\n")
          if (x$method=="MH"){
            Qdata <- cbind(round(x$Q.b.fixed, 2), df.Q.b,
                           format.p(1-pchisq(x$Q.b.fixed, df=df.Q.b)))
            dimnames(Qdata) <- list("Between groups  ",
                                    c("Q", "d.f.", "p-value"))
            prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
          }
          else{
            Qs  <- c(x$Q.b.fixed, x$Q.w.fixed)
            dfs <- c(df.Q.b, df.Q.w)
            Qdata <- cbind(round(Qs, 2),
                           dfs,
                           format.p(1-pchisq(Qs, df=dfs)))
            dimnames(Qdata) <- list(c("Between groups", "Within groups"),
                                    c("Q", "d.f.", "p-value"))
            prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
          }
        }
        ##
        if (comb.random){
          ##
          ## Subgroup analysis based on random effects model
          ##
          Q.w <- ifelse(is.na(x$Q.w),
                        "--",
                        format(round(x$Q.w, 2)))
          I2.w <- ifelse(is.na(x$I2.w$TE),
                         "--",
                         paste(round(100*x$I2.w$TE, 1), "%", sep=""))
          tau2.w <- ifelse(k.w==1, "--", format.tau(x$tau.w^2))
          ##
          Tdata <- cbind(format(k.w),
                         format(TE.random.w),
                         c(p.ci(format(lowTE.random.w), format(uppTE.random.w))),
                         Q.w, tau2.w, I2.w
                         )
          if (print.byvar)
            bylab <- paste(x$bylab,
                           " = ", 
                           format(bylevs), sep="")
          else
            bylab <- format(bylevs)
          dimnames(Tdata) <- list(bylab,
                                  c("  k", sm.lab, x$ci.lab,
                                    "Q", "tau^2", "I^2")
                                  ) #, "p-value"))
          cat("\nResults for subgroups (random effects model):\n")
          prmatrix(Tdata, quote=FALSE, right=TRUE, ...)
          ##
          cat("\nTest for subgroup differences (random effects model):\n")
          if (is.na(x$Q.w.random)){
            Qdata <- cbind(round(x$Q.b.random, 2), df.Q.b,
                           format.p(1-pchisq(x$Q.b.random, df=df.Q.b)))
          dimnames(Qdata) <- list("Between groups  ",
                                  c("Q", "d.f.", "p-value"))
          }
          else{
            Qs  <- c(x$Q.b.random, x$Q.w.random)
            dfs <- c(df.Q.b, df.Q.w)
            Qdata <- cbind(round(Qs, 2),
                           dfs, format.p(1-pchisq(Qs, df=dfs)))
            dimnames(Qdata) <- list(c("Between groups", "Within groups"),
                                    c("Q", "d.f.", "p-value"))
          }
          prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
        }
      }
    }
    ##
    ## Print information on summary method:
    ##
    catmeth(method=x$method,
            method.tau=if (comb.random) x$method.tau else "",
            sm=sm,
            k.all=x$k.all,
            hakn=!is.null(x$hakn) && (x$hakn & comb.random),
            tau.common=by & x$tau.common,
            tau.preset=x$tau.preset,
            trimfill=inherits(x, "trimfill"),
            metaprop=inherits(x, "metaprop"),
            metabin=inherits(x, "metabin"),
            metainc=inherits(x, "metainc"),
            sparse=ifelse(bip, x$sparse, FALSE),
            incr=ifelse(bip, x$incr, FALSE),
            allincr=ifelse(bip, x$allincr, FALSE),
            addincr=ifelse(bip, x$addincr, FALSE),
            MH.exact=ifelse(inherits(x, "metabin"), x$MH.exact, FALSE),
            method.ci=x$method.ci,
            metacont=inherits(x, "metacont"),
            pooledvar=x$pooledvar,
            method.smd=x$method.smd,
            sd.glass=x$sd.glass,
            exact.smd=x$exact.smd)
  }
  
  
  invisible(NULL)
}
