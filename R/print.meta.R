print.meta <- function(x,
                       sortvar,
                       comb.fixed=x$comb.fixed,
                       comb.random=x$comb.random,
                       prediction=x$prediction,
                       details=FALSE, ma=TRUE,
                       backtransf=x$backtransf,
                       digits=max(4, .Options$digits - 3),
                       ...
                       ){
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta object
  ##
  ##
  chkclass(x, "meta")
  x <- updateversion(x)
  ##
  k.all <- length(x$TE)


  ##
  ##
  ## (2) Check other arguments
  ##
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
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(prediction)
  chklogical(details)
  chklogical(ma)
  chklogical(backtransf)
  chknumeric(digits)
  ##
  ## Additional arguments / checks for metacont objects
  ##
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "print.meta"
  ##
  warnarg("level", addargs, fun, cl)
  warnarg("level.comb", addargs, fun, cl)
  warnarg("level.predict", addargs, fun, cl)
  warnarg("logscale", addargs, fun, otherarg="backtransf")
  ##
  level <- x$level
  level.comb <- x$level.comb
  level.predict <- x$level.predict
  ##  
  format.TE <- function(TE, na=FALSE){
    TE <- rmSpace(TE)
    if (na) res <- format(TE)
    else res <- ifelse(is.na(TE), "", format(TE))
    res
  }
  
  
  ##
  ##
  ## (3) Some additional settings
  ##
  ##
  metainf.metacum <- inherits(x, "metainf") | inherits(x, "metacum")
  ##  
  prediction <- prediction & x$k>=3
  ##  
  ci.lab <- paste(round(100*level, 1), "%-CI", sep="")
  ##  
  sm <- x$sm
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
  ##
  ## (4) Print title and details
  ##
  ##
  crtitle(x)
  ##  
  if (details){
    if (inherits(x, "metabin")){
      res <- data.frame(event.e=x$event.e, n.e=x$n.e,
                        event.c=x$event.c, n.c=x$n.c,
                        p.e=round(x$event.e/x$n.e, digits),
                        p.c=round(x$event.c/x$n.c, digits))
    }
    else if (inherits(x, "metacont")){
      res <- data.frame(n.e=x$n.e, mean.e=x$mean.e, sd.e=x$sd.e,
                        n.c=x$n.c, mean.c=x$mean.c, sd.c=x$sd.c)
    }
    else if (inherits(x, "metacor")){
      res <- data.frame(cor=x$cor, n=x$n)
    }
    else if (inherits(x, "metagen")){
      res <- data.frame(TE=round(x$TE, digits),
                        seTE=round(x$seTE, digits))
    }
    else if (inherits(x, "metainc")){
      res <- data.frame(event.e=x$event.e, time.e=x$time.e,
                        event.c=x$event.c, time.c=x$time.c)
    }
    else if (inherits(x, "metaprop")){
      res <- data.frame(event=x$event, n=x$n,
                        p=round(x$event/x$n, digits))
    }
    else{
      res <- data.frame(TE=x$TE, seTE=x$seTE)
    }
    dimnames(res)[[1]] <- x$studlab
    prmatrix(res[order(sortvar),])
    cat("\n\n")
  }
  
  
  ##
  ##
  ## (5) Print results for individual studies
  ##
  ##
  if (k.all == 1 & !inherits(x, "metaprop")){
    if (inherits(x, "metabin") & x$method=="MH")
      print(summary(metabin(x$event.e, x$n.e,
                            x$event.c, x$n.c,
                            sm=sm,
                            method="Inverse",
                            studlab=x$studlab,
                            incr=x$incr,
                            allincr=x$allincr,
                            allstudies=x$allstudies,
                            MH.exact=x$MH.exact,
                            warn=FALSE, level.comb=level.comb)),
            digits=digits, header=FALSE, backtransf=backtransf)
    else
      print(summary(x),
            digits=digits, header=FALSE, backtransf=backtransf)
  }
  else{
    TE <- x$TE
    seTE <- x$seTE
    lowTE <- x$lower
    uppTE <- x$upper
    ##
    if (inherits(x, "metaprop") & !backtransf){
      ciTE <- ci(TE, seTE, level=level)
      lowTE <- ciTE$lower
      uppTE <- ciTE$upper
      ##
      x$method.ci <- "NAsm"
    }
    ##
    if (backtransf){
      ## Freeman-Tukey Arcsin transformation (sm="PFT")
      if (metainf.metacum)
        npft <- x$n.harmonic.mean
      else
        npft <- x$n
      ##
      if (inherits(x, "metaprop"))
        TE <- x$event/x$n
      else{
        TE <- backtransf(TE, sm, "mean", npft)
        lowTE <- backtransf(lowTE, sm, "lower", npft)
        uppTE <- backtransf(uppTE, sm, "upper", npft)
      }
    }
    ##
    TE <- round(TE, digits)
    lowTE <- round(lowTE, digits)
    uppTE <- round(uppTE, digits)
    ##    
    if (comb.fixed)
      if (sum(x$w.fixed)>0)
        w.fixed.p <- 100*round(x$w.fixed/sum(x$w.fixed, na.rm=TRUE), 4)
      else w.fixed.p <- x$w.fixed
    ##    
    if (comb.random)
      if (sum(x$w.random)>0)
        w.random.p <- 100*round(x$w.random/sum(x$w.random, na.rm=TRUE), 4)
      else w.random.p <- x$w.random
    ##    
    if (metainf.metacum){
      is.random <- x$pooled=="random"
      ##
      sel1 <- is.na(x$I2)
      I2 <- 100*x$I2
      I2 <- ifelse(sel1, "", format(round(I2, 1)))
      ##
      sel2 <- is.na(x$p.value)
      p.value <- format.p(x$p.value)
      p.value <- ifelse(sel2, "", p.value)
      ##
      sel3 <- is.na(x$tau)
      tau2 <- x$tau^2
      tau2 <- ifelse(sel3, "", round(tau2, digits))
      ##
      res <- cbind(format.TE(TE), p.ci(format(lowTE), format(uppTE)),
                   p.value, paste("  ", format(tau2), sep=""),
                   paste("  ", I2, ifelse(sel1, "", "%"), sep=""))
      dimnames(res) <- list(paste(x$studlab, "  ", sep=""),
                            c(sm.lab, ci.lab, "p-value", "tau^2", "I^2"))
      ##
      if (inherits(x, "metainf")){
        if (!is.random)
          cat("\nInfluential analysis (Fixed effect model)\n")
        else
          cat("\nInfluential analysis (Random effects model)\n")
      }
      if (inherits(x, "metacum")){
        if (!is.random)
          cat("\nCumulative meta-analysis (Fixed effect model)\n")
        else
          cat("\nCumulative meta-analysis (Random effects model)\n")
      }
      cat("\n")
      prmatrix(res, quote=FALSE, right=TRUE, na.print="--")
      ## Print information on summary method:
      catmeth(method=x$method,
              method.tau=if (is.random) x$method.tau else "",
              sm=sm,
              k.all=k.all,
              hakn=is.random & x$hakn,
              metaprop=inherits(x, "metaprop"),
              trimfill=inherits(x, "trimfill"),
              tau.preset=x$tau.preset,
              method.smd=x$method.smd,
              sd.glass=x$sd.glass,
              exact.smd=x$exact.smd)
    }
    else{
      res <- cbind(format.TE(TE, na=TRUE),
                   p.ci(format(lowTE), format(uppTE)),
                   if (comb.fixed) format(w.fixed.p),
                   if (comb.random) format(w.random.p))
      ## Printout for a single proportion:
      if (k.all == 1){
        ##
        if (!is.null(x$method.ci)){
          if  (x$method.ci=="CP")
            method.ci.details <- "Clopper-Pearson confidence interval:\n\n"
          else if (x$method.ci=="WS")
            method.ci.details <- "Wilson Score confidence interval:\n\n"
          else if (x$method.ci=="WSCC")
            method.ci.details <- "Wilson Score confidence interval with continuity correction:\n\n"
          else if (x$method.ci=="AC")
            method.ci.details <- "Agresti-Coull confidence interval:\n\n"
          else if (x$method.ci=="SA")
            method.ci.details <- "Simple approximation confidence interval:\n\n"
          else if (x$method.ci=="SACC")
            method.ci.details <- "Simple approximation confidence interval with continuity correction:\n\n"
          if (x$method.ci!="NAsm"){
            cat(method.ci.details)
            dimnames(res) <- list("", c(sm.lab, ci.lab,
                                        if (comb.fixed) "%W(fixed)",
                                        if (comb.random) "%W(random)"))
            prmatrix(res, quote=FALSE, right=TRUE)
            cat("\n\n")
          }
        }
        cat("Normal approximation confidence interval:\n")
      }
      else{
        dimnames(res) <-
          list(x$studlab, c(sm.lab, ci.lab,
                            if (comb.fixed) "%W(fixed)",
                            if (comb.random) "%W(random)"))
        prmatrix(res[order(sortvar),], quote=FALSE, right=TRUE)
      }
    }
    
    
    ##
    ##
    ## (6) Print result for meta-analysis
    ##
    ##
    if (ma & !metainf.metacum){
      cat("\n")
      print(summary(x, warn=FALSE), digits=digits,
            comb.fixed=comb.fixed, comb.random=comb.random,
            prediction=prediction,
            header=FALSE, backtransf=backtransf)
    }
  }
  
  
  invisible(NULL)
}
