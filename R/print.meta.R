print.meta <- function(x,
                       sortvar,
                       level=x$level, level.comb=x$level.comb,
                       comb.fixed=x$comb.fixed, comb.random=x$comb.random,
                       prediction=x$prediction, level.predict=x$level.predict,
                       details=FALSE, ma=TRUE,
                       digits=max(4, .Options$digits - 3),
                       ...
                       ){
  
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")
  
  
  rmSpace <- function(x, end=FALSE, pat=" "){
    ##
    if ( !end ){
      while (any(substring(x, 1, 1) == pat, na.rm=TRUE)){
        sel <- substring(x, 1, 1) == pat
        x[sel] <- substring(x[sel], 2)
      }
    }
    else{
      last <- nchar(x)
      ##
      while ( any(substring(x, last, last) == pat, na.rm=TRUE) ){
        sel <- substring(x, last, last) == pat
          x[sel] <- substring(x[sel], 1, last[sel]-1)
        last <- nchar(x)
      }
    }
    ##
    x
  }

  p.ci <- function(lower, upper){
    lower <- rmSpace(lower)
    upper <- rmSpace(upper)
    ##
    ifelse (lower=="NA" & upper=="NA",
            "",
            paste(" [", format(lower, justify="right"),
                  "; ", format(upper, justify="right"), "]", sep=""))
  }

  format.TE <- function(TE, na=FALSE){
    TE <- rmSpace(TE)
    if (na) res <- format(TE)
    else res <- ifelse(is.na(TE), "", format(TE))
    res
  }
  
  
  k.all <- length(x$TE)
  ##
  if (missing(sortvar)) sortvar <- 1:k.all
  ##
  if (length(sortvar) != k.all)
    stop("'x' and 'sortvar' have different length")
  
  
  if (length(comb.fixed)==0){
    comb.fixed <- TRUE
  }
  ##
  if (length(comb.random)==0){
    comb.random <- TRUE
  }
  ##
  if (length(prediction)==0)
    prediction <- FALSE
  ##  
  prediction <- prediction & comb.random & x$k>=3
  
  
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
  ##
  if (length(level.predict)==0){
    if (prediction)
      warning("level.predict set to 0.95")
    level.predict <- 0.95
  }
  
  
  ci.lab <- paste(round(100*level, 1), "%-CI", sep="")

  sm <- x$sm
  
  if (sm=="ZCOR")
    sm.lab <- "COR"
  else if (sm %in% c("PFT", "PAS", "PRAW", "PLN", "PLOGIT"))
    sm.lab <- "proportion"
  else
    sm.lab <- sm
  
  crtitle(x)
  
  
  if (details){

    if (inherits(x, "metabin")){
      res <- cbind(event.e=x$event.e, n.e=x$n.e,
                   event.c=x$event.c, n.c=x$n.c,
                   p.e=round(x$event.e/x$n.e, digits),
                   p.c=round(x$event.c/x$n.c, digits))
    }
    else if (inherits(x, "metacont")){
      res <- cbind(n.e=x$n.e, mean.e=x$mean.e, sd.e=x$sd.e,
                   n.c=x$n.c, mean.c=x$mean.c, sd.c=x$sd.c)
    }
    else if (inherits(x, "metagen")){
      res <- cbind(TE=round(x$TE, digits),
                   seTE=round(x$seTE, digits))
    }
    else if (inherits(x, "metaprop")){
      res <- cbind(event=x$event, n=x$n,
                   p=round(x$event/x$n, digits))
    }
    else if (inherits(x, "metacor")){
      res <- cbind(cor=x$cor, n=x$n)
    }
    else{
      res <- cbind(TE=x$TE, seTE=x$seTE)
    }
    
    dimnames(res)[[1]] <- x$studlab
    
    prmatrix(res[order(sortvar),])
    cat("\n\n")
  }
  
  
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
                            warn=FALSE), level.comb=level.comb),
            digits=digits, header=FALSE)
    else
      print(summary(x, level.comb=level.comb),
            digits=digits, header=FALSE)
  }
  else{
    if (inherits(x, "metainf")|inherits(x, "metacum")){
      TE    <- x$TE
      lowTE <- x$lower
      uppTE <- x$upper
    }
    else{
      tsum <- summary(x, level=level, level.comb=level.comb,
                      level.predict=level.predict, warn=FALSE)
      ##
      TE    <- tsum$study$TE
      lowTE <- tsum$study$lower
      uppTE <- tsum$study$upper
    }
    ##
    if (sm == "RR" | sm == "OR" | sm == "HR"){
      TE    <- exp(TE)
      lowTE <- exp(lowTE)
      uppTE <- exp(uppTE)
    }
    ##
    if (sm=="ZCOR"){
      TE    <- z2cor(TE)
      lowTE <- z2cor(lowTE)
      uppTE <- z2cor(uppTE)
    }
    ##
    if (!inherits(x, "metaprop") & sm=="PLN"){
      TE <- exp(TE)
      lowTE <- exp(lowTE)
      uppTE <- exp(uppTE)
    }
    ##
    if (!inherits(x, "metaprop") & sm=="PLOGIT"){
      TE <- logit2p(TE)
      lowTE <- logit2p(lowTE)
      uppTE <- logit2p(uppTE)
    }
    ##
    if (!inherits(x, "metaprop") & sm %in% c("PFT", "PAS")){
      if (sm=="PAS"){
        TE    <- asin2p(TE, value="mean")
        lowTE <- asin2p(lowTE, value="lower")
        uppTE <- asin2p(uppTE, value="upper")
      }
      ##
      if (sm=="PFT"){
        if (inherits(x, "metainf")|inherits(x, "metacum")){
          TE    <- asin2p(TE, x$n.harmonic.mean, value="mean")
          lowTE <- asin2p(lowTE, x$n.harmonic.mean, value="lower")
          uppTE <- asin2p(uppTE, x$n.harmonic.mean, value="upper")
        }
        else {
          TE    <- asin2p(TE, x$n, value="mean")
          lowTE <- asin2p(lowTE, x$n, value="lower")
          uppTE <- asin2p(uppTE, x$n, value="upper")
        }
      }
    }
    ##
    TE <- round(TE, digits)
    lowTE <- round(lowTE, digits)
    uppTE <- round(uppTE, digits)
    
    if (comb.fixed)
      if (sum(x$w.fixed)>0)
        w.fixed.p <- 100*round(x$w.fixed/sum(x$w.fixed, na.rm=TRUE), 4)
      else w.fixed.p <- x$w.fixed
    
    if (comb.random)
      if (sum(x$w.random)>0)
        w.random.p <- 100*round(x$w.random/sum(x$w.random, na.rm=TRUE), 4)
      else w.random.p <- x$w.random
    
    if (inherits(x, "metainf")|inherits(x, "metacum")){
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
                            c(sm.lab, ci.lab, "p.value", "tau^2", "I^2"))
      ##
      if (inherits(x, "metainf")){
        if (!is.random)
          cat("\nInfluential analysis (Fixed effect model)\n")
        else
          cat("\nInfluential analysis (Random effects model)\n")
      }
      ##
      if (inherits(x, "metacum")){
        if (!is.random)
          cat("\nCumulative meta-analysis (Fixed effect model)\n")
        else
          cat("\nCumulative meta-analysis (Random effects model)\n")
      }
      cat("\n")
      prmatrix(res, quote=FALSE, right=TRUE, na.print="--")
      
      ## Print information on summary method:
      ##
      catmeth(method=x$method,
              method.tau=if (is.random) x$method.tau else "",
              sm=sm,
              k.all=k.all,
              hakn=is.random & x$hakn,
              metaprop=inherits(x, "metaprop"),
              trimfill=inherits(x, "trimfill"))
    }
    else{
      res <- cbind(format.TE(TE, na=TRUE),
                   p.ci(format(lowTE), format(uppTE)),
                   if (comb.fixed) format(w.fixed.p),
                   if (comb.random) format(w.random.p))
      ##
      ## Printout for a single proportion:
      ##
      if (k.all==1){
        cat("Exact CI:\n\n")
        dimnames(res) <-
          list("", c(sm.lab, ci.lab,
                     if (comb.fixed) "%W(fixed)",
                     if (comb.random) "%W(random)"))
        prmatrix(res, quote=FALSE, right=TRUE)
        cat("\n\nCI based on arcsine transformation:\n")
      }
      else{
        dimnames(res) <-
          list(x$studlab, c(sm.lab, ci.lab,
                            if (comb.fixed) "%W(fixed)",
                            if (comb.random) "%W(random)"))
        prmatrix(res[order(sortvar),], quote=FALSE, right=TRUE)
      }
      cat("\n")
    }
    
    
    if (ma&!(inherits(x, "metainf")|inherits(x, "metacum")))
      print(tsum, digits=digits,
            comb.fixed=comb.fixed, comb.random=comb.random,
            prediction=prediction,
            header=FALSE)
    
    }
  
  invisible(NULL)
}
