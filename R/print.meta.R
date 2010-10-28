print.meta <- function(x,
                       sortvar,
                       level=x$level, level.comb=x$level.comb,
                       comb.fixed=x$comb.fixed, comb.random=x$comb.random,
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
  
  
  ci.lab <- paste(round(100*level, 1), "%-CI", sep="")
  
  
  if (x$sm=="ZCOR")
    sm.lab <- "COR"
  else if (x$sm %in% c("PFT", "PAS", "PRAW", "PLN", "PLOGIT"))
    sm.lab <- "proportion"
  else
    sm.lab <- x$sm
  
  
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
  
  
  tl <- options()$width-12
  ##
  if (!is.null(x$title))
    if (x$title!="")
      if (nchar(x$title) <= tl)
        cat("Review:     ", x$title, "\n", sep="")
      else
        cat("Review:     ", substring(x$title, 1, tl-4),
            " ...\n", sep="")
  if (!is.null(x$complab))
    if (x$complab!="")
      if (nchar(x$complab) <= tl)
        cat("Comparison: ", x$complab, "\n", sep="")
      else
        cat("Comparison: ", substring(x$complab, 1, tl-4),
            " ...\n", sep="")
  if (!is.null(x$outclab))
    if (x$outclab!="")
      if (nchar(x$outclab) <= tl)
        cat("Outcome:    ", x$outclab, "\n\n", sep="")
      else
        cat("Outcome:    ", substring(x$outclab, 1, tl-4),
            " ...\n\n", sep="")
  
  
  if (k.all == 1 & !inherits(x, "metaprop")){
    if (inherits(x, "metabin") & x$method=="MH")
      print(summary(metabin(x$event.e, x$n.e,
                            x$event.c, x$n.c,
                            sm=x$sm,
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
    tsum <- summary(x, level=level, level.comb=level.comb, warn=FALSE)
    ##
    TE    <- tsum$study$TE
    lowTE <- tsum$study$lower
    uppTE <- tsum$study$upper
    ##
    if (x$sm == "RR" | x$sm == "OR" | x$sm == "HR"){
      TE    <- exp(TE)
      lowTE <- exp(lowTE)
      uppTE <- exp(uppTE)
    }
    ##
    if (x$sm=="ZCOR"){
      TE    <- z2cor(TE)
      lowTE <- z2cor(lowTE)
      uppTE <- z2cor(uppTE)
    }
    ##
    if (!inherits(x, "metaprop") & x$sm=="PLN"){
      TE <- exp(TE)
      lowTE <- exp(lowTE)
      uppTE <- exp(uppTE)
    }
    ##
    if (!inherits(x, "metaprop") & x$sm=="PLOGIT"){
      TE <- meta:::logit2p(TE)
      lowTE <- meta:::logit2p(lowTE)
      uppTE <- meta:::logit2p(uppTE)
    }
    ##
    if (!inherits(x, "metaprop") & x$sm %in% c("PFT", "PAS")){
      denum <- 1 + (x$sm=="PFT")
      ##
      TE    <- asin2p(TE, denum)
      lowTE <- asin2p(lowTE, denum)
      uppTE <- asin2p(uppTE, denum)
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
      ##res <- cbind(format.TE(TE), p.ci(format(lowTE), format(uppTE)))
      ##dimnames(res) <- list(x$studlab, c(x$sm, ci.lab))
      ##
      if (inherits(x, "metainf")){
        if (x$pooled=="fixed")
          cat("\nInfluential analysis (Fixed effect model)\n")
        ## writeLines(strwrap("Influential analysis (Fixed effect model)",
        ##                    prefix = "\t"))
        else
          cat("\nInfluential analysis (Random effects model)\n")
      }
      ##
      if (inherits(x, "metacum")){
        if (x$pooled=="fixed")
          cat("\nCumulative meta-analysis (Fixed effect model)\n")
        ## writeLines(strwrap("Influential analysis (Fixed effect model)",
        ##                    prefix = "\t"))
        else
          cat("\nCumulative meta-analysis (Random effects model)\n")
      }
      cat("\n")
      prmatrix(res, quote=FALSE, right=TRUE, na.print="--")
      ##
      method <- ifelse(x$method=="MH",
                       "Mantel-Haenszel method",
                       ifelse(x$method=="Peto", "Peto method",
                              ifelse(x$method=="Inverse",
                                     "Inverse variance method",
                                     x$method)))
      ##
      cat(paste("\nMethod:", method, "\n"))
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
            header=FALSE)
    
  }
  
  invisible(NULL)
}
