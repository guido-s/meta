print.summary.meta <- function(x,
                               digits=max(3, .Options$digits - 3),
                               print.byvar=x$print.byvar,
                               comb.fixed=x$comb.fixed,
                               comb.random=x$comb.random,
                               header=TRUE,
                               ...){

  
  k <- x$k
  sm <- x$sm
  
  
  if (length(comb.fixed)==0){
    comb.fixed <- TRUE
  }
  ##
  if (length(comb.random)==0){
    comb.random <- TRUE
  }
  ##
  if (length(print.byvar)==0){
    print.byvar <- TRUE
  }
  
  
  if (sm == "RR" | sm == "OR" | sm == "HR"){
    x$fixed$TE <- exp(x$fixed$TE)
    x$fixed$lower <- exp(x$fixed$lower)
    x$fixed$upper <- exp(x$fixed$upper)
    ##
    x$random$TE <- exp(x$random$TE)
    x$random$lower <- exp(x$random$lower)
    x$random$upper <- exp(x$random$upper)
    ##
    if (!is.null(x$bylab)){
      x$within$TE <- exp(x$within$TE)
      x$within$lower <- exp(x$within$lower)
      x$within$upper <- exp(x$within$upper)
    }
  }
  ##
  if (inherits(x, "metaprop")){
    denum <- 1 + x$freeman.tukey
    ##
    x$fixed$TE    <- sin(x$fixed$TE/denum)^2
    x$fixed$lower <- sin(x$fixed$lower/denum)^2
    x$fixed$upper <- sin(x$fixed$upper/denum)^2
    ##
    x$random$TE    <- sin(x$random$TE/denum)^2
    x$random$lower <- sin(x$random$lower/denum)^2
    x$random$upper <- sin(x$random$upper/denum)^2
    ##
    if (!is.null(x$bylab)){
      x$within$TE    <- sin(x$within$TE/denum)^2
      x$within$lower <- sin(x$within$lower/denum)^2
      x$within$upper <- sin(x$within$upper/denum)^2
    }
  }


  TE.fixed    <- round(x$fixed$TE, digits)
  lowTE.fixed <- round(x$fixed$lower, digits)
  uppTE.fixed <- round(x$fixed$upper, digits)
  pTE.fixed <- x$fixed$p
  zTE.fixed <- round(x$fixed$z, digits)
  ##
  TE.random    <- round(x$random$TE, digits)
  lowTE.random <- round(x$random$lower, digits)
  uppTE.random <- round(x$random$upper, digits)
  pTE.random <- x$random$p
  zTE.random <- round(x$random$z, digits)
  ##
  k.w <- x$k.w
  ##
  if (!is.null(x$bylab)){
    TE.w    <- round(x$within$TE, digits)
    lowTE.w <- round(x$within$lower, digits)
    uppTE.w <- round(x$within$upper, digits)
  }
  ##
  H <- x$H$TE
  lowH <- x$H$lower
  uppH <- x$H$upper
  ##
  I2 <- x$I2$TE
  lowI2 <- x$I2$lower
  uppI2 <- x$I2$upper
  
  
  if (header){
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
  }
  
  
  if (x$k.all == 1){
    res <- cbind(TE.fixed,
                 p.ci(format(lowTE.fixed), format(uppTE.fixed)),
                 format(round(zTE.fixed,4)),
                 format.p(pTE.fixed))
    
    dimnames(res) <- list("", c(sm, x$ci.lab, "z", "p.value"))
    
    prmatrix(res, quote=FALSE, right=TRUE, ...)

    method <- ifelse(x$method=="Peto",
                     "Peto method", "Inverse variance method")
    ##
    cat(paste("\nMethod:", method, "\n"))
  }
  else{

    if (comb.fixed|comb.random){
      cat(paste("Number of trials combined:", k, "\n\n"))
    
      res <- cbind(format(c(if (comb.fixed) TE.fixed,
                            if (comb.random) TE.random)),
                   p.ci(format(c(if (comb.fixed) lowTE.fixed,
                                 if (comb.random) lowTE.random)),
                        format(c(if (comb.fixed) uppTE.fixed,
                                 if (comb.random) uppTE.random))),
                   format(round(c(if (comb.fixed) zTE.fixed,
                                  if (comb.random) zTE.random),4)),
                   format.p(c(if (comb.fixed) pTE.fixed,
                              if (comb.random) pTE.random)))
      
      dimnames(res) <- list(c(if (comb.fixed) "Fixed effect model",
                              if (comb.random) "Random effects model"),  
                            c(sm, x$ci.lab, "z", "p.value"))
      
      prmatrix(res, quote=FALSE, right=TRUE, ...)
      
      
      if (inherits(x, "metabin")){
        Qdata <- cbind(round(x$Q.CMH, 2), 1,
                       format.p(1-pchisq(x$Q.CMH, df=1)))
        
        dimnames(Qdata) <- list("", c("Q", "d.f.", "p.value"))
        ##
        cat("\nCMH-test: \n")
        prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
      }
    }
    else
      cat(paste("Number of trials:", k, "\n"))
    
    
    cat(paste("\nQuantifying heterogeneity:\n",
              "tau^2 = ",
              format(round(x$tau^2, 4), 4, nsmall=4),
              paste("; H = ", round(H, 2),
                    ifelse(k>2,
                           p.ci(round(lowH, 2), round(uppH, 2)),
                           ""),
                    "; ",
                    "I^2 = ", round(100*I2, 1), "%",
                    ifelse(k>2,
                           p.ci(paste(round(100*lowI2, 1), "%", sep=""),
                                paste(round(100*uppI2, 1), "%", sep="")),
                           ""),
                           sep=""),
              "\n", sep=""))
    ##
    ##    cat(paste("\nQuantifying heterogeneity:\n",
    ##              "tau^2 = ", round(x$tau^2, 4), 
    ##              ifelse(k>2,
    ##                     paste("; H = ", round(H, 2),
    ##                           p.ci(round(lowH, 2), round(uppH, 2)), "; ",
    ##                           "I^2 = ", round(100*I2, 1), "%",
    ##                           p.ci(paste(round(100*lowI2, 1), "%", sep=""),
    ##                                paste(round(100*uppI2, 1), "%", sep="")),
    ##                           sep=""),
    ##                     ""),
    ##              "\n", sep=""))
    

    
    if (k > 1 & (comb.fixed|comb.random)){
      
      cat("\nTest of heterogeneity:")
      
      if (is.null(x$bylab)){
        Qdata <- cbind(round(x$Q, 2), k-1,
                       format.p(1-pchisq(x$Q, df=k-1)))
        
        dimnames(Qdata) <- list("", c("Q", "d.f.", "p.value"))
      }  
      else{
        if (x$method!="MH"){
          Q <- x$Q
          Q.w <- sum(x$Q.w, na.rm=TRUE)
          Q.b <- Q - Q.w
          ##
          Qs  <- c(Q, Q.b,  Q.w, x$Q.w)
          Qs <- ifelse(Qs > -0.1 & Qs < 0, 0, Qs)

          df <- k-1
          df.w <- sum((x$k.w-1)[!is.na(x$Q.w)])
          df.b <- df - df.w
          ##
          dfs <- c(df, df.b, df.w, k.w-1)
          dfs[dfs<=0] <- NA

          pval <- 1-pchisq(Qs[1:3], df=dfs[1:3])
          
          Qdata <- cbind(format(round(Qs, 2)),
                         ifelse(is.na(dfs), 0, dfs),
                         c("--", "--", "--", format(TE.w)),
                         c("--", "--", "--",
                           p.ci(format(lowTE.w),
                                format(uppTE.w))),
                         c(format.p(pval),
                           rep("--", length(x$Q.w))))
          
          if (print.byvar)
            bylab <- paste(x$bylab,
                           " = ", 
                           format(x$by.levs), sep="")
          else
            bylab <- format(x$by.levs)
          
          
          dimnames(Qdata) <- list(c("Total           ",
                                    "Between groups  ",
                                    "Within groups   ", 
                                    bylab),
                                  c("Q", "d.f.", sm, x$ci.lab,
                                    "p.value"))
        }
        else{
          Q <- x$Q
          Qs  <- c(Q, rep(NA, length(x$k.w)))
          Qs <- ifelse(Qs > -0.1 & Qs < 0, 0, Qs)
          
          df <- k-1
          dfs <- c(df, k.w-1)
          dfs[dfs<=0] <- NA
          
          pval <- 1-pchisq(Qs[1], df=dfs[1])

          fQs <- rmSpace(as.character(format(round(Qs, 2))))
          
          Qdata <- cbind(ifelse(fQs=="NA", "--", fQs),
                         ifelse(is.na(dfs), 0, dfs),
                         c("--", format(TE.w)),
                         c("--", p.ci(format(lowTE.w), format(uppTE.w))),
                         c(format.p(pval), rep("--", length(x$k.w))))
          
          if (print.byvar)
            bylab <- paste(x$bylab,
                           " = ", 
                           format(x$by.levs), sep="")
          else
            bylab <- format(x$by.levs)
          
          
          dimnames(Qdata) <- list(c("Total           ",
                                    bylab),
                                  c("Q", "d.f.", sm, x$ci.lab,
                                    "p.value"))
        }
      }
      
      cat("\n")
      prmatrix(Qdata, quote=FALSE, right=TRUE, ...)

    }
    
    method <- ifelse(x$method=="MH",
                     "Mantel-Haenszel method",
                     ifelse(x$method=="Peto", "Peto method",
                            ifelse(x$method=="Inverse",
                                   "Inverse variance method",
                                   x$method)))
    ##
    cat(paste("\nMethod:", method, "\n"))
  }
  ##
  if (inherits(x, "metaprop"))
    if (!x$freeman.tukey)
      cat("Arcsine transformation used for proportions\n")
    else
      cat("Freeman-Tukey double arcsine transformation used for proportions\n")
  
  invisible(NULL)
}
