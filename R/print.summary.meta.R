print.summary.meta <- function(x,
                               digits=max(3, .Options$digits - 3),
                               print.byvar=x$print.byvar,
                               comb.fixed=x$comb.fixed,
                               comb.random=x$comb.random,
                               header=TRUE,
                               ...){
  
  
  if (!inherits(x, "summary.meta"))
    stop("Argument 'x' must be an object of class \"summary.meta\"")
  
  
  k <- x$k
  sm <- x$sm
  sm.lab <- sm
  
  
  if (sm=="PFT"){
    sm.details <- "Freeman-Tukey double arcsine transformation"
    sm.lab <- "proportion"
  }
  else if (sm=="PAS"){
    sm.details <- "Arcsine transformation"
    sm.lab <- "proportion"
  }
  else if (sm=="PLN"){
    sm.details <- "Log transformation"
    sm.lab <- "proportion"
  }
  else if (sm=="PLOGIT"){
    sm.details <- "Logit transformation"
    sm.lab <- "proportion"
  }
  else if (sm=="PRAW"){
    sm.details <- "Untransformed proportions"
    sm.lab <- "proportion"
  }
  else if (sm=="ZCOR"){
    sm.details <- "Fisher's z transformation of correlations"
    sm.lab <- "COR"
  }
  else if (sm=="COR")
    sm.details <- "Untransformed correlations"
  else
    sm.details <- ""
  
  
  if (length(comb.fixed)==0)
    comb.fixed <- TRUE
  ##
  if (length(comb.random)==0)
    comb.random <- TRUE
  ##
  if (length(print.byvar)==0)
    print.byvar <- TRUE
  
  
  TE.fixed    <- x$fixed$TE
  lowTE.fixed <- x$fixed$lower
  uppTE.fixed <- x$fixed$upper
  ##
  TE.random    <- x$random$TE
  lowTE.random <- x$random$lower
  uppTE.random <- x$random$upper
  ##
  if (!is.null(x$bylab)){
    TE.fixed.w     <- x$within.fixed$TE
    lowTE.fixed.w  <- x$within.fixed$lower
    uppTE.fixed.w  <- x$within.fixed$upper
    TE.random.w    <- x$within.random$TE
    lowTE.random.w <- x$within.random$lower
    uppTE.random.w <- x$within.random$upper
  }
  
  
  if (sm == "RR" | sm == "OR" | sm == "HR" | sm=="PLN"){
    TE.fixed    <- exp(TE.fixed)
    lowTE.fixed <- exp(lowTE.fixed)
    uppTE.fixed <- exp(uppTE.fixed)
    ##
    TE.random <- exp(TE.random)
    lowTE.random <- exp(lowTE.random)
    uppTE.random <- exp(uppTE.random)
    ##
    if (!is.null(x$bylab)){
      TE.fixed.w     <- exp(TE.fixed.w)
      lowTE.fixed.w  <- exp(lowTE.fixed.w)
      uppTE.fixed.w  <- exp(uppTE.fixed.w)
      TE.random.w    <- exp(TE.random.w)
      lowTE.random.w <- exp(lowTE.random.w)
      uppTE.random.w <- exp(uppTE.random.w)
    }
  }
  else if (sm %in% c("PFT", "PAS")){
    denum <- 1 + (sm=="PFT")
    ##
    TE.fixed    <- asin2p(TE.fixed, denum)
    lowTE.fixed <- asin2p(lowTE.fixed, denum)
    uppTE.fixed <- asin2p(uppTE.fixed, denum)
    ##
    TE.random    <- asin2p(TE.random, denum)
    lowTE.random <- asin2p(lowTE.random, denum)
    uppTE.random <- asin2p(uppTE.random, denum)
    ##
    if (!is.null(x$bylab)){
      TE.fixed.w     <- asin2p(TE.fixed.w, denum)
      lowTE.fixed.w  <- asin2p(lowTE.fixed.w, denum)
      uppTE.fixed.w  <- asin2p(uppTE.fixed.w, denum)
      TE.random.w    <- asin2p(TE.random.w, denum)
      lowTE.random.w <- asin2p(lowTE.random.w, denum)
      uppTE.random.w <- asin2p(uppTE.random.w, denum)
    }
  }
  else if (sm=="PLOGIT"){
    TE.fixed    <- logit2p(TE.fixed)
    lowTE.fixed <- logit2p(lowTE.fixed)
    uppTE.fixed <- logit2p(uppTE.fixed)
    ##
    TE.random <- logit2p(TE.random)
    lowTE.random <- logit2p(lowTE.random)
    uppTE.random <- logit2p(uppTE.random)
    ##
    if (!is.null(x$bylab)){
      TE.fixed.w     <- logit2p(TE.fixed.w)
      lowTE.fixed.w  <- logit2p(lowTE.fixed.w)
      uppTE.fixed.w  <- logit2p(uppTE.fixed.w)
      TE.random.w    <- logit2p(TE.random.w)
      lowTE.random.w <- logit2p(lowTE.random.w)
      uppTE.random.w <- logit2p(uppTE.random.w)
    }
  }
  else if (sm=="ZCOR"){
    TE.fixed    <- z2cor(TE.fixed)
    lowTE.fixed <- z2cor(lowTE.fixed)
    uppTE.fixed <- z2cor(uppTE.fixed)
    ##
    TE.random    <- z2cor(TE.random)
    lowTE.random <- z2cor(lowTE.random)
    uppTE.random <- z2cor(uppTE.random)
    ##
    if (!is.null(x$bylab)){
      TE.fixed.w     <- z2cor(TE.fixed.w)
      lowTE.fixed.w  <- z2cor(lowTE.fixed.w)
      uppTE.fixed.w  <- z2cor(uppTE.fixed.w)
      TE.random.w    <- z2cor(TE.random.w)
      lowTE.random.w <- z2cor(lowTE.random.w)
      uppTE.random.w <- z2cor(uppTE.random.w)
    }
  }
  
  
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
  k.w <- x$k.w
  ##
  if (!is.null(x$bylab)){
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
    
    dimnames(res) <- list("", c(sm.lab, x$ci.lab, "z", "p.value"))
    
    prmatrix(res, quote=FALSE, right=TRUE, ...)
    
    if (inherits(x, "metabin")){
      method <- ifelse(x$method=="Peto",
                       "Peto method", "Inverse variance method")
      ##
      cat(paste("\nMethod:", method))
    }
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
                            c(sm.lab, x$ci.lab, "z", "p.value"))
      
      prmatrix(res, quote=FALSE, right=TRUE, ...)
      
      
      if (inherits(x, "metabin")){
        Qdata <- cbind(round(x$Q.CMH, 2), 1,
                       format.p(1-pchisq(x$Q.CMH, df=1)))
        
        dimnames(Qdata) <- list("", c("Q", "d.f.", "p.value"))
        ##
        cat("\nCochran-Mantel-Haenszel (CMH) test for overall effect: \n")
        prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
      }
    }
    else
      cat(paste("Number of trials:", k, "\n"))
    
    
    cat(paste("\nQuantifying heterogeneity:\n",
              if (x$tau^2 < 0.0001)
              "tau^2 < 0.0001"
              else
                paste("tau^2 = ",
                      format(round(x$tau^2, 4), 4, nsmall=4, scientific=FALSE), sep="")
              ,
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
        if (comb.fixed==TRUE){
        ##if (x$method!="MH" & comb.fixed==TRUE){
          if (comb.random)
            warning("Estimate from fixed effect model used in groups defined by 'byvar'")
          Q <- x$Q
          Q.w <- sum(x$Q.w, na.rm=TRUE)
          if (!is.null(x$Q.b.fixed))
            Q.b <- x$Q.b.fixed
          else
            Q.b <- Q - Q.w
          ##
          if (x$method=="MH"){
            Q.w <- NA
            Q.b <- NA
          }
          ##
          Qs  <- c(Q, Q.b,  Q.w, x$Q.w)
          Qs <- ifelse(Qs > -0.1 & Qs < 0, 0, Qs)

          df <- k-1
          df.w <- sum((x$k.w-1)[!is.na(x$Q.w)])
          df.b <- df - df.w
          ##
          if (x$method=="MH"){
            df.w <- NA
            df.b <- NA
          }
          ##
          dfs <- c(df, df.b, df.w, k.w-1)
          dfs[dfs<=0] <- NA

          pval <- 1-pchisq(Qs[1:3], df=dfs[1:3])
          
          Qdata <- cbind(ifelse(is.na(Qs),
                                "--",
                                format(round(Qs, 2))),
                         ifelse(is.na(dfs), 0, dfs),
                         c("--", "--", "--", format(TE.fixed.w)),
                         c("--", "--", "--",
                           p.ci(format(lowTE.fixed.w),
                                format(uppTE.fixed.w))),
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
                                  c("Q", "d.f.", sm.lab, x$ci.lab,
                                    "p.value"))
          
          
          if (x$method=="MH"){
            warning("Test for subgroup differences not calculated for Mantel-Haenszel method")
            Qdata <- Qdata[-(2:3),]
          }
        }
        else if (comb.random==TRUE){
          Q <- x$Q
          if (!is.null(x$Q.b.random))
            Q.b <- x$Q.b.random
          else
            Q.b <- NA
          Qs <- c(Q, Q.b, rep(NA, length(x$k.w)))
          Qs <- ifelse(Qs > -0.1 & Qs < 0, 0, Qs)
          
          df <- k-1
          df.w <- sum((x$k.w-1)[!is.na(x$Q.w)])
          df.b <- df - df.w
          dfs <- c(df, df.b, k.w-1)
          dfs[dfs<=0] <- NA
          
          pval <- 1-pchisq(Qs[1:2], df=dfs[1:2])
          
          fQs <- rmSpace(as.character(format(round(Qs, 2))))
          
          Qdata <- cbind(ifelse(fQs=="NA", "--", fQs),
                         ifelse(is.na(dfs), 0, dfs),
                         c("--", "--", format(TE.random.w)),
                         c("--", "--", p.ci(format(lowTE.random.w),
                                      format(uppTE.random.w))),
                         c(format.p(pval), rep("--", length(x$k.w))))
          
          if (print.byvar)
            bylab <- paste(x$bylab,
                           " = ", 
                           format(x$by.levs), sep="")
          else
            bylab <- format(x$by.levs)
          
          
          dimnames(Qdata) <- list(c("Total           ",
                                    "Between groups  ",
                                    bylab),
                                  c("Q", "d.f.", sm.lab, x$ci.lab,
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
    cat(paste("\nMethod:", method))
  }
  ##
  if (sm.details!="")
    cat(if (x$k.all > 1) " (",
        sm.details,
        if (x$k.all > 1) ")",
        "\n",
        sep="")
  else
    cat("\n")
  ##
  ##  if (inherits(x, "metaprop"))
  ##    if (!x$freeman.tukey)
  ##      cat("\n", if (x$k.all > 1) "        ",
  ##          "Arcsine transformation used for proportions\n",
  ##          sep="")
  ##    else
  ##      cat("\n", if (x$k.all > 1) "        ",
  ##          "Freeman-Tukey double arcsine transformation used for proportions\n",
  ##          sep="")
  ##  ##
  ##  else if (sm=="ZCOR")
  ##    cat("\n", if (x$k.all > 1) "        ",
  ##        "Fisher's z transformation of correlations\n",
  ##        sep="")
  ##  else if (sm=="COR")
  ##    cat("\n", if (x$k.all > 1) "        ",
  ##        "Untransformed correlations\n",
  ##        sep="")
  ##  else
  ##    cat("\n")
  
  invisible(NULL)
}
