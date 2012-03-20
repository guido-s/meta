print.summary.meta <- function(x,
                               digits=max(3, .Options$digits - 3),
                               print.byvar=x$print.byvar,
                               comb.fixed=x$comb.fixed,
                               comb.random=x$comb.random,
                               header=TRUE,
                               print.CMH=x$print.CMH,
                               bylab.nchar=35,
                               ...){
  
  
  if (!inherits(x, "summary.meta"))
    stop("Argument 'x' must be an object of class \"summary.meta\"")

  if (inherits(x, "metacum"))
    return(invisible(NULL))
  ##
  if (inherits(x, "metainf"))
    return(invisible(NULL))
  
  
  k <- x$k
  sm <- x$sm
  
  if (sm=="ZCOR")
    sm.lab <- "COR"
  else if (sm %in% c("PFT", "PAS", "PRAW", "PLN", "PLOGIT"))
    sm.lab <- "proportion"
  else
    sm.lab <- sm
  
  
  if (length(comb.fixed)==0)
    comb.fixed <- TRUE
  ##
  if (length(comb.random)==0)
    comb.random <- TRUE
  ##
  if (length(print.byvar)==0)
    print.byvar <- TRUE
  ##
  if (length(print.CMH)==0)
    print.CMH <- FALSE
  
  
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
    pval.fixed.w   <- x$within.fixed$p
    harmonic.mean.fixed.w <- x$within.fixed$harmonic.mean
    TE.random.w    <- x$within.random$TE
    lowTE.random.w <- x$within.random$lower
    uppTE.random.w <- x$within.random$upper
    pval.random.w   <- x$within.random$p
    harmonic.mean.random.w <- x$within.random$harmonic.mean
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
  else if (sm=="PFT"){
    TE.fixed    <- asin2p(TE.fixed, 1/mean(1/x$n))
    lowTE.fixed <- asin2p(lowTE.fixed, 1/mean(1/x$n))
    uppTE.fixed <- asin2p(uppTE.fixed, 1/mean(1/x$n))
    ##
    TE.random    <- asin2p(TE.random, 1/mean(1/x$n))
    lowTE.random <- asin2p(lowTE.random, 1/mean(1/x$n))
    uppTE.random <- asin2p(uppTE.random, 1/mean(1/x$n))
    ##
    if (!is.null(x$bylab)){
      TE.fixed.w     <- asin2p(TE.fixed.w, 1/harmonic.mean.fixed.w)
      lowTE.fixed.w  <- asin2p(lowTE.fixed.w, 1/harmonic.mean.fixed.w)
      uppTE.fixed.w  <- asin2p(uppTE.fixed.w, 1/harmonic.mean.fixed.w)
      TE.random.w    <- asin2p(TE.random.w, 1/harmonic.mean.random.w)
      lowTE.random.w <- asin2p(lowTE.random.w, 1/harmonic.mean.random.w)
      uppTE.random.w <- asin2p(uppTE.random.w, 1/harmonic.mean.random.w)
    }
  }
  else if (sm=="PAS"){
    TE.fixed    <- asin2p(TE.fixed)
    lowTE.fixed <- asin2p(lowTE.fixed)
    uppTE.fixed <- asin2p(uppTE.fixed)
    ##
    TE.random    <- asin2p(TE.random)
    lowTE.random <- asin2p(lowTE.random)
    uppTE.random <- asin2p(uppTE.random)
    ##
    if (!is.null(x$bylab)){
      TE.fixed.w     <- asin2p(TE.fixed.w)
      lowTE.fixed.w  <- asin2p(lowTE.fixed.w)
      uppTE.fixed.w  <- asin2p(uppTE.fixed.w)
      TE.random.w    <- asin2p(TE.random.w)
      lowTE.random.w <- asin2p(lowTE.random.w)
      uppTE.random.w <- asin2p(uppTE.random.w)
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
  
  
  if (!is.null(x$bylab))
    by.levs <- ifelse(nchar(x$by.levs) > bylab.nchar,
                      paste(substring(x$by.levs, 1, bylab.nchar-4), " ...", sep=""),#
                      x$by.levs)
  
  
  if (header)
    crtitle(x)
  
  
  if (x$k.all == 1){
    res <- cbind(TE.fixed,
                 p.ci(format(lowTE.fixed), format(uppTE.fixed)),
                 format(round(zTE.fixed,4)),
                 format.p(pTE.fixed))
    
    dimnames(res) <- list("", c(sm.lab, x$ci.lab, "z", "p.value"))
    
    prmatrix(res, quote=FALSE, right=TRUE, ...)
    
    ## Print information on summary method:
    ##
    catmeth(x$method, sm=sm, k.all=x$k.all,
            metaprop=inherits(x, "metaprop"))
  }
  else{

    if (comb.fixed|comb.random){
      cat(paste("Number of studies combined: k=", k, "\n\n", sep=""))
      
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
      
      dimnames(res) <- list(c(if (comb.fixed) "Fixed effect model",
                              if (comb.random) "Random effects model"),  
                            c(sm.lab, x$ci.lab, zlab, "p.value"))
      
      prmatrix(res, quote=FALSE, right=TRUE, ...)
      
      
      if (inherits(x, "metabin") & print.CMH){
        Qdata <- cbind(round(x$Q.CMH, 2), 1,
                       format.p(1-pchisq(x$Q.CMH, df=1)))
        
        dimnames(Qdata) <- list("", c("Q", "d.f.", "p.value"))
        ##
        cat("\nCochran-Mantel-Haenszel (CMH) test for overall effect: \n")
        prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
      }
    }
    else
      cat(paste("Number of studies: k=", k, "\n", sep=""))
    

    if (!is.na(x$tau))
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
      
      Qdata <- cbind(round(x$Q, 2), k-1,
                     format.p(1-pchisq(x$Q, df=k-1)))
      
      dimnames(Qdata) <- list("", c("Q", "d.f.", "p.value"))
      ##
      cat("\nTest of heterogeneity:\n")
      prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
      ##
      if (!is.null(x$bylab)){
        if (comb.fixed==TRUE){
          if (is.null(x$version) || as.numeric(strsplit(x$version, "-")[[1]][1]) < 1.7){
            ##
            ## R-version < 1.7-0:
            ##
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
            Qs  <- c(Q.b,  Q.w, x$Q.w)
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
            dfs <- c(df.b, df.w, k.w-1)
            dfs[dfs<=0] <- NA

            pval <- 1-pchisq(Qs[1:2], df=dfs[1:2])
            
            Qdata <- cbind(ifelse(is.na(Qs),
                                  "--",
                                  format(round(Qs, 2))),
                           ifelse(is.na(dfs), 0, dfs),
                           c("--", "--", format(TE.fixed.w)),
                           c("--", "--",
                             p.ci(format(lowTE.fixed.w),
                                  format(uppTE.fixed.w))),
                           c(format.p(pval),
                             rep("--", length(x$Q.w))))
            
            if (print.byvar)
              bylab <- paste(x$bylab,
                             " = ", 
                             format(by.levs), sep="")
            else
              bylab <- format(by.levs)
            
            
            dimnames(Qdata) <- list(c("Between groups  ",
                                      "Within groups   ", 
                                      bylab),
                                    c("Q", "d.f.", sm.lab, x$ci.lab,
                                      "p.value"))
            
            
            if (x$method=="MH"){
              warning("Test for subgroup differences for Mantel-Haenszel method only available in newer version (> 1.6-1) of R package meta")
              Qdata <- Qdata[-(2:3),]
            }
            ##
            cat("\nTest for subgroup differences (fixed effect model):\n")
            prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
          }
          else{
            ##
            ## R-version >= 1.7-0:
            ##
            Q.w <- ifelse(is.na(x$Q.w),
                          "--",
                          format(round(x$Q.w, 2)))
            I2.w <- ifelse(is.na(x$I2.w$TE),
                           "--",
                           paste(round(100*x$I2.w$TE, 1), "%", sep=""))
            tau2.w <- ifelse(x$k.w==1, "--", format.p(x$tau.w^2))
            ##
            Tdata <- cbind(format(x$k.w),
                           format(TE.fixed.w),
                           c(p.ci(format(lowTE.fixed.w), format(uppTE.fixed.w))),
                           Q.w, tau2.w, I2.w
                           ) #, format.p(pval.fixed.w))
            if (print.byvar)
              bylab <- paste(x$bylab,
                             " = ", 
                             format(by.levs), sep="")
            else
              bylab <- format(by.levs)
            dimnames(Tdata) <- list(bylab,
                                    c("  k", sm.lab, x$ci.lab,
                                      "Q", "tau^2", "I^2")
                                    ) #, "p.value"))
            cat("\nResults for subgroups (fixed effect model):\n")
            prmatrix(Tdata, quote=FALSE, right=TRUE, ...)
            
            cat("\nTest for subgroup differences (fixed effect model):\n")
            df.b <- (k-1) - sum((x$k.w-1)[!is.na(x$Q.w)])
            ##
            if (x$method=="MH"){
              Qdata <- cbind(round(x$Q.b.fixed, 2), df.b,
                             format.p(1-pchisq(x$Q.b.fixed, df=df.b)))
              dimnames(Qdata) <- list("Between groups  ",
                                      c("Q", "d.f.", "p.value"))
              prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
            }
            else{
              Qs  <- round(c(x$Q.b.fixed, sum(x$Q.w, na.rm=TRUE)), 2)
              dfs <- c(df.b, sum(k.w-1))
              Qdata <- cbind(Qs,
                             dfs,
                             format.p(1-pchisq(Qs, df=dfs)))
              dimnames(Qdata) <- list(c("Between groups", "Within groups"),
                                      c("Q", "d.f.", "p.value"))
              prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
            }
          }
        }
        ##
        ##
        if (comb.random==TRUE){
          if (is.null(x$version) || as.numeric(strsplit(x$version, "-")[[1]][1]) < 1.7){
            ##
            ## R-version < 1.7-0:
            ##
            Q <- x$Q
            if (!is.null(x$Q.b.random))
              Q.b <- x$Q.b.random
            else
              Q.b <- NA
            Qs <- c(Q.b, rep(NA, length(x$k.w)))
            Qs <- ifelse(Qs > -0.1 & Qs < 0, 0, Qs)
            
            df <- k-1
            df.w <- sum((x$k.w-1)[!is.na(x$Q.w)])
            df.b <- df - df.w
            dfs <- c(df.b, k.w-1)
            dfs[dfs<=0] <- NA
            
            pval <- 1-pchisq(Qs[1], df=dfs[1])
            
            fQs <- rmSpace(as.character(format(round(Qs, 2))))
            
            Qdata <- cbind(ifelse(fQs=="NA", "--", fQs),
                           ifelse(is.na(dfs), 0, dfs),
                           c("--", format(TE.random.w)),
                           c("--", p.ci(format(lowTE.random.w),
                                        format(uppTE.random.w))),
                           c(format.p(pval), rep("--", length(x$k.w))))
            
            if (print.byvar)
              bylab <- paste(x$bylab,
                             " = ", 
                             format(by.levs), sep="")
            else
              bylab <- format(by.levs)
            
            
            dimnames(Qdata) <- list(c("Between groups  ", bylab),
                                    c("Q", "d.f.", sm.lab, x$ci.lab,
                                      "p.value"))
            ##
            cat("\nTest for subgroup differences (random effects model):\n")
            prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
          }
          else{
            ##
            ## R-version >= 1.7-0:
            ##
            Q.w <- ifelse(is.na(x$Q.w),
                          "--",
                          format(round(x$Q.w, 2)))
            I2.w <- ifelse(is.na(x$I2.w$TE),
                           "--",
                           paste(round(100*x$I2.w$TE, 1), "%", sep=""))
            tau2.w <- ifelse(x$k.w==1, "--", format.p(x$tau.w^2))
            ##
            Tdata <- cbind(format(x$k.w),
                           format(TE.random.w),
                           c(p.ci(format(lowTE.random.w), format(uppTE.random.w))),
                           Q.w, tau2.w, I2.w
                           ) #, format.p(pval.random.w))
            if (print.byvar)
              bylab <- paste(x$bylab,
                             " = ", 
                             format(by.levs), sep="")
            else
              bylab <- format(by.levs)
            dimnames(Tdata) <- list(bylab,
                                    c("  k", sm.lab, x$ci.lab,
                                      "Q", "tau^2", "I^2")
                                    ) #, "p.value"))
            cat("\nResults for subgroups (random effects model):\n")
            prmatrix(Tdata, quote=FALSE, right=TRUE, ...)
            
            cat("\nTest for subgroup differences (random effects model):\n")
            df.b <- (k-1) - sum((x$k.w-1)[!is.na(x$Q.w)])
            ##
            Qdata <- cbind(round(x$Q.b.random, 2), df.b,
                           format.p(1-pchisq(x$Q.b.random, df=df.b)))
            dimnames(Qdata) <- list("Between groups  ",
                                    c("Q", "d.f.", "p.value"))
            prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
          }
        }
      }
    }
    
    ## Print information on summary method:
    ##
    catmeth(x$method,
            if (comb.random) x$method.tau else "",
            sm,
            x$k.all,
            !is.null(x$hakn) && (x$hakn & comb.random),
            metaprop=inherits(x, "metaprop"))
  }
  
  invisible(NULL)
}
