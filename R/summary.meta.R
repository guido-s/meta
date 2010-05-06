summary.meta <- function(object,
                         byvar=object$byvar,
                         bylab=object$bylab,
                         print.byvar=object$print.byvar,
                         bystud=FALSE,
                         level=object$level,
                         level.comb=object$level.comb,
                         comb.fixed=object$comb.fixed,
                         comb.random=object$comb.random,
                         warn=TRUE,
                         ...){
  
  
  if (!inherits(object, "meta"))
    stop("Argument 'object' must be an object of class \"meta\"")

  if (warn){
    if (inherits(object, "metainf"))
      warning("Summary method not defined for objects of class  \"metainf\"")
    ##
    if (inherits(object, "metacum"))
      warning("Summary method not defined for objects of class \"metacum\"")
  }
  
  
  k <- object$k
  Q <- object$Q
  
  
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
  ## Higgins & Thompson (2002), Statistics in Medicine, 21, 1539-58
  ##
  H <- sqrt(Q/(k-1))
  ##
  selogH <- ifelse(k>2,
                   ifelse(Q<=k,
                          sqrt(1/(2*(k-2))*(1-1/(3*(k-2)^2))),
                          0.5*(log(Q)-log(k-1))/(sqrt(2*Q)-sqrt(2*k-3))),
                   NA)
  ##
  tres <- ci(log(H), selogH, level.comb)
  ##
  ci.H <- list(TE=max(exp(tres$TE),1),
               lower=max(exp(tres$lower),1),
               upper=max(exp(tres$upper),1))

  
  func.t <- function(x) (x^2-1)/x^2
  ##
  ci.I2 <- list(TE=func.t(ci.H$TE),
                lower=func.t(ci.H$lower),
                upper=func.t(ci.H$upper))


  ci.lab <- paste(round(100*level.comb, 1), "%-CI", sep="")
  ##
  ci.study <- ci(object$TE, object$seTE, level)
  ci.f <- ci(object$TE.fixed , object$seTE.fixed , level.comb)
  ci.r <- ci(object$TE.random, object$seTE.random, level.comb)
  ##
  ## Calculate exact confidence intervals for individual studies
  ##
  ##print(!(inherits(object, "metainf")|inherits(object, "metacum")) &
  ##      inherits(object, "metaprop"))
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
      ##
      ci.f$z    <- NA
      ci.f$p    <- NA
      ##
      ci.r$z    <- NA
      ci.r$p    <- NA
    }
  }
  
  
  if (!missing(byvar) & length(object$byvar)==0){
    byvar.name <- deparse(substitute(byvar))  
    if (!is.null(object[[byvar.name]]))
      byvar <- object[[byvar.name]]
  }
  
  
  if (length(byvar)>0){
    
    if (any(is.na(byvar))) stop("Missing values in 'byvar'")
    
    by.levs <- unique(byvar)
    
    if (length(bylab)==0) bylab <- byvar.name
    
    res.w <- matrix(NA, ncol=12, nrow=length(by.levs))
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
                         warn=object$warn)
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
                          comb.random=comb.random)
      }
      ##
      if (inherits(object, "metagen")){
        meta1 <- metagen(object$TE[sel], object$seTE[sel],
                         sm=object$sm,
                         studlab=object$studlab[sel],
                         level=level, level.comb=level.comb,
                         comb.fixed=comb.fixed,
                         comb.random=comb.random)
      }
      ##
      if (inherits(object, "metaprop")){
        meta1 <- metaprop(object$event[sel], object$n[sel],
                          sm=object$sm,
                          studlab=object$studlab[sel],
                          level=level, level.comb=level.comb,
                          comb.fixed=comb.fixed,
                          comb.random=comb.random,
                          warn=object$warn)
      }
      ##
      if (inherits(object, "metacor")){
        meta1 <- metacor(object$cor[sel], object$n[sel],
                         sm=object$sm,
                         studlab=object$studlab[sel],
                         level=level, level.comb=level.comb,
                         comb.fixed=comb.fixed,
                         comb.random=comb.random)
      }
      ##
      if (bystud){
        if (print.byvar)
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
      res.w[j,] <- c(meta1$TE.fixed, meta1$seTE.fixed,
                     meta1$Q, meta1$k,
                     meta1$TE.random, meta1$seTE.random,
                     unlist(summary(meta1)$H),
                     unlist(summary(meta1)$I2))
      ##
      if (object$method=="MH")
        res.w[j,3] <- NA
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
    ##
    ci.fixed.w  <- ci(TE.fixed.w, seTE.fixed.w, level.comb)
    ci.random.w <- ci(TE.random.w, seTE.random.w, level.comb)
    ##
    if (object$method!="MH"){
      Q.b <- summary(metagen(TE.fixed.w, seTE.fixed.w))$Q
      ##
      if (comb.fixed)
        if ((round(Q-sum(Q.w, na.rm=TRUE),2) - round(Q.b,2)) != 0)
          warning(paste("Q-sum(Q.w) != Q.b\nQ.b =", round(Q.b,2)))
    }
    else
      Q.b <- NA
    ##
    if (bystud) cat("\n")
  }
  
  
  res <- list(study=ci.study,
              fixed=ci.f, random=ci.r,
              k=k, Q=Q, tau=object$tau, H=ci.H, I2=ci.I2,
              k.all=length(object$TE),
              Q.CMH=object$Q.CMH,
              sm=object$sm, method=object$method,
              call=match.call(),
              ci.lab=ci.lab,
              comb.fixed=comb.fixed,
              comb.random=comb.random)

  
  if (length(byvar)>0){
    res$within.fixed  <- ci.fixed.w
    res$within.random <- ci.random.w
    res$k.w           <- k.w
    res$Q.w           <- Q.w
    res$Q.b.fixed     <- Q.b
    res$Q.b.random    <- summary(metagen(TE.random.w, seTE.random.w))$Q
    res$H.w           <- list(TE=H.w, lower=H.w.low, upper=H.w.upp)
    res$I2.w          <- list(TE=I2.w, lower=I2.w.low, upper=I2.w.upp)
    res$bylab         <- bylab
    res$by.levs       <- by.levs
    res$within        <- "Returned list 'within' replaced by lists 'within.fixed' and 'within.random'."
  }
  
  
  class(res) <- "summary.meta"
  ##
  if (inherits(object, "metaprop")){
    res$event  <- object$event
    res$n      <- object$n
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
  }
  ##
  if (inherits(object, "trimfill")){
    res$object <- object
    ##
    class(res) <- c(class(res), "trimfill")
  }
  
  res$complab <- object$complab
  res$outclab <- object$outclab
  res$title   <- object$title
  ##
  res$print.byvar <- print.byvar
  
  res$version <- packageDescription("meta")$Version
  
  res
}
