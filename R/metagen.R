metagen <- function(TE, seTE,
                    studlab, data=NULL, subset=NULL,
                    sm="",
                    level=0.95, level.comb=level,
                    comb.fixed=TRUE, comb.random=TRUE,
                    hakn=FALSE,
                    method.tau="DL", tau.preset=NULL, TE.tau=NULL,
                    tau.common=FALSE,
                    prediction=FALSE, level.predict=level,
                    method.bias="linreg",
                    n.e=NULL, n.c=NULL,
                    title="", complab="", outclab="",
                    label.e="Experimental", label.c="Control",
                    label.left="", label.right="",
                    byvar, bylab, print.byvar=TRUE,
                    keepdata=TRUE,
                    warn=TRUE
                    ){
  
  
  Ccalc <- function(x){
    res <- (sum(x  , na.rm=TRUE) -
            sum(x^2, na.rm=TRUE)/
            sum(x  , na.rm=TRUE))
    ##
    res
  }
  
  
  ##if (missing(data)) data <- NULL
  nulldata <- is.null(data)
  ##
  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch TE, seTE, studlab, byvar, subset from data:
  ##
  TE <- eval(mf[[match("TE", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  ##
  seTE <- eval(mf[[match("seTE", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  ##
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  ##
  byvar <- eval(mf[[match("byvar", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  
  
  missing.subset <- is.null(subset)
  if (!missing.subset)
    if ((is.logical(subset) & (sum(subset) > length(TE))) ||
        (length(subset) > length(TE)))
      stop("Length of subset is larger than number of studies.")
  
  
  missing.byvar <- is.null(byvar)
  if (!missing.byvar){
    byvar.name <- as.character(mf[[match("byvar", names(mf))]])
    if (length(byvar.name)>1 & byvar.name[1]=="$")
      byvar.name <- byvar.name[length(byvar.name)]
    if (length(byvar.name)>1)
      byvar.name <- "byvar"
  }
  
  
  if (!is.null(studlab))
    studlab <- as.character(studlab)
  else
    studlab <- seq(along=TE)
  
  
  if (keepdata){
    if (nulldata){
      data <- data.frame(.TE=TE, .seTE=seTE, .studlab=studlab)
      if (!missing.byvar)
        data$.byvar <- byvar
      ##
      if (!missing.subset){
        if (length(subset) == dim(data)[1])
          data$.subset <- subset
        else{
          data$.subset <- FALSE
          data$.subset[subset] <- TRUE
        }
      }
    }
    else{
      data$.TE <- TE
      data$.seTE <- seTE
      ##
      data$.studlab <- studlab
      ##
      if (!missing.byvar)
        data$.byvar <- byvar
      ##
      if (!missing.subset){
        if (length(subset) == dim(data)[1])
          data$.subset <- subset
        else{
          data$.subset <- FALSE
          data$.subset[subset] <- TRUE
        }
      }
    }
  }
  
  
  if (!missing.subset){
    TE <- TE[subset]
    seTE <- seTE[subset]
    studlab <- studlab[subset]
    if (!missing.byvar)
      byvar <- byvar[subset]
    if (!is.null(n.e))
      if (length(n.e) == length(subset))
        n.e <- n.e[subset]
      else
        warning("No subsetting for argument 'n.e' due to different length of argument 'subset'")
    if (!is.null(n.c))
      if (length(n.c) == length(subset))
        n.c <- n.c[subset]
      else
        warning("No subsetting for argument 'n.c' due to different length of argument 'subset'")
  }
  
  
  k.all <- length(TE)
  ##
  if ( k.all == 0 ) stop("No studies to combine in meta-analysis.")
  ##
  ## No meta-analysis for a single study
  ##
  if (k.all == 1){
    comb.fixed <- FALSE
    comb.random <- FALSE
    prediction <- FALSE
  }
  
  
  if ( length(seTE) != k.all )
    stop("TE and seTE must have the same length")
  ##
  if (!(is.numeric(TE) & is.numeric(seTE)))
    stop("Non-numeric value for TE or seTE")
  ##
  ## Studies with zero standard error will be included in
  ## meta-analysis, however with zero weight
  ## (if-statement commented out by sc, 3.6.2008):
  ##
  ##if ( any(seTE[!is.na(seTE)] <= 0) )
  ##  stop("seTE must be larger than zero")
  ##
  if ( length(studlab) != k.all )
    stop("Number of studies and labels differ")
  
  
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
  
  
  imethod.tau <- charmatch(tolower(method.tau),
                     c("dl", "reml", "ml", "hs", "sj", "he", "eb"), nomatch = NA)
  ##
  if (is.na(imethod.tau) || imethod.tau==0)
    stop('method.tau should be "DL", "REML", "ML", "HS", "SJ", "HE", or "EB"')
  ##
  method.tau <- c("DL", "REML", "ML", "HS", "SJ", "HE", "EB")[imethod.tau]
  
  
  ##
  ## Replace zero standard errors with NAs
  ## (added by sc, 3.6.2008):
  ##
  if (any(seTE[!is.na(seTE)] <= 0)){
    if (warn)
      warning("Zero values in seTE replaced by NAs")
    seTE[!is.na(seTE) & seTE==0] <- NA
  }
  
  
  ##
  ## Recode integer as numeric:
  ##
  if (is.integer(TE))   TE   <- as.numeric(TE)
  if (is.integer(seTE)) seTE <- as.numeric(seTE)
  
  
  k <- sum(!is.na(seTE))

  
  tau2 <- NA
  
  if (k==0){
    TE.fixed <- NA
    seTE.fixed <- NA
    zval.fixed <- NA
    pval.fixed <- NA
    lower.fixed <- NA
    upper.fixed <- NA
    w.fixed <- rep(0, k.all)
    ##
    TE.random <- NA
    seTE.random <- NA
    zval.random <- NA
    pval.random <- NA
    lower.random <- NA
    upper.random <- NA
    w.random <- rep(0, k.all)
    ##
    Q <- NA
    df.Q <- NA
    tau2 <- NA
    se.tau2 <- NA
    ##
    Cval <- NA
  }
  else{
    ##
    ## Subgroup analysis with equal tau^2:
    ##
    if (!missing.byvar & tau.common){
      sm1 <- summary(metagen(TE, seTE,
                             byvar=byvar,
                             method.tau=method.tau))
      sQ.w <- sum(sm1$Q.w)
      sk.w <- sum(sm1$k.w-1)
      sC.w <- sum(sm1$C.w)
      ##
      if (round(sQ.w, digits=18)<=sk.w) tau2 <- 0
      else tau2 <- (sQ.w-sk.w)/sC.w
      ##
      if (!is.null(tau.preset)){
        warning("Value for argument 'tau.preset' overwritten as argument 'tau.common=TRUE'")
        tau.preset <- sqrt(tau2)
      }
    }
    ##
    if (method.tau=="DL" & hakn==FALSE){
      ##
      ## Fixed effects estimate
      ## (Cooper & Hedges, 1994, p. 265-6)
      ##
      w.fixed <- 1/seTE^2
      w.fixed[is.na(w.fixed)] <- 0
      ##
      TE.fixed   <- weighted.mean(TE, w.fixed, na.rm=TRUE)
      seTE.fixed <- sqrt(1/sum(w.fixed, na.rm=TRUE))
      ##
      ci.f <- ci(TE.fixed, seTE.fixed, level=level.comb)
      zval.fixed <- ci.f$z
      pval.fixed <- ci.f$p
      lower.fixed <- ci.f$lower
      upper.fixed <- ci.f$upper
      
      ## (Cooper & Hedges (1994), p. 274-5)
      ##
      if (is.null(TE.tau))
        Q <- sum(w.fixed * (TE - TE.fixed)^2, na.rm=TRUE)
      else
        Q <- sum(w.fixed * (TE - TE.tau  )^2, na.rm=TRUE)
      ##
      df.Q <- k-1
      ##
      ## Calculate scaling factor C
      ##
      Cval <- Ccalc(w.fixed)
      ##
      ## Calculate between-study heterogeneity tau^2 
      ##
      if (is.null(tau.preset)){
        ##
        ## Following check is necessary because of argument tau.common
        ##
        if (is.na(tau2)){
          if (round(Q, digits=18)<=df.Q) tau2 <- 0
          else tau2 <- (Q-df.Q)/Cval
        }
      }
      else
        tau2 <- tau.preset^2
      ##
      se.tau2 <- NULL
      ##
      ## Random effects estimate
      ## (Cooper & Hedges (1994), p. 265, 274-5)
      ##
      w.random <- 1/(seTE^2 + tau2)
      w.random[is.na(w.random)] <- 0
      ##
      TE.random   <- weighted.mean(TE, w.random, na.rm=TRUE)
      seTE.random <- sqrt(1/sum(w.random, na.rm=TRUE))
      ##
      ci.r <- ci(TE.random, seTE.random, level=level.comb)
      zval.random <- ci.r$z
      pval.random <- ci.r$p
      lower.random <- ci.r$lower
      upper.random <- ci.r$upper
    }
    else{
      ##
      ## Check whether R package metafor is installed
      ##
      is.installed.metafor("'metagen' with argument 'method.tau' unequal to 'DL' or 'hakn=TRUE'")
      
      ##
      ## Calculate fixed effect and random effects estimates
      ##
      tres.f <- metafor::rma.uni(yi=TE, vi=seTE^2, method="FE")
      TE.fixed <- as.numeric(tres.f$b[,1])
      seTE.fixed <- tres.f$se
      ## Studies with missing treatment effect or standard error are
      ## dropped by rma.uni function (metafor package)
      ##
      ## These studies are included in meta package (however with
      ## weight zero)
      w.fixed <- rep(0, length(seTE))
      w.fixed[!(is.na(TE)|is.na(seTE))] <- 1 / tres.f$vi
      ##
      Cval <- Ccalc(w.fixed)
      ##
      zval.fixed <- tres.f$zval
      pval.fixed <- tres.f$pval
      lower.fixed <- tres.f$ci.lb
      upper.fixed <- tres.f$ci.ub
      ##
      Q <- tres.f$QE
      df.Q <- tres.f$k-tres.f$p
      ##
      ##
      if (missing(tau.preset) | is.null(tau.preset))
        tres.r <- metafor::rma.uni(yi=TE, vi=seTE^2, method=method.tau, knha=hakn)
      else
        tres.r <- metafor::rma.uni(yi=TE, vi=seTE^2, method=method.tau, knha=hakn, tau2=tau.preset^2)
      ##
      TE.random <- as.numeric(tres.r$b[,1])
      seTE.random <- tres.r$se
      ## Studies with missing standard error are dropped by rma.uni
      ## function (metafor package)
      ##
      ## These studies are included in meta package (however with
      ## weight zero)
      w.random <- rep(0, length(seTE))
      w.random[!(is.na(TE)|is.na(seTE))] <- 1 / (tres.r$vi + tres.r$tau2)
      ##
      zval.random <- tres.r$zval
      pval.random <- tres.r$pval
      lower.random <- tres.r$ci.lb
      upper.random <- tres.r$ci.ub
      ##
      df.hakn <- tres.r$k-tres.r$p
      ##
      tau2 <- tres.r$tau2
      if (method.tau %in% c("REML", "ML"))
        se.tau2 <- tres.r$se.tau2
      else
        se.tau2 <- NULL
    }
  }

  if (k>=3){
    seTE.predict <- sqrt(seTE.random^2 + tau2)
    ci.p <- ci(TE.random, seTE.predict, level.predict, k-2)
    p.lower <- ci.p$lower
    p.upper <- ci.p$upper
  }
  else{
    seTE.predict <- NA
    p.lower <- NA
    p.upper <- NA
  }
  
  
  ##
  ## Heterogeneity statistic
  ##
  if (!missing.byvar & tau.common){
    Q <- sQ.w
    df.Q <- sk.w
    Cval <- sC.w
  }
  
  
  ##
  ## Calculate H and I-Squared
  ##
  Hres  <- calcH(Q, df.Q, level.comb)
  I2res <- isquared(Q, df.Q, level.comb)
  
  
  if (!missing.byvar & tau.common)
    tau.preset <- NULL
  
  
  res <- list(TE=TE, seTE=seTE,
              studlab=studlab,
              w.fixed=w.fixed, w.random=w.random,
              ##
              TE.fixed=TE.fixed, seTE.fixed=seTE.fixed,
              lower.fixed=lower.fixed, upper.fixed=upper.fixed,
              zval.fixed=zval.fixed, pval.fixed=pval.fixed,
              ##
              TE.random=TE.random, seTE.random=seTE.random,
              lower.random=lower.random, upper.random=upper.random,
              zval.random=zval.random, pval.random=pval.random,
              ##
              seTE.predict=seTE.predict,
              lower.predict=p.lower, upper.predict=p.upper,
              level.predict=level.predict,
              ##
              k=k, Q=Q, df.Q=df.Q,
              tau=sqrt(tau2), se.tau2=se.tau2,
              C=Cval,
              ##
              H=Hres$TE,
              lower.H=Hres$lower,
              upper.H=Hres$upper,
              ##
              I2=I2res$TE,
              lower.I2=I2res$lower,
              upper.I2=I2res$upper,
              ##
              sm=sm, method="Inverse",
              level=level,
              level.comb=level.comb,
              comb.fixed=comb.fixed,
              comb.random=comb.random,
              hakn=hakn,
              df.hakn=if (hakn) df.hakn else NULL,
              method.tau=method.tau,
              tau.preset=tau.preset,
              TE.tau=if (!missing(TE.tau) & method.tau=="DL") TE.tau else NULL,
              tau.common=tau.common,
              prediction=prediction,
              method.bias=method.bias,
              n.e=n.e,
              n.c=n.c,
              title=title,
              complab=complab,
              outclab=outclab,
              label.e=label.e,
              label.c=label.c,
              label.left=label.left,
              label.right=label.right,
              data=if (keepdata) data else NULL,
              subset=if (keepdata) subset else NULL,
              warn=warn,
              call=match.call())
  ##
  if (!missing.byvar){
    res$byvar <- byvar
    res$bylab <- if (!missing(bylab) && !is.null(bylab)) bylab else byvar.name
  }
  res$print.byvar <- print.byvar
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metagen", "meta")

  res
}
