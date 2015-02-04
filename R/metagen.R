metagen <- function(TE, seTE, studlab,
                    ##
                    data=NULL, subset=NULL,
                    ##
                    sm="",
                    ##
                    level=.settings$level, level.comb=.settings$level.comb,
                    comb.fixed=.settings$comb.fixed,
                    comb.random=.settings$comb.random,
                    ##
                    hakn=.settings$hakn,
                    method.tau=.settings$method.tau,
                    tau.preset=NULL, TE.tau=NULL,
                    tau.common=.settings$tau.common,
                    ##
                    prediction=.settings$prediction,
                    level.predict=.settings$level.predict,
                    ##
                    method.bias=.settings$method.bias,
                    ##
                    n.e=NULL, n.c=NULL,
                    ##
                    backtransf=.settings$backtransf,
                    title=.settings$title, complab=.settings$complab,
                    outclab="",
                    label.e=.settings$label.e, label.c=.settings$label.c,
                    label.left=.settings$label.left,
                    label.right=.settings$label.right,
                    ##
                    byvar, bylab, print.byvar=.settings$print.byvar,
                    ##
                    keepdata=.settings$keepdata,
                    warn=.settings$warn
                    ){
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chklevel(level)
  chklevel(level.comb)
  chklogical(comb.fixed)
  chklogical(comb.random)
  ##
  chklogical(hakn)
  method.tau <- setchar(method.tau,
                        c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB"))
  chklogical(tau.common)
  ##
  chklogical(prediction)
  chklevel(level.predict)
  ##
  method.bias <- setchar(method.bias,
                         c("rank", "linreg", "mm", "count", "score", "peters"))
  ##
  chklogical(backtransf)
  chklogical(keepdata)
  ##
  ## Additional arguments / checks for metacont objects
  ##
  fun <- "metagen"
  chklogical(warn)
  chkmetafor(method.tau, fun)
  ##
  if (tau.common & method.tau=="PM"){
    warning("Argument 'method.tau' set to \"DL\" as argument tau.common=TRUE.")
    method.tau <- "DL"
  }
  ##
  if (!is.null(tau.preset) & method.tau=="PM"){
    warning("Argument 'tau.preset' not considered as argument method.tau=\"PM\".")
    tau.preset <- NULL
  }
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  nulldata <- is.null(data)
  ##
  if (nulldata)
    data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch TE, seTE, n.e, n.c from data:
  ##
  TE <- eval(mf[[match("TE", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  chknull(TE)
  k.All <- length(TE)
  ##
  seTE <- eval(mf[[match("seTE", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  chknull(seTE)
  ##
  n.e <- eval(mf[[match("n.e", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  ##
  n.c <- eval(mf[[match("n.c", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  ##
  ## Catch studlab, byvar, subset from data:
  ##
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  studlab <- setstudlab(studlab, k.All)
  ##
  byvar <- eval(mf[[match("byvar", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  missing.byvar <- is.null(byvar)
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  missing.subset <- is.null(subset)
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  chklength(seTE, k.All, fun)
  chklength(studlab, k.All, fun)
  ##
  if (!missing.byvar)
    chklength(byvar, k.All, fun)
  ##
  ## Additional checks
  ##
  if (missing.byvar & tau.common){
    warning("Value for argument 'tau.common' set to FALSE as argument 'byvar' is missing.")
    tau.common <- FALSE
  }
  if (!missing.byvar & !tau.common & !is.null(tau.preset)){
    warning("Argument 'tau.common' set to TRUE as argument tau.preset is not NULL.")
    tau.common <- TRUE
  }
  if (!is.null(n.e))
    chklength(n.e, k.All, fun)
  if (!is.null(n.c))
    chklength(n.c, k.All, fun)
  if (!missing.byvar)
    chklogical(print.byvar)
  
  
  ##
  ##
  ## (4) Subset and subgroups
  ##
  ##
  if (!missing.subset)
    if ((is.logical(subset) & (sum(subset) > k.All)) ||
        (length(subset) > k.All))
      stop("Length of subset is larger than number of studies.")
  ##  
  if (!missing.byvar){
    chkmiss(byvar)
    byvar.name <- byvarname(mf[[match("byvar", names(mf))]])
    bylab <- if (!missing(bylab) && !is.null(bylab)) bylab else byvar.name
  }
  
  
  ##
  ##
  ## (5) Store complete dataset in list object data
  ##     (if argument keepdata is TRUE)
  ##
  ##
  if (keepdata){
    if (nulldata)
      data <- data.frame(.TE=TE)
    else
      data$.TE <- TE
    ##
    data$.seTE <- seTE
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
  
  
  ##
  ##
  ## (6) Use subset for analysis
  ##
  ##
  if (!missing.subset){
    TE <- TE[subset]
    seTE <- seTE[subset]
    studlab <- studlab[subset]
    ##
    if (!missing.byvar)
      byvar <- byvar[subset]
    ##
    if (!is.null(n.e))
      n.e <- n.e[subset]
    if (!is.null(n.c))
      n.c <- n.c[subset]
  }
  ##
  ## Determine total number of studies
  ##
  k.all <- length(TE)
  ##
  if (k.all == 0)
    stop("No studies to combine in meta-analysis.")
  ##
  ## No meta-analysis for a single study
  ##
  if (k.all == 1){
    comb.fixed  <- FALSE
    comb.random <- FALSE
    prediction  <- FALSE
  }
  ##
  ## Check variable values
  ##
  chknumeric(TE)
  chknumeric(seTE, 0)
  ##
  ## Recode integer as numeric:
  ##
  TE   <- int2num(TE)
  seTE <- int2num(seTE)
  
  
  ##
  ##
  ## (7) Check standard errors
  ##
  ##
  if (any(seTE[!is.na(seTE)] <= 0)){
    if (warn)
      warning("Zero values in seTE replaced by NAs.")
    seTE[!is.na(seTE) & seTE==0] <- NA
  }
  
  
  ##
  ##
  ## (8) Do meta-analysis
  ##
  ##
  k <- sum(!is.na(seTE))
  ##
  tau2 <- NA
  ##
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
    se.tau2 <- NA
    ##
    Cval <- NA
  }
  else{
    ## Estimate tau-squared
    hc <- hetcalc(TE, seTE, method.tau, TE.tau)
    ##
    if (!missing.byvar & tau.common){
      ## Estimate common tau-squared across subgroups
      hcc <- hetcalc(TE, seTE, method.tau, TE.tau, byvar)
    }
    ##
    if (is.null(tau.preset)){
      tau2 <- hc$tau^2
      se.tau2 <- hc$se.tau2
    }
    else{
      tau2 <- tau.preset^2
      se.tau2 <- NULL
    }
    ##
    Q    <- hc$Q
    df.Q <- hc$df.Q
    Cval <- hc$Cval
    ##
    ## Fixed effect estimate (Cooper & Hedges, 1994, p. 265-6)
    ##
    w.fixed <- 1/seTE^2
    w.fixed[is.na(w.fixed) | is.na(TE)] <- 0
    ##
    TE.fixed   <- weighted.mean(TE, w.fixed, na.rm=TRUE)
    seTE.fixed <- sqrt(1/sum(w.fixed, na.rm=TRUE))
    ##
    ci.f <- ci(TE.fixed, seTE.fixed, level=level.comb)
    zval.fixed <- ci.f$z
    pval.fixed <- ci.f$p
    lower.fixed <- ci.f$lower
    upper.fixed <- ci.f$upper
    ##
    ## Random effects estimate
    ##
    if (method.tau=="PM"){
      pm <- paulemandel(TE, seTE)
      TE.random <- pm$TE.random
      seTE.random <- pm$seTE.random
      w.random <- pm$w.random
      ##
      tau2 <- pm$tau^2
    }
    else{
      ##
      ## Cooper & Hedges (1994), p. 265, 274-5
      ##
      w.random <- 1/(seTE^2 + tau2)
      w.random[is.na(w.random) | is.na(TE)] <- 0
      ##
      TE.random   <- weighted.mean(TE, w.random, na.rm=TRUE)
      seTE.random <- sqrt(1/sum(w.random, na.rm=TRUE))
    }
    ##
    ## Hartung-Knapp adjustment
    ##
    if (hakn){
      seTE.random <- sqrt(1/(k-1)*sum(w.random*(TE-TE.random)^2/
                                      sum(w.random), na.rm=TRUE))
      df.hakn <- k-1
      ci.r <- ci(TE.random, seTE.random, level=level.comb, df=df.hakn)
    }
    else
      ci.r <- ci(TE.random, seTE.random, level=level.comb)
    ##
    zval.random <- ci.r$z
    pval.random <- ci.r$p
    lower.random <- ci.r$lower
    upper.random <- ci.r$upper
  }
  ##
  ## Individual study results
  ##
  ci.study <- ci(TE, seTE, level=level)
  ##
  ## Prediction interval
  ##
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
  ## Calculate H and I-Squared
  ##
  Hres  <- calcH(Q, df.Q, level.comb)
  I2res <- isquared(Q, df.Q, level.comb)
  
  
  ##
  ##
  ## (9) Generate R object
  ##
  ##
  res <- list(studlab=studlab,
              ##
              TE=TE, seTE=seTE,
              lower=ci.study$lower, upper=ci.study$upper,
              zval=ci.study$z, pval=ci.study$p,
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
              title=title, complab=complab, outclab=outclab,
              label.e=label.e,
              label.c=label.c,
              label.left=label.left,
              label.right=label.right,
              data=if (keepdata) data else NULL,
              subset=if (keepdata) subset else NULL,
              print.byvar=print.byvar,
              warn=warn,
              call=match.call(),
              backtransf=backtransf,
              version=packageDescription("meta")$Version)
  ##
  class(res) <- c(fun, "meta")
  ##
  ## Add results from subgroup analysis
  ##
  if (!missing.byvar){
    res$byvar <- byvar
    res$bylab <- bylab
    ##
    if (!tau.common)
      res <- c(res, subgroup(res))
    else if (!is.null(tau.preset))
      res <- c(res, subgroup(res, tau.preset))
    else{
      res <- c(res, subgroup(res, hcc$tau))
      res$Q.w.random <- hcc$Q
      res$df.Q.w.random <- hcc$df.Q
    }
    ##
    res$event.e.w <- NULL
    res$event.c.w <- NULL
    res$event.w <- NULL
    res$n.w <- NULL
  }
  ##
  class(res) <- c(fun, "meta")
  
  
  res
}
