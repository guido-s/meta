metacor <- function(cor, n, studlab,
                    data=NULL, subset=NULL,
                    sm=.settings$smcor,
                    level=.settings$level, level.comb=.settings$level.comb,
                    comb.fixed=.settings$comb.fixed, comb.random=.settings$comb.random,
                    hakn=.settings$hakn,
                    method.tau=.settings$method.tau, tau.preset=NULL, TE.tau=NULL,
                    tau.common=.settings$tau.common,
                    prediction=.settings$prediction, level.predict=.settings$level.predict,
                    method.bias=.settings$method.bias,
                    ##
                    title=.settings$title, complab=.settings$complab, outclab="",
                    byvar, bylab, print.byvar=.settings$print.byvar,
                    keepdata=.settings$keepdata
                    ){
  

  ##if (missing(data)) data <- NULL
  nulldata <- is.null(data)
  ##
  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch cor, n, studlab, byvar, subset from data:
  ##
  cor <- eval(mf[[match("cor", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  ##
  n <- eval(mf[[match("n", names(mf))]],
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
    if ((is.logical(subset) & (sum(subset) > length(cor))) ||
        (length(subset) > length(cor)))
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
    studlab <- seq(along=cor)
  
  
  if (keepdata){
    if (nulldata){
      data <- data.frame(.cor=cor, .n=n, .studlab=studlab)
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
      data$.cor <- cor
      data$.n <- n
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
    cor <- cor[subset]
    n   <- n[subset]
    studlab <- studlab[subset]
    if (!missing.byvar)
      byvar <- byvar[subset]
  }
  
  
  k.all <- length(cor)
  ##
  if (k.all == 0) stop("No studies to combine in meta-analysis.")
  ##
  ## No meta-analysis for a single study
  ##
  if (k.all == 1){
    comb.fixed <- FALSE
    comb.random <- FALSE
    prediction <- FALSE
  }
  
  
  if (!(is.numeric(cor) & is.numeric(n)))
    stop("Non-numeric value for cor or n")
  ##
  if (any(n <= 0)) stop("n must be positive")
  ##
  if (any(cor < -1) | any(cor > 1))
    stop("cor must be between -1 and 1")
  ##
  if (length(studlab) != k.all)
    stop("Number of studies and labels are different")
  ##
  imeth <- charmatch(tolower(sm), c("zcor", "cor"), nomatch = NA)
  ##
  if(is.na(imeth) || imeth==0)
    stop("sm should be \"ZCOR\" or \"COR\"")
  ##
  sm <- c("ZCOR", "COR")[imeth]
  
  
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
  
  
  if (sm=="ZCOR"){
    TE   <- 0.5 * log((1 + cor)/(1 - cor))
    seTE <- sqrt(1/(n - 3))
  }
  if (sm=="COR"){
    TE <- cor
    seTE <- sqrt((1-cor^2)^2/(n-1))
  }
  
  
  ##
  ## Subgroup analysis with equal tau^2:
  ##
  if (!missing.byvar & tau.common){
    if (!is.null(tau.preset))
      warning("Value for argument 'tau.preset' not considered as argument 'tau.common=TRUE'")
    ##
    sm1 <- summary(metagen(TE, seTE, byvar=byvar,
                           method.tau=method.tau))
    sQ.w <- sum(sm1$Q.w)
    sk.w <- sum(sm1$k.w-1)
    sC.w <- sum(sm1$C.w)
    ##
    if (round(sQ.w, digits=18)<=sk.w) tau2 <- 0
    else tau2 <- (sQ.w-sk.w)/sC.w
    tau.preset <- sqrt(tau2)
  }
  
  
  if (!is.null(tau.preset))
    m <- metagen(TE, seTE,
                 hakn=hakn, method.tau=method.tau,
                 tau.preset=tau.preset, TE.tau=TE.tau,
                 level=level,
                 level.comb=level.comb,
                 prediction=prediction,
                 level.predict=level.predict)
  else
    m <- metagen(TE, seTE,
                 hakn=hakn, method.tau=method.tau,
                 TE.tau=TE.tau,
                 level=level,
                 level.comb=level.comb,
                 prediction=prediction,
                 level.predict=level.predict)
  
  
  if (m$k>=3){
    seTE.predict <- sqrt(m$seTE.random^2 + m$tau^2)
    ci.p <- ci(m$TE.random, seTE.predict, level.predict, m$k-2)
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
  else{
    Q <- m$Q
    df.Q <- m$df.Q
    Cval <- m$C
  }
  
  
  ##
  ## Calculate H and I-Squared
  ##
  Hres  <- calcH(Q, df.Q, level.comb)
  I2res <- isquared(Q, df.Q, level.comb)
  
  
  if (!missing.byvar & tau.common)
    tau.preset <- NULL
  
  
  res <- list(cor=cor, n=n,
              studlab=studlab,
              TE=TE, seTE=seTE,
              w.fixed=m$w.fixed, w.random=m$w.random,
              ##
              TE.fixed=m$TE.fixed, seTE.fixed=m$seTE.fixed,
              lower.fixed=m$lower.fixed, upper.fixed=m$upper.fixed,
              zval.fixed=m$zval.fixed, pval.fixed=m$pval.fixed,
              TE.random=m$TE.random, seTE.random=m$seTE.random,
              lower.random=m$lower.random, upper.random=m$upper.random,
              zval.random=m$zval.random, pval.random=m$pval.random,
              ##
              seTE.predict=seTE.predict,
              lower.predict=p.lower, upper.predict=p.upper,
              level.predict=level.predict,
              ##
              k=m$k, Q=Q, df.Q=df.Q,
              tau=m$tau, se.tau2=m$se.tau2,
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
              sm=sm,
              method=m$method,
              level=level, level.comb=level.comb,
              comb.fixed=comb.fixed,
              comb.random=comb.random,
              hakn=hakn,
              df.hakn=if (hakn) m$df.hakn else NULL,
              method.tau=method.tau,
              tau.preset=tau.preset,
              TE.tau=if (!missing(TE.tau) & method.tau=="DL") TE.tau else NULL,
              tau.common=tau.common,
              prediction=prediction,
              method.bias=method.bias,
              title="", complab="", outclab="",
              data=if (keepdata) data else NULL,
              subset=if (keepdata) subset else NULL,
              call=match.call())
  ##
  if (!missing.byvar){
    res$byvar <- byvar
    res$bylab <- if (!missing(bylab) && !is.null(bylab)) bylab else byvar.name
  }
  res$print.byvar <- print.byvar
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metacor", "meta")

  res
}
