metainc <- function(event.e, time.e, event.c, time.c, studlab,
                    data=NULL, subset=NULL, method="MH",
                    sm="IRR",
                    incr=0.5, allincr=FALSE, addincr=FALSE,
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
  
  
  ##if (missing(data)) data <- NULL
  nulldata <- is.null(data)
  ##
  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch event.e, time.e, event.c, time.c,
  ## n.e, n.c,
  ## studlab, byvar, subset from data:
  ##
  event.e <- eval(mf[[match("event.e", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  time.e <- eval(mf[[match("time.e", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  event.c <- eval(mf[[match("event.c", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  time.c <- eval(mf[[match("time.c", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  n.e <- eval(mf[[match("n.e", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  n.c <- eval(mf[[match("n.c", names(mf))]],
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
    if ((is.logical(subset) & (sum(subset) > length(event.e))) ||
        (length(subset) > length(event.e)))
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
    studlab <- seq(along=event.e)
  
  
  if (keepdata){
    if (nulldata){
      data <- data.frame(.event.e=event.e, .time.e=time.e,
                         .event.c=event.c, .time.c=time.c,
                         .studlab=studlab)
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
      ##
      if (!is.null(n.e))
        data$.n.e <- n.e
      if (!is.null(n.e))
        data$.n.c <- n.c
    }
    else{
      data$.event.e <- event.e
      data$.time.e <- time.e
      data$.n.e <- n.e
      data$.event.c <- event.c
      data$.time.c <- time.c
      data$.n.c <- n.c
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
    event.e <- event.e[subset]
    time.e <- time.e[subset]
    event.c <- event.c[subset]
    time.c <- time.c[subset]
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
  
  
  k.all <- length(event.e)
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
  
  
  if (!(is.numeric(event.e) & is.numeric(time.e)))
    stop("Non-numeric value for event.e or time.e")
  if (!(is.numeric(event.c) & is.numeric(time.c)))
    stop("Non-numeric value for event.c or time.c")
  ##
  if (any(time.e <= 0)) warning("time.e must be positive")
  if (any(time.c <= 0)) warning("time.c must be positive")
  ##
  if (any(event.e < 0))
    stop("event.e must be larger equal zero")
  if (any(event.c < 0))
    stop("event.c must be larger equal zero")
  ##
  if (length(studlab) != k.all)
    stop("Number of studies and labels are different")
  ##
  imeth <- charmatch(tolower(sm),
                     c("irr", "ird"), nomatch = NA)
  ##
  if(is.na(imeth) || imeth==0)
    stop("sm should be \"IRR\", \"IRD\"")
  ##
  sm <- c("IRR", "IRD")[imeth]
  
  
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
  
  
  imeth <- charmatch(tolower(method),
                     c("inverse", "mh", "cochran"), nomatch = NA)
  ##
  if(is.na(imeth))
    stop("method should be \"Inverse\", \"MH\", or \"Cochran\"")
  ##
  method <- c("Inverse", "MH", "Cochran")[imeth]
  
  
  ##
  ## Recode integer as numeric:
  ##
  if (is.integer(event.e)) event.e <- as.numeric(event.e)
  if (is.integer(time.e))  time.e  <- as.numeric(time.e)
  if (is.integer(event.c)) event.c <- as.numeric(event.c)
  if (is.integer(time.c))  time.c  <- as.numeric(time.c)
  
  
  ##
  ## Sparse computation
  ##
  sel <- switch(sm,
                IRD=rep(FALSE, length(event.e)),
                IRR=event.e == 0 | event.c == 0)
  ##
  sparse <- any(sel)
  ##
  if (addincr)
    incr.event <- rep(incr, k.all)
  else
    if (sparse)
      if (allincr)
        incr.event <- rep(incr, k.all)
      else
        incr.event <- incr*sel
    else
      incr.event <- rep(0, k.all)
  
  
  if (sm=="IRR"){
    TE <- log(((event.e+incr.event)/time.e) / ((event.c+incr.event)/time.c))
    seTE <- sqrt(1/(event.e+incr.event) + 1/(event.c+incr.event))
  }
  else if (sm=="IRD"){
    TE <- event.e/time.e - event.c/time.c
    seTE <- sqrt((event.e+incr.event)/time.e^2 + (event.c+incr.event)/time.c^2)
  }
  
  
  if (method=="Inverse"){
    ## Subgroup analysis with equal tau^2:
    if (!missing.byvar & tau.common){
      if (!is.null(tau.preset))
        warning("Value for argument 'tau.preset' not considered as argument 'tau.common=TRUE'.")
      ##
      sm1 <- summary(metagen(TE, seTE, byvar=byvar,
                             method.tau=method.tau,
                             tau.common=tau.common))
      sQ.w <- sum(sm1$Q.w)
      sk.w <- sum(sm1$k.w-1)
      sC.w <- sum(sm1$C.w)
      ##
      if (round(sQ.w, digits=18)<=sk.w) tau2 <- 0
      else tau2 <- (sQ.w-sk.w)/sC.w
      tau.preset <- sqrt(tau2)
    }
    ##
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
    ##
    TE.fixed <- m$TE.fixed
    seTE.fixed <- m$seTE.fixed
    w.fixed <- m$w.fixed
    zval.fixed <- m$zval.fixed
    pval.fixed <- m$pval.fixed
    lower.fixed <- m$lower.fixed
    upper.fixed <- m$upper.fixed
    ##
    TE.random <- m$TE.random
    seTE.random <- m$seTE.random
    w.random <- m$w.random
    zval.random <- m$zval.random
    pval.random <- m$pval.random
    lower.random <- m$lower.random
    upper.random <- m$upper.random
    ##
    df.hakn <- m$df.hakn
    ##
    tau2 <- m$tau^2
    se.tau2 <- m$se.tau2
  }
  ##
  else if (method == "MH"){
    ##
    if (method.tau!="DL"){
      if (warn)
        warning("DerSimonian-Laird method used to estimate between-study variance for Mantel-Haenszel method.")
      method.tau <- "DL"
    }
    ##
    if (hakn){
      if (warn)
        warning("Hartung-Knapp method not available for Mantel-Haenszel method.")
      hakn <- FALSE
    }
    ##
    if (!missing.byvar & tau.common){
      if (warn)
        warning("Argument 'tau.common' not considered for Mantel-Haenszel method.")
      tau.common <- FALSE
    }
    ##
    if (!is.null(TE.tau)){
      if (warn)
        warning("Argument 'TE.tau' not considered for Mantel-Haenszel method.")
      TE.tau <- NULL
    }
    ##
    if (!is.null(tau.preset)){
      if (warn)
        warning("Argument 'tau.preset' not considered for Mantel-Haenszel method.")
      tau.preset <- NULL
    }
    ##
    ## Greenland, Robins (1985)
    ## 
    x.k <- event.e
    y.k <- event.c
    n.k <- time.e
    m.k <- time.c
    ##
    N.k <- n.k + m.k
    t.k <- x.k + y.k
    ##
    if (sm == "IRR"){
      D <- n.k*m.k*t.k/N.k^2
      R <- x.k*m.k/N.k
      S <- y.k*n.k/N.k
      ##
      w.fixed <- S
      TE.fixed <- log(sum(R, na.rm=TRUE)/sum(S, na.rm=TRUE))
      seTE.fixed <- sqrt(sum(D, na.rm=TRUE)/(sum(R, na.rm=TRUE)*
                                  sum(S, na.rm=TRUE)))
    }
    else if (sm == "IRD"){
      L <- (x.k*m.k^2 + y.k*n.k^2)/N.k^2
      ##
      S <- n.k*m.k/N.k
      ##
      w.fixed <- S
      TE.fixed <- weighted.mean(TE, w.fixed, na.rm=TRUE)
      seTE.fixed <- sqrt(sum(L, na.rm=TRUE)/sum(S, na.rm=TRUE)^2)
    }
    ##
    ci.f <- ci(TE.fixed, seTE.fixed, level=level.comb)
    zval.fixed <- ci.f$z
    pval.fixed <- ci.f$p
    lower.fixed <- ci.f$lower
    upper.fixed <- ci.f$upper
    
    m <- metagen(TE, seTE, method.tau="DL", TE.tau=TE.fixed,
                 level=level,
                 level.comb=level.comb,
                 prediction=prediction,
                 level.predict=level.predict)
    ##
    TE.random <- m$TE.random
    seTE.random <- m$seTE.random
    w.random <- m$w.random
    ##
    zval.random <- m$zval.random
    pval.random <- m$pval.random
    lower.random <- m$lower.random
    upper.random <- m$upper.random
    ##
    tau2 <- m$tau^2
    se.tau2 <- m$se.tau2
  }
  ##
  else if (method == "Cochran"){
    ##
    if (method.tau!="DL"){
      if (warn)
        warning("DerSimonian-Laird method used to estimate between-study variance for Cochran method.")
      method.tau <- "DL"
    }
    ##
    if (hakn){
      if (warn)
        warning("Hartung-Knapp method not available for Cochran method.")
      hakn <- FALSE
    }
    ##
    if (!missing.byvar & tau.common){
      if (warn)
        warning("Argument 'tau.common' not considered for Cochran method.")
      tau.common <- FALSE
    }
    ##
    if (!is.null(TE.tau)){
      if (warn)
        warning("Argument 'TE.tau' not considered for Cochran method.")
      TE.tau <- NULL
    }
    ##
    if (!is.null(tau.preset)){
      if (warn)
        warning("Argument 'tau.preset' not considered for Cochran method.")
      tau.preset <- NULL
    }
    ##
    ## Smoking and Health - Report of the Advisory Committee to the
    ## Surgeon General of the Public Health Service,
    ## Chapter 8
    ## 
    if (sm == "IRR"){
      w.fixed <- event.c*time.e/time.c
      TE.fixed <- weighted.mean(TE, w.fixed)
      seTE.fixed <- sqrt(1/sum(event.e) + 1/sum(event.c))
    }
    else if (sm == "IRD"){
      warning("Cochran method only available for Incidence Rate Ratio (sm=\"IRR\")")
      return(NULL)
    }
    ##
    ci.f <- ci(TE.fixed, seTE.fixed, level=level.comb)
    zval.fixed <- ci.f$z
    pval.fixed <- ci.f$p
    lower.fixed <- ci.f$lower
    upper.fixed <- ci.f$upper
    
    m <- metagen(TE, seTE, method.tau="DL", TE.tau=TE.fixed,
                 level=level,
                 level.comb=level.comb,
                 prediction=prediction,
                 level.predict=level.predict)
    ##
    TE.random <- m$TE.random
    seTE.random <- m$seTE.random
    w.random <- m$w.random
    ##
    zval.random <- m$zval.random
    pval.random <- m$pval.random
    lower.random <- m$lower.random
    upper.random <- m$upper.random
    ##
    tau2 <- m$tau^2
    se.tau2 <- m$se.tau2
  }
  
  
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
  
  
  res <- list(event.e=event.e, time.e=time.e,
              event.c=event.c, time.c=time.c,
              studlab=studlab,
              TE=TE, seTE=seTE,
              w.fixed=w.fixed, w.random=w.random,
              ##
              TE.fixed=TE.fixed, seTE.fixed=seTE.fixed,
              lower.fixed=lower.fixed, upper.fixed=upper.fixed,
              zval.fixed=zval.fixed, pval.fixed=pval.fixed,
              TE.random=TE.random, seTE.random=seTE.random,
              lower.random=lower.random, upper.random=upper.random,
              zval.random=zval.random, pval.random=pval.random,
              ##
              seTE.predict=seTE.predict,
              lower.predict=p.lower, upper.predict=p.upper,
              level.predict=level.predict,
              ##
              k=m$k, Q=Q, df.Q=df.Q,
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
              sm=sm,
              method=method,
              incr=incr,
              sparse=sparse,
              allincr=allincr,
              addincr=addincr,
              incr.event=incr.event,
              level=level, level.comb=level.comb,
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
              n.e=n.e, n.c=n.c,
              title="", complab="", outclab="",
              label.e=label.e,
              label.c=label.c,
              label.left=label.left,
              label.right=label.right,
              keepdata=keepdata,
              data=if (keepdata) data else NULL,
              subset=if (keepdata) subset else NULL,
              call=match.call(),
              warn=warn)
  ##
  if (!missing.byvar){
    res$byvar <- byvar
    res$bylab <- if (!missing(bylab) && !is.null(bylab)) bylab else byvar.name
  }
  res$print.byvar <- print.byvar
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metainc", "meta")
  
  res
}
