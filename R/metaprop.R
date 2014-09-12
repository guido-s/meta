metaprop <- function(event, n, studlab,
                     data=NULL, subset=NULL,
                     sm=.settings$smprop,
                     incr=.settings$incr, allincr=.settings$allincr,
                     addincr=.settings$addincr,
                     method.ci=.settings$method.ci,
                     level=.settings$level, level.comb=.settings$level.comb,
                     comb.fixed=.settings$comb.fixed, comb.random=.settings$comb.random,
                     hakn=.settings$hakn,
                     method.tau=.settings$method.tau, tau.preset=NULL, TE.tau=NULL,
                     tau.common=.settings$tau.common,
                     prediction=.settings$prediction, level.predict=.settings$level.predict,
                     method.bias=.settings$method.bias,
                     ##
                     backtransf=.settings$backtransf,
                     title=.settings$title, complab=.settings$complab, outclab="",
                     byvar, bylab, print.byvar=.settings$print.byvar,
                     keepdata=.settings$keepdata,
                     warn=.settings$warn
                     ){
  
  
  ##if (missing(data)) data <- NULL
  nulldata <- is.null(data)
  ##
  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch event, n, studlab, byvar, subset from data:
  ##
  event <- eval(mf[[match("event", names(mf))]],
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
    if ((is.logical(subset) & (sum(subset) > length(event))) ||
        (length(subset) > length(event)))
      stop("Length of subset is larger than number of studies.")
  
  
  missing.byvar <- is.null(byvar)
  if (!missing.byvar){
    byvar.name <- as.character(mf[[match("byvar", names(mf))]])
    if (length(byvar.name)>1 & byvar.name[1]=="$")
      byvar.name <- byvar.name[length(byvar.name)]
    if (length(byvar.name)>1)
      byvar.name <- "byvar"
  }
  
  
  if (is.null(studlab))
    studlab <- seq(along=event)
  ##
  if (is.factor(studlab))
    studlab <- as.character(studlab)
  
  
  if (keepdata){
    if (nulldata){
      data <- data.frame(.event=event, .n=n, .studlab=studlab,
                         stringsAsFactors=FALSE)
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
      data$.event <- event
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
    event <- event[subset]
    n   <- n[subset]
    studlab <- studlab[subset]
    if (!missing.byvar)
      byvar <- byvar[subset]
  }
  
  
  if (is.null(method.ci))
    method.ci <- "CP"
  
  
  k.all <- length(event)
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
  
  
  if (!(is.numeric(event) & is.numeric(n)))
    stop("Non-numeric value for event or n")
  ##
  if (any(n <= 0)) stop("n must be positive")
  ##
  if (any(event < 0))
    stop("event must be larger equal zero")
  ##
  if (any(event > n)) stop("event > n")
  ##
  if (length(studlab) != k.all)
    stop("Number of studies and labels are different")
  ##
  imeth <- charmatch(tolower(sm),
                     c("pft", "pas", "praw", "pln", "plogit"), nomatch = NA)
  ##
  if(is.na(imeth) || imeth==0)
    stop("sm should be \"PLOGIT\", \"PLN\", \"PFT\", \"PAS\", or \"PRAW\"")
  ##
  sm <- c("PFT", "PAS", "PRAW", "PLN", "PLOGIT")[imeth]
  ##
  imci <- charmatch(tolower(method.ci),
                     c("cp", "ws", "wscc", "ac", "sa", "sacc", "nasm"), nomatch = NA)
  ##
  if(is.na(imci) || imci==0)
    stop("method.ci should be \"CP\", \"WS\", \"WSCC\", \"AC\", \"SA\", \"SACC\", or \"NAsm\"")
  ##
  method.ci <- c("CP", "WS", "WSCC", "AC", "SA", "SACC", "NAsm")[imci]
  ##
  if (any(n < 10) & sm=="PFT")
    warning("Sample size very small (below 10) in at least one study. Accordingly, back transformation for pooled effect may be misleading for Freeman-Tukey double arcsine transformation. Please look at results for other transformations (e.g. sm='PAS' or sm='PLOGIT'), too.")
  
  
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
  
  
  ##
  ## Recode integer as numeric:
  ##
  if (is.integer(event)) event <- as.numeric(event)
  if (is.integer(n))     n     <- as.numeric(n)
  
  
  ##
  ## Sparse computation
  ##
  sel <- switch(sm,
                PFT=rep(FALSE, length(event)),
                PAS=rep(FALSE, length(event)),
                PRAW=event == 0 | (n-event) == 0,
                PLN=event == 0 | (n-event) == 0,
                PLOGIT=event == 0 | (n-event) == 0)
  ##
  sparse <- any(sel)
  ##
  ## No need to add anything to cell counts for arcsine transformation
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
  
  
  if (sm=="PFT"){
    TE <- asin(sqrt(event/(n+1))) + asin(sqrt((event+1)/(n+1)))
    seTE <- sqrt(1/(n+0.5))
  }
  else if (sm=="PAS"){
    TE <- asin(sqrt(event/n))
    seTE <- sqrt(0.25*(1/n))
  }
  else if (sm=="PRAW"){
    TE <- event/n
    seTE <- sqrt((event+incr.event)*(n-event+incr.event)/
                 (n+2*incr.event)^3)
  }
  else if (sm=="PLN"){
    TE <- log((event+incr.event)/(n+incr.event))
    ## Hartung, Knapp (2001), p. 3880, formula (18):
    seTE <- ifelse(event == n,
                   sqrt(1/event - 1/(n+incr.event)),
                   sqrt(1/(event+incr.event) - 1/(n+incr.event))
                   )
  }
  else if (sm=="PLOGIT"){
    TE <- log((event+incr.event)/(n-event+incr.event))
    seTE <- sqrt(1/(event+incr.event) +
                 1/((n-event+incr.event)))
  }
  
  
  ##
  ## Subgroup analysis with equal tau^2:
  ##
  if (!missing.byvar & tau.common){
    if (!is.null(tau.preset))
      warning("Value for argument 'tau.preset' not considered as argument 'tau.common=TRUE'.")
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
  
  
  NAs <- rep(NA, k.all)
  ##
  if (method.ci=="CP"){
    lower.study <- upper.study <- NAs
    for (i in 1:k.all){
      cint <- binom.test(event[i], n[i], conf.level=level)
      ##
      lower.study[i] <- cint$conf.int[[1]]
      upper.study[i] <- cint$conf.int[[2]]
    }
  }
  ##
  else if (method.ci=="WS")
    ci.study <- ciWilsonScore(event, n, level=level)
  ##
  else if (method.ci=="WSCC")
    ci.study <- ciWilsonScore(event, n, level=level, correct=TRUE)
  ##
  else if (method.ci=="AC")
    ci.study <- ciAgrestiCoull(event, n, level=level)
  ##
  else if (method.ci=="SA")
    ci.study <- ciSimpleAsymptotic(event, n, level=level)
  ##
  else if (method.ci=="SACC")
    ci.study <- ciSimpleAsymptotic(event, n, level=level, correct=TRUE)
  ##
  else if (method.ci=="NAsm"){
    ci.study <- ci(TE, seTE, level=level)
  }
  
  
  if (method.ci!="CP"){
    lower.study <- ci.study$lower
    upper.study <- ci.study$upper
  }
  
  
  if (method.ci=="NAsm"){
    if (sm=="PLN"){
      lower.study <- exp(lower.study)
      upper.study <- exp(upper.study)
    }
    ##
    else if (sm=="PLOGIT"){
      lower.study <- logit2p(lower.study)
      upper.study <- logit2p(upper.study)
    }
    ##
    else if (sm=="PAS"){
      lower.study <- asin2p(lower.study, value="lower")
      upper.study <- asin2p(upper.study, value="upper")
    }
    ##
    else if (sm=="PFT"){
      lower.study <- asin2p(lower.study, n, value="lower")
      upper.study <- asin2p(upper.study, n, value="upper")
    }
    ##
    lower.study[lower.study<0] <- 0
    upper.study[upper.study>1] <- 1
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
  
  
  res <- list(event=event, n=n,
              studlab=studlab,
              TE=TE, seTE=seTE,
              lower=lower.study, upper=upper.study,
              zval=NAs, pval=NAs,
              w.fixed=m$w.fixed, w.random=m$w.random,
              ##
              TE.fixed=m$TE.fixed, seTE.fixed=m$seTE.fixed,
              lower.fixed=m$lower.fixed, upper.fixed=m$upper.fixed,
              zval.fixed=NA, pval.fixed=NA,
              TE.random=m$TE.random, seTE.random=m$seTE.random,
              lower.random=m$lower.random, upper.random=m$upper.random,
              zval.random=NA, pval.random=NA,
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
              sparse=sparse,
              incr=incr,
              allincr=allincr,
              addincr=addincr,
              incr.event=incr.event,
              method.ci=method.ci,
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
              title=title, complab=complab, outclab=outclab,
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
  
  res$backtransf <- backtransf
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metaprop", "meta")
  
  res
}
