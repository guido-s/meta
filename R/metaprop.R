metaprop <- function(event, n, studlab,
                     ##
                     data=NULL, subset=NULL,
                     ##
                     sm=.settings$smprop,
                     ##
                     incr=.settings$incr, allincr=.settings$allincr,
                     addincr=.settings$addincr,
                     method.ci=.settings$method.ci,
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
                     backtransf=.settings$backtransf,
                     title=.settings$title, complab=.settings$complab,
                     outclab="",
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
  ## Additional arguments / checks for metainc objects
  ##
  fun <- "metaprop"
  sm <- setchar(sm, c("PFT", "PAS", "PRAW", "PLN", "PLOGIT"))
  chklogical(allincr)
  chklogical(addincr)
  method.ci <- setchar(method.ci,
                      c("CP", "WS", "WSCC", "AC", "SA", "SACC", "NAsm"))
  chklogical(warn)
  chkmetafor(method.tau, fun)
  
  
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
  ## Catch event, n from data:
  ##
  event <- eval(mf[[match("event", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  chknull(event)
  k.All <- length(event)
  ##
  n <- eval(mf[[match("n", names(mf))]],
            data, enclos = sys.frame(sys.parent()))
  chknull(n)
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
  chklength(n, k.All, fun)
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
      data <- data.frame(.event=event)
    else
      data$.event <- event
    ##
    data$.n <- n
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
    event <- event[subset]
    n   <- n[subset]
    studlab <- studlab[subset]
    if (!missing.byvar)
      byvar <- byvar[subset]
  }
  ##
  ## Determine total number of studies
  ##
  k.all <- length(event)
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
  chknumeric(event, 0)
  chknumeric(n, 0, zero=TRUE)
  ##
  if (any(n < 10, na.rm=TRUE) & sm=="PFT")
    warning("Sample size very small (below 10) in at least one study. Accordingly, back transformation for pooled effect may be misleading for Freeman-Tukey double arcsine transformation. Please look at results for other transformations (e.g. sm='PAS' or sm='PLOGIT'), too.")
  ##
  if (any(event > n, na.rm=TRUE))
    stop("Number of events must not be larger than number of observations")
  ##
  ## Recode integer as numeric:
  ##
  event <- int2num(event)
  n     <- int2num(n)
  
  
  ##
  ##
  ## (7) Calculate results for individual studies
  ##
  ##
  sel <- switch(sm,
                PFT=rep(FALSE, length(event)),
                PAS=rep(FALSE, length(event)),
                PRAW=event == 0 | (n-event) == 0,
                PLN=event == 0 | (n-event) == 0,
                PLOGIT=event == 0 | (n-event) == 0)
  ##
  sparse <- any(sel, na.rm=TRUE)
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
  ##  
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
    ## Hartung, Knapp (2001), p. 4.00, formula (18):
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
  ## Calculate confidence intervals
  ##
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
  else{
    if (method.ci=="WS")
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
    else if (method.ci=="NAsm")
      ci.study <- ci(TE, seTE, level=level)
    ##
    lower.study <- ci.study$lower
    upper.study <- ci.study$upper
  }
  ##  
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
  
  
  ##
  ##
  ## (8) Do meta-analysis
  ##
  ##
  m <- metagen(TE, seTE, studlab,
               ##
               sm=sm,
               level=level,
               level.comb=level.comb,
               comb.fixed=comb.fixed,
               comb.random=comb.random,
               ##
               hakn=hakn,
               method.tau=method.tau,
               tau.preset=tau.preset,
               TE.tau=TE.tau,
               tau.common=FALSE,
               ##
               prediction=prediction,
               level.predict=level.predict,
               ##
               method.bias=method.bias,
               ##
               backtransf=backtransf,
               title=title, complab=complab, outclab=outclab,
               ##
               keepdata=FALSE,
               warn=warn)
  ##
  if (!missing.byvar & tau.common){
    ## Estimate common tau-squared across subgroups
    hcc <- hetcalc(TE, seTE, method.tau,
                   TE.tau,
                   byvar)
  }
  
  
  ##
  ##
  ## (9) Generate R object
  ##
  ##
  res <- list(event=event, n=n,
              incr=incr, sparse=sparse,
              allincr=allincr, addincr=addincr,
              method.ci=method.ci,
              incr.event=incr.event)
  ##
  ## Add meta-analysis results
  ## (after removing unneeded list elements)
  ##
  m$n.e <- NULL
  m$n.c <- NULL
  m$label.e <- NULL
  m$label.c <- NULL
  m$label.left <- NULL
  m$label.right <- NULL
  ##
  res <- c(res, m)
  ##
  ## Add data
  ##
  res$lower <- lower.study
  res$upper <- upper.study
  res$zval <- NAs
  res$pval <- NAs
  ##
  res$zval.fixed <- NA
  res$pval.fixed <- NA
  res$zval.random <- NA
  res$pval.random <- NA
  ##
  res$call <- match.call()
  ##
  if (keepdata){
    res$data <- data
    if (!missing.subset)
      res$subset <- subset
  }
  ##
  class(res) <- c(fun, "meta")
  ##
  ## Add results from subgroup analysis
  ##
  if (!missing.byvar){
    res$byvar <- byvar
    res$bylab <- bylab
    res$print.byvar <- print.byvar
    res$tau.common <- tau.common
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
    res$n.e.w <- NULL
    res$n.c.w <- NULL
  }
  ##
  class(res) <- c(fun, "meta")

  
  res
}
