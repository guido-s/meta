metainc <- function(event.e, time.e, event.c, time.c, studlab,
                    ##
                    data=NULL, subset=NULL, method="MH",
                    ##
                    sm=.settings$sminc,
                    ##
                    incr=.settings$incr, allincr=.settings$allincr,
                    addincr=.settings$addincr,
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
  ## Additional arguments / checks for metainc objects
  ##
  fun <- "metainc"
  sm <- setchar(sm, c("IRR", "IRD"))
  method <- setchar(method, c("Inverse", "MH", "Cochran"))
  chklogical(allincr)
  chklogical(addincr)
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
  ## Catch event.e, time.e, event.c, time.c, n.e, n.c from data:
  ##
  event.e <- eval(mf[[match("event.e", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  chknull(event.e)
  k.All <- length(event.e)
  ##
  time.e <- eval(mf[[match("time.e", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknull(time.e)
  ##
  event.c <- eval(mf[[match("event.c", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  chknull(event.c)
  ##
  time.c <- eval(mf[[match("time.c", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknull(time.c)
  ##
  n.e <- eval(mf[[match("n.e", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
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
  chklength(time.e, k.All, fun)
  chklength(event.c, k.All, fun)
  chklength(time.c, k.All, fun)
  chklength(studlab, k.All, fun)
  ##
  if (!missing.byvar)
    chklength(byvar, k.All, fun)
  ##
  if (!is.null(n.e))
    chklength(n.e, k.All, fun)
  if (!is.null(n.c))
    chklength(n.c, k.All, fun)
  ##
  ## Additional checks
  ##
  if (method %in% c("MH", "Cochran")){
    ##
    mtl <- if (method=="MH") "Mantel-Haenszel method" else "Cochran method"
    ##
    if (method.tau!="DL"){
      if (warn)
        warning("DerSimonian-Laird method used to estimate between-study variance for ",
                mtl, ".")
      method.tau <- "DL"
    }
    ##
    if (hakn){
      if (warn)
        warning("Hartung-Knapp method not available for ", mtl, ".")
      hakn <- FALSE
    }
    ##
    if (!missing.byvar & tau.common){
      if (warn)
        warning("Argument 'tau.common' not considered for ", mtl, ".")
      tau.common <- FALSE
    }
    ##
    if (!is.null(TE.tau)){
      if (warn)
        warning("Argument 'TE.tau' not considered for ", mtl, ".")
      TE.tau <- NULL
    }
    ##
    if (!is.null(tau.preset)){
      if (warn)
        warning("Argument 'tau.preset' not considered for ", mtl, ".")
      tau.preset <- NULL
    }
  }
  ##
  if (missing.byvar & tau.common){
    warning("Value for argument 'tau.common' set to FALSE as argument 'byvar' is missing.")
    tau.common <- FALSE
  }
  if (!(method %in% c("MH", "Cochran"))){
    if (!missing.byvar & !tau.common & !is.null(tau.preset)){
      warning("Argument 'tau.common' set to TRUE as argument tau.preset is not NULL.")
      tau.common <- TRUE
    }
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
      data <- data.frame(.event.e=event.e)
    else
      data$.event.e <- event.e
    ##
    data$.time.e <- time.e
    data$.event.c <- event.c
    data$.time.c <- time.c
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
    ##
    if (!is.null(n.e))
      data$.n.e <- n.e
    if (!is.null(n.e))
      data$.n.c <- n.c
  }  
  
  
  ##
  ##
  ## (6) Use subset for analysis
  ##
  ##
  if (!missing.subset){
    event.e <- event.e[subset]
    time.e <- time.e[subset]
    event.c <- event.c[subset]
    time.c <- time.c[subset]
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
  k.all <- length(event.e)
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
  chknumeric(event.e, 0)
  chknumeric(time.e, 0, zero=TRUE)
  chknumeric(event.c, 0)
  chknumeric(time.c, zero=TRUE)
  ##
  ## Recode integer as numeric:
  ##
  event.e <- int2num(event.e)
  time.e  <- int2num(time.e)
  event.c <- int2num(event.c)
  time.c  <- int2num(time.c)
  
  
  ##
  ##
  ## (7) Calculate results for individual studies
  ##
  ##
  sel <- switch(sm,
                IRD=rep(FALSE, length(event.e)),
                IRR=event.e == 0 | event.c == 0)
  ##
  ## Sparse computation
  ##
  sparse <- any(sel, na.rm=TRUE)
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
  if (sm=="IRR"){
    TE <- log(((event.e+incr.event)/time.e) / ((event.c+incr.event)/time.c))
    seTE <- sqrt(1/(event.e+incr.event) + 1/(event.c+incr.event))
  }
  else if (sm=="IRD"){
    TE <- event.e/time.e - event.c/time.c
    seTE <- sqrt((event.e+incr.event)/time.e^2 + (event.c+incr.event)/time.c^2)
  }
  
  
  ##
  ##
  ## (8) Do meta-analysis
  ##
  ##
  if (method == "MH"){
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
  }
  ##
  else if (method == "Cochran"){
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
  }
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
               TE.tau=if (method=="Inverse") TE.tau else TE.fixed,
               tau.common=FALSE,
               ##
               prediction=prediction,
               level.predict=level.predict,
               ##
               method.bias=method.bias,
               ##
               backtransf=backtransf,
               title=title, complab=complab, outclab=outclab,
               label.e=label.e, label.c=label.c,
               label.left=label.left, label.right=label.right,
               ##
               keepdata=FALSE,
               warn=warn)
  ##
  if (!missing.byvar & tau.common){
    ## Estimate common tau-squared across subgroups
    hcc <- hetcalc(TE, seTE, method.tau,
                   if (method=="Inverse") TE.tau else TE.fixed,
                   byvar)
  }
  
  
  ##
  ##
  ## (9) Generate R object
  ##
  ##
  res <- list(event.e=event.e, time.e=time.e,
              event.c=event.c, time.c=time.c,
              method=method,
              incr=incr, sparse=sparse,
              allincr=allincr, addincr=addincr,
              incr.event=incr.event)
  ##
  ## Add meta-analysis results
  ## (after removing unneeded list elements)
  ##
  m$n.e <- NULL
  m$n.c <- NULL
  m$method <- NULL
  ##
  res <- c(res, m)
  ##
  ## Add data
  ##
  res$n.e <- n.e
  res$n.c <- n.c
  res$TE.tau <- TE.tau
  res$call <- match.call()
  ##
  if (method %in% c("MH", "Cochran")){
    ##
    ci.f <- ci(TE.fixed, seTE.fixed, level=level.comb)
    ##
    res$TE.fixed <- TE.fixed
    res$seTE.fixed <- seTE.fixed
    res$w.fixed <- w.fixed
    res$lower.fixed <- ci.f$lower
    res$upper.fixed <- ci.f$upper
    res$zval.fixed <- ci.f$z
    res$pval.fixed <- ci.f$p
  }
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
    res$event.w <- NULL
    res$n.w <- NULL
  }
  ##
  class(res) <- c(fun, "meta")
  
  
  res
}
