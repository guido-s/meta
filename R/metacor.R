metacor <- function(cor, n, studlab,
                    ##
                    data=NULL, subset=NULL,
                    ##
                    sm=.settings$smcor,
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
                    keepdata=.settings$keepdata
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
  ## Additional arguments / checks for metacor objects
  ##
  fun <- "metacor"
  sm <- setchar(sm, c("ZCOR", "COR"))
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
  ## Catch cor, n from data:
  ##
  cor <- eval(mf[[match("cor", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  chknull(cor)
  k.All <- length(cor)
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
      data <- data.frame(.cor=cor)
    else
      data$.cor <- cor
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
    cor <- cor[subset]
    n   <- n[subset]
    studlab <- studlab[subset]
    ##
    if (!missing.byvar)
      byvar <- byvar[subset]
  }
  ##
  ## Determine total number of studies
  ##
  k.all <- length(cor)
  ##
  if (k.all == 0)
    stop("No studies to combine in meta-analysis.")
  ##
  ## No meta-analysis for a single study
  ##
  if (k.all == 1){
    comb.fixed <- FALSE
    comb.random <- FALSE
    prediction <- FALSE
  }
  ##
  ## Check variable values
  ##
  chknumeric(cor, -1, 1)
  chknumeric(n, 0, zero=TRUE)
  
  
  ##
  ##
  ## (7) Calculate results for individual studies
  ##
  ##
  if (sm=="ZCOR"){
    TE   <- 0.5 * log((1 + cor)/(1 - cor))
    seTE <- sqrt(1/(n - 3))
  }
  if (sm=="COR"){
    TE <- cor
    seTE <- sqrt((1-cor^2)^2/(n-1))
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
               warn=FALSE)
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
  res <- list(cor=cor, n=n)
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
  m$warn <- NULL
  ##
  res <- c(res, m)
  ##
  ## Add data
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
    res$event.w <- NULL
    res$n.e.w <- NULL
    res$n.c.w <- NULL
  }
  ##
  class(res) <- c(fun, "meta")
  
  
  res
}
