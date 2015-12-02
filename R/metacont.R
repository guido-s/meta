metacont <- function(n.e, mean.e, sd.e, n.c, mean.c, sd.c, studlab,
                     ##
                     data=NULL, subset=NULL,
                     ##
                     sm=.settings$smcont,
                     ##
                     pooledvar=.settings$pooledvar,
                     method.smd=.settings$method.smd,
                     sd.glass=.settings$sd.glass,
                     exact.smd=.settings$exact.smd,
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
  backtransf <- FALSE
  chklogical(keepdata)
  ##
  ## Additional arguments / checks for metacont objects
  ##
  fun <- "metacont"
  sm <- setchar(sm, c("MD", "SMD"))
  chklogical(pooledvar)
  method.smd <- setchar(method.smd, c("Hedges", "Cohen", "Glass"))
  sd.glass <- setchar(sd.glass, c("control", "experimental"))
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
  ## Catch n.e, mean.e, sd.e, n.c, mean.c, sd.c from data:
  ##
  n.e <- eval(mf[[match("n.e", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  chknull(n.e)
  k.All <- length(n.e)
  ##
  mean.e <- eval(mf[[match("mean.e", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknull(mean.e)
  ##
  sd.e <- eval(mf[[match("sd.e", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  chknull(sd.e)
  ##
  n.c <- eval(mf[[match("n.c", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  chknull(n.c)
  ##
  mean.c <- eval(mf[[match("mean.c", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknull(mean.c)
  ##
  sd.c <- eval(mf[[match("sd.c", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  chknull(sd.c)
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
  chklength(mean.e, k.All, fun)
  chklength(sd.e, k.All, fun)
  chklength(n.c, k.All, fun)
  chklength(mean.c, k.All, fun)
  chklength(sd.c, k.All, fun)
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
      data <- data.frame(.n.e=n.e)
    else
      data$.n.e <- n.e
    ##
    data$.mean.e <- mean.e
    data$.sd.e <- sd.e
    data$.n.c <- n.c
    data$.mean.c <- mean.c
    data$.sd.c <- sd.c
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
    n.e <- n.e[subset]
    mean.e <- mean.e[subset]
    sd.e <- sd.e[subset]
    n.c <- n.c[subset]
    mean.c <- mean.c[subset]
    sd.c <- sd.c[subset]
    studlab <- studlab[subset]
    ##
    if (!missing.byvar)
      byvar <- byvar[subset]
  }
  ##
  ## Determine total number of studies
  ##
  k.all <- length(n.e)
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
  chknumeric(n.e)
  chknumeric(mean.e)
  chknumeric(sd.e)
  chknumeric(n.c)
  chknumeric(mean.c)
  chknumeric(sd.c)
  ##
  ## Recode integer as numeric:
  ##
  n.e    <- int2num(n.e)
  mean.e <- int2num(mean.e)
  sd.e   <- int2num(sd.e)
  n.c    <- int2num(n.c)
  mean.c <- int2num(mean.c)
  sd.c   <- int2num(sd.c)
  
  
  ##
  ##
  ## (7) Calculate results for individual studies
  ##
  ##
  npn1 <- npn(n.e)
  npn2 <- npn(n.c)
  npn <- npn1 | npn2
  ##
  N <- n.e+n.c
  ##
  if (any(npn) & warn)
    warning("Studies with non-positive values for n.e and/or n.c get no weight in meta-analysis.")
  ##  
  if (sm=="MD"){
    TE    <- ifelse(npn, NA, mean.e - mean.c)
    ##
    if (pooledvar)
      seTE <- ifelse(npn, NA,
                     sqrt((1/n.e + 1/n.c) *
                            ((n.e-1)*sd.e^2 + (n.c-1)*sd.c^2) / (N-2)))
    else
      seTE <- ifelse(npn, NA,
                     sqrt(sd.e^2/n.e + sd.c^2/n.c))
    ##
    seTE[is.na(TE)] <- NA
  }
  else if (sm=="SMD"){
    J <- function(x) gamma(x/2)/(sqrt(x/2)*gamma((x-1)/2))
    K <- function(x) 1 - (x-2) / (x*J(x)^2)
    ##
    if (method.smd %in% c("Hedges", "Cohen"))
      S.within <- sqrt(((n.e-1)*sd.e^2 + (n.c-1)*sd.c^2)/(N-2))
    else
      S.within <- if (sd.glass=="control") sd.c else sd.e
    ##
    smd <- ifelse(npn, NA, (mean.e - mean.c) / S.within)
    ##
    if (method.smd=="Cohen"){
      ##
      ## Borenstein et al. (2009), p. 26-27;
      ## White and Thomas (2005), p. 143
      ##
      TE <- smd
      if (exact.smd){
        J <- function(x) gamma(x/2)/(sqrt(x/2)*gamma((x-1)/2))
        K <- function(x) 1 - (x-2) / (x*J(x)^2)
        seTE <- ifelse(npn, NA,
                       sqrt(N / (n.e*n.c) + (J(N-2)*smd)^2 * K(N-2)))
      }
      else
        seTE <- ifelse(npn, NA,
                       sqrt(N / (n.e*n.c) + TE^2/(2*N)))
    }
    else if (method.smd=="Hedges"){
      ##
      ## Hedges and Olkin (1985); White and Thomas (2005), p. 143;
      ## formulae used in RevMan 5 (exact.smd=FALSE)
      ##
      if (exact.smd){
        J <- function(x) gamma(x/2)/(sqrt(x/2)*gamma((x-1)/2))
        K <- function(x) 1 - (x-2) / (x*J(x)^2)
      }
      else{
        J <- function(x) 1-3/(4*x-1)
        K <- function(x) 1/(2*(x-1.94))
      }
      ##
      TE   <- J(N-2)*smd
      seTE <- ifelse(npn, NA,
                     sqrt(N / (n.e*n.c) + TE^2 * K(N-2)))
    }
    else if (method.smd=="Glass"){
      ##
      ## see Cooper & Hedges (1994), p. 238
      ##
      n.g  <- if (sd.glass=="control") n.c else n.e
      ##
      TE <- smd
      seTE <- ifelse(npn, NA,
                     sqrt(N / (n.e*n.c) + TE^2/(2*n.g-1)))
    }
    ##
    seTE[is.na(TE)] <- NA
  }
  ##
  ## Studies with non-positive variance get zero weight in meta-analysis
  ##
  sel <- sd.e <= 0 | sd.c <= 0
  ##
  if (any(sel, na.rm=TRUE) & warn)
    warning("Studies with non-positive values for sd.e or sd.c get no weight in meta-analysis.")
  ##
  seTE[sel] <- NA
  ##
  if (sm == "SMD")
    TE[sel] <- NA
  
  
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
               label.e=label.e, label.c=label.c,
               label.left=label.left, label.right=label.right,
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
  res <- list(n.e=n.e, mean.e=mean.e, sd.e=sd.e,
              n.c=n.c, mean.c=mean.c, sd.c=sd.c,
              pooledvar=pooledvar,
              method.smd=method.smd, sd.glass=sd.glass,
              exact.smd=exact.smd)
  ##
  ## Add meta-analysis results
  ## (after removing unneeded list elements)
  ##
  m$n.e <- NULL
  m$n.c <- NULL
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
    res$n.w <- NULL
  }
  ##
  class(res) <- c(fun, "meta")
  
  
  res
}
