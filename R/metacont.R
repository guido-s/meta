metacont <- function(n.e, mean.e, sd.e, n.c, mean.c, sd.c, studlab,
                     data=NULL, subset=NULL,
                     sm=.settings$smcont,
                     level=.settings$level, level.comb=.settings$level.comb,
                     comb.fixed=.settings$comb.fixed, comb.random=.settings$comb.random,
                     hakn=.settings$hakn,
                     method.tau=.settings$method.tau, tau.preset=NULL, TE.tau=NULL,
                     tau.common=.settings$tau.common,
                     prediction=.settings$prediction, level.predict=.settings$level.predict,
                     method.bias=.settings$method.bias,
                     ##
                     title=.settings$title, complab=.settings$complab, outclab="",
                     label.e=.settings$label.e, label.c=.settings$label.c,
                     label.left=.settings$label.left, label.right=.settings$label.right,
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
  ## Catch n.e, mean.e, sd.e, n.c, mean.c, sd.c,
  ## studlab, byvar, subset from data:
  ##
  n.e <- eval(mf[[match("n.e", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  ##
  mean.e <- eval(mf[[match("mean.e", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  sd.e <- eval(mf[[match("sd.e", names(mf))]],
               data, enclos = sys.frame(sys.parent()))
  ##
  n.c <- eval(mf[[match("n.c", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  ##
  mean.c <- eval(mf[[match("mean.c", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  sd.c <- eval(mf[[match("sd.c", names(mf))]],
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
    if ((is.logical(subset) & (sum(subset) > length(n.e))) ||
        (length(subset) > length(n.e)))
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
    studlab <- seq(along=n.e)
  
  
  if (keepdata){
    if (nulldata){
      data <- data.frame(.n.e=n.e, .mean.e=mean.e, .sd.e=sd.e,
                         .n.c=n.c, .mean.c=mean.c, .sd.c=sd.c,
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
    }
    else{
      data$.n.e <- n.e
      data$.mean.e <- mean.e
      data$.sd.e <- sd.e
      data$.n.c <- n.c
      data$.mean.c <- mean.c
      data$.sd.c <- sd.c
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
    n.e <- n.e[subset]
    mean.e <- mean.e[subset]
    sd.e <- sd.e[subset]
    n.c <- n.c[subset]
    mean.c <- mean.c[subset]
    sd.c <- sd.c[subset]
    studlab <- studlab[subset]
    if (!missing.byvar)
      byvar <- byvar[subset]
  }
  
  
  k.all <- length(n.e)
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
  
  
  if (sm == "WMD"|sm=="wmd"){
    if (warn)
      warning("Effect measure '", sm, "' renamed as 'MD'.")
    sm <- "MD"
  }
  
  if (match(sm, c("MD", "SMD"), nomatch=0) == 0)
    stop("possible summary measures are \"MD\" and \"SMD\"")
  ##
  if ((length(n.c) != k.all) ||
      (length(mean.e) != k.all) || (length(sd.e) != k.all) ||
      (length(mean.c) != k.all) || (length(sd.c) != k.all))
    stop("n.e, mean.e, sd.e, n.c, mean.c and sd.c must have the same length")
  ##
  if (!(is.numeric(n.e) & is.numeric(mean.e) & is.numeric(sd.e) &
        is.numeric(n.c) & is.numeric(mean.c) & is.numeric(sd.c)))
    stop("Non-numeric value for n.e, mean.e, sd.e, n.c, mean.c or sd.c")
  ##
  npn <- n.e <= 0 | n.c <= 0
  if (any(npn) & warn)
    warning("Studies with non-positive values for n.e and/or n.c get no weight in meta-analysis.")
  ##
  ## Studies with zero sd.e and/or sd.c will be included in
  ## meta-analysis, however with zero weight
  ## (if-statement commented out by sc, 3.6.2008):
  ##
  ##if (any(sd.e <= 0 | sd.c <= 0))
  ##  stop("sd.e and sd.c must be larger than zero")
  ##
  if (length(studlab) != k.all)
    stop("Number of studies and labels are different")
  
  
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
  
  
  if (sm == "MD"){
    TE    <- ifelse(npn, NA,
                    mean.e - mean.c)
    seTE <- ifelse(npn, NA,
                   sqrt(sd.e^2/n.e + sd.c^2/n.c))
    seTE[is.na(TE)] <- NA
  }
  else if (sm == "SMD"){
    N <- n.e+n.c
    TE    <- ifelse(npn, NA,
                    (1-3/(4*N-9)) * (mean.e - mean.c) /
                    sqrt(((n.e-1)*sd.e^2 + (n.c-1)*sd.c^2)/(N-2)))
    seTE <- ifelse(npn, NA,
                   sqrt(N / (n.e*n.c) + TE^2/(2*(N-3.94))))
    seTE[is.na(TE)] <- NA
  }
  ##
  ## Studies with zero variance get zero weight in meta-analysis
  ## (added by sc, 3.6.2008):
  ##
  sel <- sd.e==0 | sd.c == 0
  ##
  if (any(sel[!is.na(sel)]) & warn)
    warning("Studies with zero values for sd.e or sd.c get no weight in meta-analysis.")
  ##
  seTE[sel] <- NA
  ##
  if (sm == "SMD")
    TE[sel] <- NA
  
  
  ##
  ## Recode integer as numeric:
  ##
  if (is.integer(n.e))    n.e    <- as.numeric(n.e)
  if (is.integer(mean.e)) mean.e <- as.numeric(mean.e)
  if (is.integer(sd.e))   sd.e   <- as.numeric(sd.e)
  if (is.integer(n.c))    n.c    <- as.numeric(n.c)
  if (is.integer(mean.c)) mean.c <- as.numeric(mean.c)
  if (is.integer(sd.c))   sd.c   <- as.numeric(sd.c)
  
  
  ##
  ## Subgroup analysis with equal tau^2:
  ##
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
  
  
  res <- list(n.e=n.e, mean.e=mean.e, sd.e=sd.e,
              n.c=n.c, mean.c=mean.c, sd.c=sd.c,
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
              sm=sm, method=m$method,
              level=level,
              level.comb=level.comb,
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
  
  class(res) <- c("metacont", "meta")
  
  res
}
