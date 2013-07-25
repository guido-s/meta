metacont <- function(n.e, mean.e, sd.e, n.c, mean.c, sd.c, studlab,
                     data=NULL, subset=NULL,
                     sm="MD",
                     level=0.95, level.comb=level,
                     comb.fixed=TRUE, comb.random=TRUE,
                     hakn=FALSE,
                     method.tau="DL", tau.preset=NULL, TE.tau=NULL,
                     tau.common=FALSE,
                     prediction=FALSE, level.predict=level,
                     method.bias="linreg",
                     title="", complab="", outclab="",
                     label.e="Experimental", label.c="Control",
                     label.left="", label.right="",
                     byvar, bylab, print.byvar=TRUE,
                     warn=TRUE
                     ){
  
  
  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  ## Catch n.e, mean.e, sd.e, n.c, mean.c, sd.c, studlab (possibly), byvar (possibly) from data:
  ##
  mf <- match.call()
  mf$data <- mf$subset <- mf$sm <- NULL
  mf$level <- mf$level.comb <- mf$level.predict <- mf$prediction <- NULL
  mf$hakn <- mf$method.tau <- mf$tau.preset <- mf$TE.tau <- mf$method.bias <- NULL
  mf[[1]] <- as.name("data.frame")
  mf <- eval(mf, data)
  ##
  ## Catch subset (possibly) from data:
  ##
  mf2 <- match.call()
  mf2$n.e <- mf2$mean.e <- mf2$sd.e <- NULL
  mf2$n.c <- mf2$mean.c <- mf2$sd.c <- NULL
  mf2$studlab <- NULL 
  mf2$data <- mf2$sm <- NULL
  mf2$level <- mf2$level.comb <- mf2$level.predict <- mf2$prediction <- NULL
  mf2$hakn <- mf2$method.tau <- mf2$tau.preset <- mf2$TE.tau <- mf2$method.bias <- NULL
  mf2$byvar <- NULL
  mf2[[1]] <- as.name("data.frame")
  ##
  mf2 <- eval(mf2, data)
  ##
  if (!is.null(mf2$subset))
    if ((is.logical(mf2$subset) & (sum(mf2$subset) > length(mf$n.e))) ||
        (length(mf2$subset) > length(mf$n.e)))
      stop("Length of subset is larger than number of studies.")
    else
      mf <- mf[mf2$subset,]
  ##if (!is.null(mf2$subset))
  ##  if (length(mf2$subset) > length(mf$n.e))
  ##    stop("Length of subset is larger than number of studies.")
  ##  else
  ##    mf <- mf[mf2$subset,]
  ##
  n.e     <- mf$n.e
  mean.e  <- mf$mean.e
  sd.e    <- mf$sd.e
  n.c     <- mf$n.c
  mean.c  <- mf$mean.c
  sd.c    <- mf$sd.c
  ##
  missing.byvar <- missing(byvar)
  if (!missing.byvar){
    byvar.name <- deparse(substitute(byvar))
    byvar <- mf$byvar
  }
  ##
  if (!missing(studlab))
    studlab <- as.character(mf$studlab)
  else
    studlab <- row.names(mf)
  
  
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
      warning("Effect measure '", sm, "' renamed as 'MD'")
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
    warning("Studies with non-positive values for n.e and/or n.c get no weight in meta-analysis")
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
  }
  else if (sm == "SMD"){
    N <- n.e+n.c
    TE    <- ifelse(npn, NA,
                    (1-3/(4*N-9)) * (mean.e - mean.c) /
                    sqrt(((n.e-1)*sd.e^2 + (n.c-1)*sd.c^2)/(N-2)))
    seTE <- ifelse(npn, NA,
                   sqrt(N / (n.e*n.c) + TE^2/(2*(N-3.94))))
  }
  ##  
  ## Studies with zero variance get zero weight in meta-analysis
  ## (added by sc, 3.6.2008):
  ##
  sel <- sd.e==0 | sd.c == 0
  ##
  if (any(sel[!is.na(sel)]) & warn)
    warning("Studies with zero values for sd.e or sd.c get no weight in meta-analysis")
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
      warning("Value for argument 'tau.preset' not considered as argument 'tau.common=TRUE'")
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
              k=m$k, Q=m$Q, tau=m$tau, se.tau2=m$se.tau2,
              C=m$C,
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
              warn=warn,
              call=match.call())
  ##
  if (!missing(byvar)){
    res$byvar <- byvar
    res$bylab <- if (!missing(bylab)) bylab else byvar.name
  }
  res$print.byvar <- print.byvar
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metacont", "meta")
  
  res
}
