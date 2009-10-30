metacont <- function(n.e, mean.e, sd.e, n.c, mean.c, sd.c, studlab,
                     data=NULL, subset=NULL,
                     sm="MD",
                     level=0.95, level.comb=level,
                     comb.fixed=TRUE, comb.random=TRUE,
                     title="", complab="", outclab="",
                     label.e="Experimental", label.c="Control",
                     byvar, bylab, print.byvar=TRUE
                     ){
  
  
  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  ## Catch n.e, mean.e, sd.e, n.c, mean.c, sd.c, studlab (possibly) from data:
  ##
  mf <- match.call()
  mf$data <- mf$subset <- mf$sm <- NULL
  mf$level <- mf$level.comb <- NULL
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
  mf2$level <- mf2$level.comb <- NULL
  mf2[[1]] <- as.name("data.frame")
  ##
  mf2 <- eval(mf2, data)
  ##
  if (!is.null(mf2$subset))
    if ((is.logical(mf2$subset) & (sum(mf2$subset) > length(mf$n.e))) ||
        (length(mf2$subset) > length(mf$n.e)))
      stop("Length of subset is larger than number of trials.")
    else
      mf <- mf[mf2$subset,]
  ##if (!is.null(mf2$subset))
  ##  if (length(mf2$subset) > length(mf$n.e))
  ##    stop("Length of subset is larger than number of trials.")
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
  if (!missing(studlab))
    studlab <- as.character(mf$studlab)
  else
    studlab <- row.names(mf)


  k.all <- length(n.e)
  ##
  if (k.all == 0) stop("No trials to combine in meta-analysis.")

  if (sm == "WMD"|sm=="wmd"){
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
  if (any(n.e <= 0 | n.c <= 0))
    stop("n.e and n.c must be positive")
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
  
  
  if (sm == "MD"){
    TE    <- mean.e - mean.c
    varTE <- sd.e^2/n.e + sd.c^2/n.c
  }
  else if (sm == "SMD"){
    N <- n.e+n.c
    TE    <- (1-3/(4*N-9)) * (mean.e - mean.c) /
      sqrt(((n.e-1)*sd.e^2 + (n.c-1)*sd.c^2)/(N-2))
    varTE <- N / (n.e*n.c) + TE^2/(2*(N-3.94))
  }
  ##  
  ## Studies with zero variance get zero weight in meta-analysis
  ## (added by sc, 3.6.2008):
  ##
  sel <- sd.e==0 | sd.c == 0
  ##
  if (any(sel[!is.na(sel)]))
    warning("Studies with zero values for sd.e or sd.c get no weight in meta-analysis")
  ##
  varTE[sel] <- NA
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
  
  
  m <- metagen(TE, sqrt(varTE))
  

  res <- list(n.e=n.e, mean.e=mean.e, sd.e=sd.e,
              n.c=n.c, mean.c=mean.c, sd.c=sd.c,
              studlab=studlab,
              TE=TE, seTE=sqrt(varTE),
              w.fixed=m$w.fixed, w.random=m$w.random,
              TE.fixed=m$TE.fixed, seTE.fixed=m$seTE.fixed,
              TE.random=m$TE.random, seTE.random=m$seTE.random,
              k=m$k, Q=m$Q, tau=m$tau,
              sm=sm, method=m$method,
              level=level,
              level.comb=level.comb,
              comb.fixed=comb.fixed,
              comb.random=comb.random,
              title=title,
              complab=complab,
              outclab=outclab,
              label.e=label.e,
              label.c=label.c,
              call=match.call())
  ##
  if (!missing(byvar)) res$byvar <- byvar
  if (!missing(bylab)) res$bylab <- bylab
  res$print.byvar <- print.byvar
  
  class(res) <- c("metacont", "meta")
  
  res
}
