metagen <- function(TE, seTE,
                    studlab, data=NULL, subset=NULL,
                    sm="",
                    level=0.95, level.comb=level,
                    comb.fixed=TRUE, comb.random=TRUE,
                    title="", complab="", outclab="",
                    label.e="Experimental", label.c="Control",
                    byvar, bylab, print.byvar=TRUE
                    ){

  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  ## Catch TE, seTE, studlab (possibly) from data:
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
  mf2$TE <- mf2$seTE <- NULL
  mf2$studlab <- NULL
  mf2$data <- mf2$sm <- NULL
  mf2$level <- mf2$level.comb <- NULL
  mf2[[1]] <- as.name("data.frame")
  ##
  mf2 <- eval(mf2, data)
  ##
  if (!is.null(mf2$subset))
    if ((is.logical(mf2$subset) & (sum(mf2$subset) > length(mf$TE))) ||
        (length(mf2$subset) > length(mf$TE)))
      stop("Length of subset is larger than number of trials.")
    else
      mf <- mf[mf2$subset,]
  ##if (!is.null(mf2$subset))
  ##  if (length(mf2$subset) > length(mf$TE))
  ##    stop("Length of subset is larger than number of trials.")
  ##  else
  ##    mf <- mf[mf2$subset,]
  ##
  TE   <- mf$TE
  seTE <- mf$seTE
  ##
  if (!missing(studlab))
    studlab <- as.character(mf$studlab)
  else
    studlab <- row.names(mf)

  
  k.all <- length(TE)
  ##
  if ( k.all == 0 ) stop("No trials to combine in meta-analysis.")

  
  if ( length(seTE) != k.all )
    stop("TE and seTE must have the same length")
  ##
  if (!(is.numeric(TE) & is.numeric(seTE)))
    stop("Non-numeric value for TE or seTE")
  ##
  ## Studies with zero standard error will be included in
  ## meta-analysis, however with zero weight
  ## (if-statement commented out by sc, 3.6.2008):
  ##
  ##if ( any(seTE[!is.na(seTE)] <= 0) )
  ##  stop("seTE must be larger than zero")
  ##
  if ( length(studlab) != k.all )
    stop("Number of studies and labels differ")
  
  
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
  ## Replace zero standard errors with NAs
  ## (added by sc, 3.6.2008):
  ##
  if (any(seTE[!is.na(seTE)] <= 0)){
    warning("Zero values in seTE replaced by NAs")
    seTE[!is.na(seTE) & seTE==0] <- NA
  }
  
  
  ##
  ## Recode integer as numeric:
  ##
  if (is.integer(TE))   TE   <- as.numeric(TE)
  if (is.integer(seTE)) seTE <- as.numeric(seTE)
  
  
  k <- sum(!is.na(seTE))

  
  if (k==0){
    TE.fixed <- NA
    seTE.fixed <- NA
    w.fixed <- NA
    ##
    TE.random <- NA
    seTE.random <- NA
    w.random <- rep(0, k.all)
    ##
    Q <- NA
    tau2 <- NA
  }
  else{

    varTE <- seTE^2

    ##
    ## Fixed effects estimate
    ## (Cooper & Hedges, 1994, p. 265-6)
    ##
    w.fixed <- 1/varTE
    w.fixed[is.na(w.fixed)] <- 0
    ##
    TE.fixed   <- weighted.mean(TE, w.fixed, na.rm=TRUE)
    seTE.fixed <- sqrt(1/sum(w.fixed, na.rm=TRUE))

    ##
    ## Heterogeneity statistic
    ## (Cooper & Hedges (1994), p. 274-5)
    ##
    Q <- sum(w.fixed * (TE - TE.fixed)^2, na.rm=TRUE)
    if (Q<=(k-1)) tau2 <- 0
    else tau2 <- (Q-(k-1))/(sum(w.fixed  , na.rm=TRUE) -
                            sum(w.fixed^2, na.rm=TRUE)/
                            sum(w.fixed  , na.rm=TRUE))
    
    ##
    ## Random effects estimate
    ## (Cooper & Hedges (1994), p. 265, 274-5)
    ##
    w.random <- 1/(varTE + tau2)
    w.random[is.na(w.random)] <- 0
    ##
    TE.random   <- weighted.mean(TE, w.random, na.rm=TRUE)
    seTE.random <- sqrt(1/sum(w.random, na.rm=TRUE))
  }


  res <- list(TE=TE, seTE=seTE,
              studlab=studlab,
              w.fixed=w.fixed, w.random=w.random,
              TE.fixed=TE.fixed, seTE.fixed=seTE.fixed, 
              TE.random=TE.random, seTE.random=seTE.random,
              k=k, Q=Q, tau=sqrt(tau2),
              sm=sm, method="Inverse",
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
  
  class(res) <- c("metagen", "meta")

  res
}
