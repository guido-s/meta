metaprop <- function(event, n, studlab,
                     data=NULL, subset=NULL,
                     sm="PFT",
                     freeman.tukey,
                     incr=0.5, allincr=FALSE, addincr=FALSE,
                     level=0.95, level.comb=level,
                     comb.fixed=TRUE, comb.random=TRUE,
                     title="", complab="", outclab="",
                     byvar, bylab, print.byvar=TRUE,
                     warn=TRUE
                     ){
  
  
  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  ## Catch event, n, studlab (possibly) from data:
  ##
  mf <- match.call()
  mf$data <- mf$subset <- mf$sm <- NULL
  mf$freeman.tukey <- mf$level <- mf$level.comb <- NULL
  mf$incr <- mf$allincr <- mf$addincr <- NULL
  mf[[1]] <- as.name("data.frame")
  mf <- eval(mf, data)
  ##
  ## Catch subset (possibly) from data:
  ##
  mf2 <- match.call()
  mf2$event <- mf2$n <- NULL
  mf2$studlab <- NULL
  mf2$data <- mf2$sm <- mf2$freeman.tukey <- NULL
  mf2$level <- mf2$level.comb <- NULL
  mf2$incr <- mf2$allincr <- mf2$addincr <- NULL
  mf2[[1]] <- as.name("data.frame")
  ##
  mf2 <- eval(mf2, data)
  ##
  if (!is.null(mf2$subset))
    if ((is.logical(mf2$subset) & (sum(mf2$subset) > length(mf$event))) ||
        (length(mf2$subset) > length(mf$event)))
      stop("Length of subset is larger than number of trials.")
    else
      mf <- mf[mf2$subset,]
  ##
  event <- mf$event
  n     <- mf$n
  ##
  if (!missing(byvar)){
    byvar.name <- deparse(substitute(byvar))
    byvar <- mf$byvar
  }
  ##
  if (!missing(studlab))
    studlab <- as.character(mf$studlab)
  else
    studlab <- row.names(mf)

  
  k.all <- length(event)
  ##
  if (k.all == 0) stop("No trials to combine in meta-analysis.")

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
    stop("sm should be \"PFT\", \"PAS\", \"PRAW\", \"PLN\", or \"PLOGIT\"")
  ##
  sm <- c("PFT", "PAS", "PRAW", "PLN", "PLOGIT")[imeth]
  ##
  if (!missing(freeman.tukey))
    warning(paste("Use of parameter freeman.tukey is deprecated.",
                  "Effect measure used:", sm))
  ##if (!is.logical(freeman.tukey))
  ##  stop("Parameter freeman.tukey must be of type 'logical'")
  
  
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
                PLN=event == 0,
                PLOGIT=event == 0 | (n-event) == 0)
  ##
  sparse <- any(sel)
  ##
  ## No need to add anything to cell counts for arcsine transformation
  ## Accordingly, no warning will be printed.
  ##
  warn2 <- !(sm %in% c("PFT", "PAS"))
  ##
  if (addincr){
    ##
    incr.event <- rep(incr, k.all)
    ##
    if (warn & warn2 & incr > 0)
      warning(paste("Increment", incr, "added to each cell frequency of all studies"))
  }
  else{
    if (sparse){
      if (allincr){
        ##
        incr.event <- rep(incr, k.all)
        ##
        if (warn & warn2 & incr > 0)
          warning(paste("Increment", incr, "added to each cell frequency of all studies"))
      }
      else{
        ##
        incr.event <- incr*sel
        ##
        if (warn & warn2 & incr > 0)
          warning(paste("Increment", incr, "added to each cell in 2x2 tables with zero cell frequencies"))
      }
    }
    else
      incr.event <- rep(0, k.all)
  }
  
  
  if (sm=="PFT"){
    TE <- asin(sqrt(event/(n+1))) + asin(sqrt((event+1)/(n+1)))
    ##varTE <- 1/(n+0.5)
    varTE <- 1/(n+1)
  }
  else if (sm=="PAS"){
    TE <- asin(sqrt(event/n))
    varTE <- 0.25*(1/n)
  }
  else if (sm=="PRAW"){
    TE <- (event+incr.event)/(n+incr.event)
    varTE <- (((event+incr.event)/(n+incr.event)) *
             (1-(event+incr.event)/(n+incr.event)) / (n+incr.event))
  }
  else if (sm=="PLN"){
    TE <- log((event+incr.event)/(n+incr.event))
    varTE <- 1/(event+incr.event) - 1/(n+incr.event)
  }
  else if (sm=="PLOGIT"){
    TE <- log(((event+incr.event)/(n+incr.event)) /
              (1-(event+incr.event)/(n+incr.event)))
    varTE <- 1/(event+incr.event) + 1/((n+incr.event)-(event+incr.event))
  }
  
  
  m <- metagen(TE, sqrt(varTE))
  
  
  res <- list(event=event, n=n,
              studlab=studlab,
              TE=TE, seTE=sqrt(varTE),
              w.fixed=m$w.fixed, w.random=m$w.random,
              TE.fixed=m$TE.fixed, seTE.fixed=m$seTE.fixed,
              TE.random=m$TE.random, seTE.random=m$seTE.random,
              k=m$k, Q=m$Q, tau=m$tau,
              sm=sm,
              method=m$method,
              level=level, level.comb=level.comb,
              comb.fixed=comb.fixed,
              comb.random=comb.random,
              title="", complab="", outclab="",
              call=match.call(),
              warn=warn)
  ##
  if (!missing(byvar)){
    res$byvar <- byvar
    res$bylab <- if (!missing(bylab)) bylab else byvar.name
  }
  res$print.byvar <- print.byvar
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metaprop", "meta")
  
  res
}
