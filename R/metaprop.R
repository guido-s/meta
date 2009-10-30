metaprop <- function(event, n, studlab,
                     data=NULL, subset=NULL,
                     freeman.tukey=TRUE,
                     level=0.95, level.comb=level,
                     comb.fixed=TRUE, comb.random=TRUE,
                     title="", complab="", outclab="",
                     byvar, bylab, print.byvar=TRUE
                     ){
  
  
  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  ## Catch event, n, studlab (possibly) from data:
  ##
  mf <- match.call()
  mf$data <- mf$subset <- mf$freeman.tukey <- NULL
  mf$level <- mf$level.comb <- NULL
  mf[[1]] <- as.name("data.frame")
  mf <- eval(mf, data)
  ##
  ## Catch subset (possibly) from data:
  ##
  mf2 <- match.call()
  mf2$event <- mf2$n <- NULL
  mf2$studlab <- NULL
  mf2$data <- mf2$freeman.tukey <- NULL
  mf2$level <- mf2$level.comb <- NULL
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
  if (!is.logical(freeman.tukey))
    stop("Parameter freeman.tukey must be of type 'logical'")
  
  
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
  
  
  method <- "Inverse"
  
  
  ##
  ## Recode integer as numeric:
  ##
  if (is.integer(event)) event <- as.numeric(event)
  if (is.integer(n))     n     <- as.numeric(n)
  
  
  if (freeman.tukey)
    TE <- asin(sqrt(event/(n+1))) + asin(sqrt((event+1)/(n+1)))
  else
    TE <- asin(sqrt(event/n))
  ##
  seTE <- 1/sqrt(n+freeman.tukey)
  
  
  m <- metagen(TE, seTE)
  
  
  res <- list(event=event, n=n,
              studlab=studlab,
              TE=TE, seTE=seTE,
              w.fixed=m$w.fixed, w.random=m$w.random,
              TE.fixed=m$TE.fixed, seTE.fixed=m$seTE.fixed,
              TE.random=m$TE.random, seTE.random=m$seTE.random,
              k=m$k, Q=m$Q, tau=m$tau,
              sm="proportion", method=m$method,
              freeman.tukey=freeman.tukey,
              level=level, level.comb=level.comb,
              comb.fixed=comb.fixed,
              comb.random=comb.random,
              title="", complab="", outclab="",
              call=match.call())
  ##
  if (!missing(byvar)) res$byvar <- byvar
  if (!missing(bylab)) res$bylab <- bylab
  res$print.byvar <- print.byvar
  
  class(res) <- c("metaprop", "meta")

  res
}
