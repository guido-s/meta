metacor <- function(cor, n, studlab,
                    data=NULL, subset=NULL,
                    sm="ZCOR",
                    level=0.95, level.comb=level,
                    comb.fixed=TRUE, comb.random=TRUE,
                    title="", complab="", outclab="",
                    byvar, bylab, print.byvar=TRUE
                    ){


  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  ## Catch cor, n, studlab (possibly) from data:
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
  mf2$cor <- mf2$n <- NULL
  mf2$studlab <- NULL
  mf2$data <- mf2$sm <- NULL
  mf2$level <- mf2$level.comb <- NULL
  mf2[[1]] <- as.name("data.frame")
  ##
  mf2 <- eval(mf2, data)
  ##
  if (!is.null(mf2$subset))
    if ((is.logical(mf2$subset) & (sum(mf2$subset) > length(mf$cor))) ||
        (length(mf2$subset) > length(mf$cor)))
      stop("Length of subset is larger than number of trials.")
    else
      mf <- mf[mf2$subset,]
  ##
  cor <- mf$cor
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
  
  
  k.all <- length(cor)
  ##
  if (k.all == 0) stop("No trials to combine in meta-analysis.")

  if (!(is.numeric(cor) & is.numeric(n)))
    stop("Non-numeric value for cor or n")
  ##
  if (any(n <= 0)) stop("n must be positive")
  ##
  if (any(cor < -1) | any(cor > 1))
    stop("cor must be between -1 and 1")
  ##
  if (length(studlab) != k.all)
    stop("Number of studies and labels are different")
  ##
  imeth <- charmatch(tolower(sm), c("zcor", "cor"), nomatch = NA)
  ##
  if(is.na(imeth) || imeth==0)
    stop("sm should be \"ZCOR\" or \"COR\"")
  ##
  sm <- c("ZCOR", "COR")[imeth]
  

  if (sm=="ZCOR"){
    TE   <- 0.5 * log((1 + cor)/(1 - cor))
    seTE <- sqrt(1/(n - 3))
  }
  if (sm=="COR"){
    TE <- cor
    seTE <- sqrt((1-cor^2)^2/(n-1))
  }
  
  
  m <- metagen(TE, seTE)
  
  res <- list(cor=cor, n=n,
              studlab=studlab,
              TE=TE, seTE=seTE,
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
              call=match.call())
  ##
  if (!missing(byvar)){
    res$byvar <- byvar
    res$bylab <- if (!missing(bylab)) bylab else byvar.name
  }
  res$print.byvar <- print.byvar
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metacor", "meta")

  res
}
