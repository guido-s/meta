metacor <- function(cor, n, studlab,
                    data=NULL, subset=NULL,
                    sm="ZCOR",
                    level=0.95, level.comb=level,
                    comb.fixed=TRUE, comb.random=TRUE,
                    hakn=FALSE,
                    method.tau="DL", tau.preset=NULL, TE.tau=NULL,
                    method.bias="linreg",
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
  mf$hakn <- mf$method.tau <- mf$tau.preset <- mf$TE.tau <- mf$method.bias <- NULL
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
  mf2$hakn <- mf2$method.tau <- mf2$tau.preset <- mf2$TE.tau <- mf2$method.bias <- NULL
  mf2[[1]] <- as.name("data.frame")
  ##
  mf2 <- eval(mf2, data)
  ##
  if (!is.null(mf2$subset))
    if ((is.logical(mf2$subset) & (sum(mf2$subset) > length(mf$cor))) ||
        (length(mf2$subset) > length(mf$cor)))
      stop("Length of subset is larger than number of studies.")
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
  if (k.all == 0) stop("No studies to combine in meta-analysis.")

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
  
  
  if (!is.null(tau.preset))
    m <- metagen(TE, seTE,
                 hakn=hakn, method.tau=method.tau,
                 tau.preset=tau.preset, TE.tau=TE.tau)
  else
    m <- metagen(TE, seTE,
                 hakn=hakn, method.tau=method.tau,
                 TE.tau=TE.tau)
  
  
  res <- list(cor=cor, n=n,
              studlab=studlab,
              TE=TE, seTE=seTE,
              w.fixed=m$w.fixed, w.random=m$w.random,
              TE.fixed=m$TE.fixed, seTE.fixed=m$seTE.fixed,
              lower.fixed=m$lower.fixed, upper.fixed=m$upper.fixed,
              zval.fixed=m$zval.fixed, pval.fixed=m$pval.fixed,
              TE.random=m$TE.random, seTE.random=m$seTE.random,
              lower.random=m$lower.random, upper.random=m$upper.random,
              zval.random=m$zval.random, pval.random=m$pval.random,
              k=m$k, Q=m$Q, tau=m$tau, se.tau2=m$se.tau2,
              sm=sm,
              method=m$method,
              level=level, level.comb=level.comb,
              comb.fixed=comb.fixed,
              comb.random=comb.random,
              hakn=hakn,
              df.hakn=if (hakn) m$df.hakn else NULL,
              method.tau=method.tau,
              tau.preset=tau.preset,
              TE.tau=if (!missing(TE.tau) & method.tau=="DL") TE.tau else NULL,
              method.bias=method.bias,
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
