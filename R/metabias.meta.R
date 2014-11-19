metabias.meta <- function(x, method.bias=x$method.bias,
                          plotit=FALSE, correct=FALSE,
                          k.min=10, ...){
  
  ##
  ##
  ## (1) Check for meta object
  ##
  ##
  chkclass(x, "meta")
  ##  
  if (inherits(x, "metacum"))
    stop("Test for funnel plot asymmetry not meaningful for object of class \"metacum\".")
  ##
  if (inherits(x, "metainf"))
    stop("Test for funnel plot asymmetry not meaningful for object of class \"metainf\".")
  
  
  TE <- x$TE
  seTE <- x$seTE
  n.e <- x$n.e
  n.c <- x$n.c
  data.name <- deparse(substitute(x))
  ##
  if (inherits(x, "metabin")){
    event.e <- x$event.e
    event.c <- x$event.c
    ##
    if (x$method=="Peto"){
      warning(paste("Inverse variance method used for pooling",
                    "to perform test of funnel plot asymmetry."))
      m <- metabin(event.e, n.e, event.c, n.c, sm=x$sm, method="Inverse",
                   incr=x$incr, allincr=x$allincr,
                   allstudies=x$allstudies, MH.exact=x$MH.exact,
                   RR.cochrane=x$RR.cochrane)
      TE <- m$TE
      seTE <- m$seTE
    }
  }
  
  if(length(TE) != length(seTE))
    stop("length of argument TE and seTE must be equal")
  ##
  sel <- !is.na(TE) & !is.na(seTE)
  if (length(TE) != sum(sel))
    warning(paste(length(TE) - sum(sel),
                  "observation(s) dropped due to missing values"))
  ##
  TE <- TE[sel]
  seTE <- seTE[sel]
  ##
  if (inherits(x, "metabin")){
    n.e <- n.e[sel]
    n.c <- n.c[sel]
    event.e <- event.e[sel]
    event.c <- event.c[sel]
  }
  ##
  k <- length(TE)
  
  ##
  ## Set default for method.bias (see Sterne et al., BMJ 2011;343:d4002):
  ## - "score" for OR as effect measure
  ## - "linreg" otherwise
  ##
  if (length(method.bias)==0)
    if (inherits(x, "metabin") & x$sm=="OR")
      method.bias <- "score"
    else
      method.bias <- "linreg"
  
  imeth <- charmatch(method.bias,
                     c("rank", "linreg", "mm", "count", "score", "peters"),
                     nomatch = NA)
  if(is.na(imeth) | imeth==0)
    stop("method.bias should be \"rank\", \"linreg\", \"mm\", \"count\", \"score\", or \"peters\"")
  ##
  method.bias <- c("rank", "linreg", "mm", "count", "score", "peters")[imeth]
  
  
  if (k < k.min | k < 3)
    res <- list(k=k, k.min=k.min)
  else if (length(x$byvar)!=0)
    res <- list(subgroup=TRUE)
  else{
    if (method.bias=="rank"){
      if (length(unique(seTE))==1)
        stop("Test for small-study effects not feasible as all studies have same precision")
      ##
      ## Begg und Mazumdar (1994), Biometrics, 50, 1088-1101
      ##
      m <- metagen(TE, seTE)
      TE.fixed <- m$TE.fixed
      seTE.fixed <- m$seTE.fixed
      ##
      varTE.s <- seTE^2 - seTE.fixed^2
      TE.s <- (TE - TE.fixed) / sqrt(varTE.s)
      ##
      ktau <- kentau(TE.s, seTE^2, correct=correct)
      ##
      res <- list(estimate=c(ktau$ks, ktau$se.ks),
                  statistic=ktau$ks/ktau$se.ks,
                  p.value=ktau$p.value)
      names(res$statistic) <- "z"
      names(res$estimate) <- c("ks", "se.ks")
    }
    else if (method.bias=="linreg" | method.bias=="mm" | method.bias=="score" |
             method.bias=="peters"){
      
      if (method.bias=="linreg"){
        if (length(unique(seTE))==1)
          stop("Test for small-study effects not feasible as all studies have same precision")
        ##
        ## Egger, Smith, Schneider, Minder (1997), BMJ, 315, 629-34
        ##
        lreg <- linregcore(seTE, TE, 1/seTE^2)
        se.bias <- lreg$se.slope
      }
      else if (method.bias=="mm"){
        if (length(unique(seTE))==1)
          stop("Test for small-study effects not feasible as all studies have same precision")
        ##
        ## Thompson und Sharp (1999), Stat Med, 18, 2693-2708
        ##
        fit1 <- linregcore(1/seTE, TE/seTE)
        Q <- sum((TE/seTE-fit1$intercept - fit1$slope/seTE)^2)
        ##
        x <- seTE
        y <- TE
        w <- 1/seTE^2
        ##
        tau2 <- (Q-(k-2))/ (sum(w) - (sum(w^2)*sum(w*x^2) -
                                      2*sum(w^2*x)*sum(w*x)+
                                      sum(w)*sum(w^2*x^2)) /
                            (sum(w)*sum(w*x^2)-(sum(w*x))^2))
        ##
        tau2 <- ifelse(tau2<0, 0, tau2)
        ##
        lreg <- linregcore(seTE, TE, 1/(seTE^2+tau2))
        se.bias <- lreg$se.slope/sqrt(lreg$MSE.w)
      }
      else if (method.bias=="score"){
        ##
        ## Harbord et al. (2006), Stat Med
        ##
        if (inherits(x, "metabin")){
          Z <- event.e - (event.e+event.c)*n.e/(n.e+n.c)
          V <- (n.e*n.c*(event.e+event.c)*
                (n.e-event.e+n.c-event.c)/((n.e+n.c)^2*((n.e+n.c)-1)))
          ##
          TE.score <- Z/V
          seTE.score <- 1/sqrt(V)
          ##
          lreg <- linregcore(seTE.score, TE.score, 1/seTE.score^2)
          se.bias <- lreg$se.slope
        }
        else{
          stop(paste("method.bias '", method.bias, "' only defined for meta-analysis with binary outcome data (function 'metabin')", sep=""))
        }
      }
      else if (method.bias=="peters"){
        ##
        ## Peters et al. (2006), JAMA
        ##
        if (inherits(x, "metabin")){
          seTE.peters <- sqrt(1/(event.e+event.c) +
                              1/(n.e-event.e + n.c-event.c))
          ##
          lreg <- linregcore(1/(n.e+n.c), TE, 1/seTE.peters^2)
          se.bias <- lreg$se.slope
        }
        else{
          stop(paste("method.bias '", method.bias, "' only defined for meta-analysis with binary outcome data (function 'metabin')", sep=""))
        }
      }
      ##
      ##
      bias <- lreg$slope
      df <- lreg$df
      slope <- lreg$intercept
      ##
      statistic <- bias / se.bias
      p.value <- 2*(1 - pt(abs(statistic), df=df))[[1]]
      ##
      res <- list(estimate=c(bias, se.bias, slope),
                  parameters=df, statistic=statistic,
                  p.value=p.value)
      
      names(res$statistic) <- "t"
      names(res$parameters) <- "df"
      names(res$estimate) <- c("bias", "se.bias", "slope")
    }
    else if (method.bias=="count"){
      ##
      ## Schwarzer, Antes, Schumacher (2006), Stat Med
      ##
      if (inherits(x, "metabin")){
        TE.MH <- metabin(x$event.e, x$n.e,
                         x$event.c, x$n.c,
                         sm="OR", method="MH", warn=FALSE)$TE.fixed
        
        n.. <- n.e + n.c
        n11 <- event.e
        n1. <- n.e
        n.1 <- event.e + event.c
        
        n11.e <- rep(NA, k)
        n11.v <- rep(NA, k)
        ##
        for ( i in seq(1:k) ){
          obj <- hypergeometric(n1.[i], n.1[i], n..[i], exp(TE.MH))
          n11.e[i] <- obj$mean()
          n11.v[i] <- obj$var()
        }
        
        ktau <- kentau((n11-n11.e)/sqrt(n11.v), 1/n11.v, correct=correct)
        ##
        res <- list(estimate=c(ktau$ks, ktau$se.ks),
                    statistic=ktau$ks/ktau$se.ks,
                    p.value=ktau$p.value)
        names(res$statistic) <- "z"
        names(res$estimate) <- c("ks", "se.ks")
      }
      else{
          stop(paste("method.bias '", method.bias, "' only defined for meta-analysis with binary outcome data (function 'metabin')", sep=""))
      }
    }
  
    res$alternative <- "asymmetry in funnel plot"
    
    res$method <- c(paste("Rank correlation test of funnel plot asymmetry",
                          ifelse(correct==TRUE, " (with continuity correction)", ""),
                          sep=""),
                    "Linear regression test of funnel plot asymmetry",
                    "Linear regression test of funnel plot asymmetry (methods of moment)",
                    paste("Rank correlation test of funnel plot asymmetry (based on counts)",
                          ifelse(correct==TRUE, " (with continuity correction)", ""),
                          sep=""),
                    "Linear regression test of funnel plot asymmetry (efficient score)",
                    "Linear regression test of funnel plot asymmetry (based on sample size)")[imeth]
    
    res$data.name <- data.name
    
    if (plotit){
      ##
      if (method.bias=="linreg"|method.bias=="mm"){
        radial(TE, seTE, comb.fixed=FALSE)
        abline(lreg$slope, lreg$intercept)
      }
      else if (method.bias=="rank"){
        ##
        if (plotit){
          plot(TE.s, seTE^2,
               xlab="Standardised treatment effect",
               ylab="Estimated variance of treatment estimate")
        }
      }
      else if (method.bias=="score"){
        ##
        radial(TE.score, seTE.score, comb.fixed=FALSE)
        abline(lreg$slope, lreg$intercept)
      }
    }

    if (inherits(x, "meta")){
      if (!is.null(x$title))
        res$title <- x$title
      if (!is.null(x$complab))
        res$complab <- x$complab
      if (!is.null(x$outclab))
        res$outclab <- x$outclab
    }
  }
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metabias", "htest")
  
  res
}
