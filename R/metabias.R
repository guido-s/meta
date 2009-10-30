metabias <- function(x, seTE, TE.fixed, seTE.fixed,
                     method="rank",
                     plotit=FALSE, correct=FALSE){
  
  
  kentau <- function(x, y, correct=FALSE, keep.data=FALSE){
    ##
    ## Check:
    ##
    if(length(x) != length(y))
      stop("length of argument x and y must be equal")
    ##
    sel <- !is.na(x) & !is.na(y)
    if (length(x) != sum(sel))
      warning(paste(length(x) - sum(sel),
                    "observation(s) dropped due to missing values"))
    ##
    x <- x[sel]
    y <- y[sel]
    n <- length(x)
    ##
    ks <- .C("kenscore",
             kenscore = as.double(0),
             x = as.double(x),
             y = as.double(y),
             n = as.integer(n),
             PACKAGE="meta")$kenscore
    ##
    ## Calculate S and s.e(S) according to
    ## Stata, release 5, Reference P-Z, p.239-240
    ##
    ## see also Kendall, Gibbons (1990), Rank Correlation Methods
    ## p. 66-68
    ##
    t <- rle(sort(x))$lengths
    u <- rle(sort(y))$lengths
    ##
    ##
    se.ks <- sqrt(1/18 * (n*(n-1)*(2*n+5) -
                          sum(t*(t-1)*(2*t+5)) -
                          sum(u*(u-1)*(2*u+5))) +
                  1/(9*n*(n-1)*(n-2)) *
                  sum(t*(t-1)*(t-2))*
                  sum(u*(u-1)*(u-2)) +
                  1/(2*n*(n-1)) *
                  sum(t*(t-1)) *
                  sum(u*(u-1)))
    ##
    if (as.logical(correct) &
        any(c(length(unique(x)), length(unique(y)))==2))
      warning(paste("Continuity corrected statistic may be inappropriate,\n",
                    "see Kendall, Gibbons (1990), Rank Correlation Methods, ",
                    "p.67\n ", sep=""))
    ##
    statistic <- (ks - sign(ks) * as.logical(correct)) / se.ks
    p.value <- 2*(1-pnorm(abs(statistic)))
    ##
    N <- 0.5 * n * (n-1)
    N1 <- N-sum(t*(t-1)/2)
    N2 <- N-sum(u*(u-1)/2)
    ##
    res <- list(tau.a=ks/N,
                tau.b=ks/(sqrt(N1)*sqrt(N2)),
                ks=ks - sign(ks) * as.logical(correct),
                se.ks=se.ks,
                statistic=statistic,
                p.value=p.value,
                correction=as.logical(correct))
    ##
    if (keep.data){
      res$x <- x
      res$y <- y
    }
    ##
    res
  }
  
  
  hypergeometric <- function(n1, m1, N, psi){    
    ##
    ## R program for computing the mean, variance, density, cumulative
    ## distribution and generating random deviates.
    ##
    ## Based on
    ## Liao and Rosen (2001): Fast and Stable Algorithms for Computing
    ## and Sampling from the Noncentral Hypergeometric Distribution,
    ## The American Statistician, 55, 236-369. 
    ##
    ## this is how to use the function
    ##
    ## n1 <- 100
    ## n2 <- 100
    ## m1 <- 100
    ## N <- n1+n2
    ## odds.ratio <- 3;
    ## obj <- hypergeometric(n1, m1, N, odds.ratio)
    ## obj$mean()
    ## obj$var()
    ## obj$d(40)
    ## obj$p(40)
    ## obj$r()
    
    n2 <- N - n1;
    
    if(n1<0 | n2<0 | m1<0 | m1>N | psi<=0)
      stop("wrong argument in hypergeometric");
    
    mode.compute <- function(){
      a <- psi - 1;
      b <- -( (m1+n1+2)*psi + n2-m1 ) ;    
      c <- psi*(n1+1)*(m1+1);
      q <- b + sign(b)*sqrt(b*b-4*a*c);
      q <- -q/2;
      
      mode <- trunc(c/q); 
      if(uu>=mode && mode>=ll) return(mode)
      else return( trunc(q/a) );      
    }       
    
    r.function <- function(i) (n1-i+1)*(m1-i+1)/i/(n2-m1+i)*psi;
    
    ##
    mean <- function() sum( prob[(ll:uu)+shift]*(ll:uu) ); 
    
    var <-  function() sum( prob[(ll:uu)+shift]*(ll:uu)^2 ) - mean()^2;          
    
    d <- function(x) return(prob[x + shift]);
    
    p <- function(x, lower.tail=TRUE){   
      if(lower.tail) return( sum(prob[ll:(x+shift)]) )
      else return( sum( prob[(x+shift):uu] ) );
    }
    
    ##
    
    sample.low.to.high <- function(lower.end, ran){ 
      for(i in lower.end:uu){                                
        if(ran <= prob[i+shift]) return(i);
        ran <- ran - prob[i+shift];
      }                                
    }
    
    sample.high.to.low <- function(upper.end, ran){           
      for(i in upper.end:ll){                              
        if(ran <= prob[i+shift]) return(i);
        ran <- ran - prob[i+shift];
      } 
    }  
    
    
    r <- function(){
      ran <- runif(1); 
      
      if(mode==ll) return( sample.low.to.high(ll, ran) );            
      if(mode==uu) return( sample.high.to.low(uu, ran) );                                         
      
      if(ran < prob[mode+shift]) return(mode);             
      ran <- ran - prob[mode+shift];
      
      lower <- mode - 1;                                                                            
      upper <- mode + 1;
      
      repeat{                                     
        if(prob[upper + shift] >= prob[lower + shift]){              
          if(ran < prob[upper+shift]) return(upper);
          ran <- ran - prob[upper+shift];
          if(upper==uu) return( sample.high.to.low(lower, ran) );
          upper <- upper + 1;                            
        }
        
        else{
          if(ran < prob[lower+shift]) return(lower);
          ran <- ran - prob[lower+shift];
          if(lower==ll) return( sample.low.to.high(upper, ran) );
          lower <- lower - 1;                   
        }      
      } 
    }
    
    ##
    
    ll <- max(0, m1-n2);
    uu <- min(n1, m1);              
    mode <- mode.compute();  
    
    prob <- array( 1, uu-ll+1 );
    
    shift <- 1-ll; 
    if(mode<uu) #note the shift of location
      {  
        r1 <- r.function( (mode+1):uu );       
        prob[ (mode+1 + shift):(uu + shift) ] <- cumprod(r1);       
      }
    
    if(mode>ll){
      r1 <- 1/r.function( mode:(ll+1) );
      prob[ (mode-1 + shift):(ll + shift) ] <- cumprod(r1);
    }
    
    prob <- prob/sum(prob); 
    
    return(list(mean=mean, var=var, d=d, p=p, r=r));
  }
  
  
  linregcore <- function(x, y, w=NULL){
    ##
    ## core function for method linreg and mm
    ##
    ##
    ## Check:
    ##
    if(length(x) != length(y))
      stop("length of argument x and y must be equal")
    ##
    if (!is.null(w)){
      if(length(x) != length(w))
        stop("length of argument x and w must be equal")
    }
    ##
    sel <- !is.na(x) & !is.na(y)
    if (length(x) != sum(sel))
      warning(paste(length(x) - sum(sel),
                    "observation(s) dropped due to missing values"))
    ##
    x <- x[sel]
    y <- y[sel]
    n <- length(x)
    ##
    if (is.null(w))
      w <- rep(1, n)
    else
      w <- w[sel]
    ##
    W <- diag(w)
    ##
    if (n > 2){
      X <- cbind(rep(1, n), x)
      XWX <- solve(t(X) %*% W %*% X)
      ##
      coefs <- XWX %*% t(X) %*% W %*% y
      MSE.w <- 1/(n-2) * sum(w*(y - X %*% coefs)^2)
      df <- n - 2
      se.coefs <- sqrt(MSE.w * diag(XWX))
      ##
      intercept <- coefs[1]
      se.intercept <- se.coefs[1]
      slope <- coefs[2]
      se.slope <- se.coefs[2]
    }
    else {
      intercept <- NA
      se.intercept <- NA
      slope <- NA
      se.slope <- NA
    }
    ##
    res <- list(intercept=intercept,
                se.intercept=se.intercept,
                slope=slope,
                se.slope=se.slope,
                df=df,
                MSE.w=MSE.w)
    res
  }
  
  
  if (inherits(x, "meta")){
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
      if (x$meth=="Peto"){
        warning(paste("Inverse variance method used for pooling",
                      "to perform test of funnel plot asymmetry"))
        m <- metabin(event.e, n.e, event.c, n.c, sm=x$sm, meth="Inv",
                     incr=x$incr, allincr=x$allincr,
                     allstudies=x$allstudies, MH.exact=x$MH.exact,
                     RR.cochrane=x$RR.cochrane)
        TE <- m$TE
        seTE <- m$seTE
      }
    }
  }
  else{
    TE <- x
    data.name <- paste(deparse(substitute(x)),
                       deparse(substitute(seTE)),
                       sep=", ")
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
  
  
  imeth <- charmatch(method,
                     c("rank", "linreg", "mm", "count", "score", "peters"),
                     nomatch = NA)
  if(is.na(imeth) | imeth==0)
    stop("method should be \"rank\", \"linreg\", \"mm\", \"count\", \"score\", or \"peters\"")
  ##
  method <- c("rank", "linreg", "mm", "count", "score", "peters")[imeth]
  
  
  if (k <= 2){
    warning("number of trials is too small")
    res <- NULL
  }
  else{
    if (method=="rank"){
      ##
      ## Begg und Mazumdar (1994), Biometrics, 50, 1088-1101
      ##
      m <- metagen(TE, seTE)
      TE.fixed <- m$TE.fixed
      seTE.fixed <- m$seTE.fixed
      ##
      varTE <- seTE^2
      varTE.s <- varTE - seTE.fixed^2
      TE.s <- (TE - TE.fixed) / sqrt(varTE.s)
      ##
      ktau <- kentau(TE.s, varTE, correct=correct)
      ##
      res <- list(estimate=c(ktau$ks, ktau$se.ks),
                  statistic=ktau$ks/ktau$se.ks,
                  p.value=ktau$p.value)
      names(res$statistic) <- "z"
      names(res$estimate) <- c("ks", "se.ks")
    }
    else if (method=="linreg" | method=="mm" | method=="score" |
             method=="peters"){
      
      if (method=="linreg"){
        ##
        ## Egger, Smith, Schneider, Minder (1997), BMJ, 315, 629-34
        ##
        lreg <- linregcore(seTE, TE, 1/seTE^2)
        se.bias <- lreg$se.slope
      }
      else if (method=="mm"){
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
      else if (method=="score"){
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
          stop(paste("method '", method, "' only defined for meta-analysis with binary outcome data (function 'metabin')", sep=""))
        }
      }
      else if (method=="peters"){
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
          stop(paste("method '", method, "' only defined for meta-analysis with binary outcome data (function 'metabin')", sep=""))
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
    else if (method=="count"){
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
          stop(paste("method '", method, "' only defined for meta-analysis with binary outcome data (function 'metabin')", sep=""))
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
      if (method=="linreg"|method=="mm"){
        radial(TE, seTE, comb.fixed=FALSE)
        abline(lreg$slope, lreg$intercept)
      }
      else if (method=="rank"){
        ##
        if (plotit){
          plot(TE.s, varTE,
               xlab="Standardised treatment effect",
               ylab="Estimated variance of treatment estimate")
        }
      }
      else if (method=="score"){
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
    
    class(res) <- c("metabias", "htest")
  }
  
  res
}
