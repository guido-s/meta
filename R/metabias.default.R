metabias.default <- function(x, seTE,
                             method.bias="linreg",
                             plotit=FALSE, correct=FALSE,
                             k.min=10, ...){
  
  TE <- x
  data.name <- paste(deparse(substitute(x)),
                     deparse(substitute(seTE)),
                     sep=", ")
  
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
  k <- length(TE)
  
  imeth <- charmatch(method.bias,
                     c("rank", "linreg", "mm"),
                     nomatch = NA)
  if(is.na(imeth) | imeth==0)
    stop("method.bias should be \"rank\", \"linreg\", or \"mm\"")
  ##
  method.bias <- c("rank", "linreg", "mm")[imeth]
  
  
  if (k < k.min){
    res <- list(k=k, k.min=k.min)
  }
  else{
    if (method.bias=="rank"){
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
    else if (method.bias=="linreg" | method.bias=="mm"){
      
      if (method.bias=="linreg"){
        ##
        ## Egger, Smith, Schneider, Minder (1997), BMJ, 315, 629-34
        ##
        lreg <- linregcore(seTE, TE, 1/seTE^2)
        se.bias <- lreg$se.slope
      }
      else if (method.bias=="mm"){
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
    
    res$alternative <- "asymmetry in funnel plot"
    
    res$method <- c(paste("Rank correlation test of funnel plot asymmetry",
                          ifelse(correct==TRUE, " (with continuity correction)", ""),
                          sep=""),
                    "Linear regression test of funnel plot asymmetry",
                    "Linear regression test of funnel plot asymmetry (methods of moment)")[imeth]
    
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
    }
  }
  
  res$version <- packageDescription("meta")$Version
  
  class(res) <- c("metabias", "htest")
  
  res
}
