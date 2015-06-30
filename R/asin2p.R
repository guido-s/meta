asin2p <- function(x, n=NULL, value="mean", warn=TRUE){
  
  ##
  ## Do nothing if all values are NA
  ## 
  if (all(is.na(x)))
    return(x)
  
  
  ##
  ## Calculate possible minimum and maximum
  ## for each transformation
  ##
  if (is.null(n)){
    minimum <- asin(sqrt(0))
    maximum <- asin(sqrt(1))
  }
  else{
    minimum <- asin(sqrt(0/(n+1))) + asin(sqrt((0+1)/(n+1)))
    maximum <- asin(sqrt(n/(n+1))) + asin(sqrt((n+1)/(n+1)))
  }
  ##
  sel0 <- x < minimum
  sel1 <- x > maximum
  
  
  ##
  ## Check for (impossible) negative values
  ##
  if (any(sel0, na.rm=TRUE)){
    if (is.null(n)){
      if (warn)
        warning("Negative value for ",
                if (length(x)>1) "at least one ",
                if (value=="mean") "transformed proportion using arcsine transformation.\n  Proportion set to 0.",
                if (value=="lower") "lower confidence limit using arcsine transformation.\n  Lower confidence limit set to 0.",
                if (value=="upper") "upper confidence limit using arcsine transformation.\n  Upper confidence limit set to 0.",
                sep="")
    }
    else{
      if (warn)
        warning("Too small value for ",
                if (length(x)>1) "at least one ",
                if (value=="mean") "transformed proportion using Freeman-Tukey double arcsine transformation.\n  Proportion set to 0.",
                if (value=="lower") "lower confidence limit using Freeman-Tukey double arcsine transformation.\n  Lower confidence limit set to 0.",
                if (value=="upper") "upper confidence limit using Freeman-Tukey double arcsine transformation.\n  Upper confidence limit set to 0.",
                sep="")
    }
  }
  
  ##
  ## Check for (impossible) large values
  ##
  if (any(sel1, na.rm=TRUE)){
    if (is.null(n)){
      if (warn)
        warning("Too large value for ",
                if (length(x)>1) "at least one ",
                if (value=="mean") "transformed proportion using arcsine transformation.\n  Proportion set to 1.",
                if (value=="lower") "lower confidence limit using arcsine transformation.\n  Lower confidence limit set to 1.",
                if (value=="upper") "upper confidence limit using arcsine transformation.\n  Upper confidence limit set to 1.",
                sep="")
    }
    else{
      if (warn)
        warning("Too large value for ",
                if (length(x)>1) "at least one ",
                if (value=="mean") "transformed proportion using Freeman-Tukey double arcsine transformation.\n  Proportion set to 1.",
                if (value=="lower") "lower confidence limit using Freeman-Tukey double arcsine transformation.\n  Lower confidence limit set to 1.",
                if (value=="upper") "upper confidence limit using Freeman-Tukey double arcsine transformation.\n  Upper confidence limit set to 1.",
                sep="")
    }
  }
  
  
  res <- rep(NA, length(x))
  ##
  sel <- !(sel0 | sel1)
  sel <- !is.na(sel) & sel
  ##
  res[sel0] <- 0
  res[sel1] <- 1
  ##
  if (is.null(n)){
    ##
    ## Back transformation of arcsine transformation:
    ##
    res[sel] <- sin(x[sel])^2
  }
  else{
    ##
    ## Back transformation of Freeman-Tukey double arcsine transformation:
    ##
    res[sel] <- 0.5*(1-sign(cos(x[sel]))*
                     sqrt(1-(sin(x[sel])+
                             (sin(x[sel])-1/sin(x[sel]))/n[sel])**2))
  }
  res
}
