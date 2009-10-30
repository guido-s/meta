nnt.metabin <- function(x, baseline.risk){
  
  
  res <- list()
  
  if (x$sm=="RD"){
    res$fixed  <- meta::ci(x$TE.fixed, x$seTE.fixed)[c(1,4,3)]
    res$fixed$TE <- -1/res$fixed$TE
    res$fixed$lower <- -1/res$fixed$lower
    res$fixed$upper <- -1/res$fixed$upper
    ##
    res$random <- meta::ci(x$TE.random, x$seTE.random)[c(1,4,3)]
    res$random$TE <- -1/res$random$TE
    res$random$lower <- -1/res$random$lower
    res$random$upper <- -1/res$random$upper
  }
  ##
  if (x$sm=="RR"){
    res$fixed <- meta::ci(x$TE.fixed, x$seTE.fixed)[c(1,3,4)]
    ##
    res$fixed$TE    <- -1/(baseline.risk*exp(res$fixed$TE)-baseline.risk)
    res$fixed$lower <- -1/(baseline.risk*exp(res$fixed$lower)-baseline.risk)
    res$fixed$upper <- -1/(baseline.risk*exp(res$fixed$upper)-baseline.risk)
    ##
    res$random <- meta::ci(x$TE.random, x$seTE.random)[c(1,3,4)]
    ##
    res$random$TE    <- -1/(baseline.risk*exp(res$random$TE)-baseline.risk)
    res$random$lower <- -1/(baseline.risk*exp(res$random$lower)-baseline.risk)
    res$random$upper <- -1/(baseline.risk*exp(res$random$upper)-baseline.risk)
  }
  
  res  
}
