xlab <- function(sm, backtransf){

  res <- NULL
  
  if (sm=="SMD")
    res <- "Standardised mean difference"
  ##
  else if (sm=="WMD"|sm=="MD")
    res <- "Mean difference"
  ##
  else if (sm=="COR")
    res <- "Correlation"
  ##
  else if (sm=="RD" )
    res <- "Risk Difference"
  ##
  else if (sm=="ASD" )
    res <- "Arcus Sinus Difference"
  ##
  else if (sm=="IRD" )
    res <- "Incidence Rate Difference"
  ##
  else if (backtransf){
    if (sm=="OR" )
      res <- "Odds Ratio"
    else if (sm=="RR" )
      res <- "Risk Ratio"
    ##
    else if (sm=="HR" )
      res <- "Hazard Ratio"
    ##
    else if (sm=="IRR" )
      res <- "Incidence Rate Ratio"
    ##
    else if (sm=="ZCOR")
      res <- "Correlation"
    ##
    else if (sm %in% c("PFT", "PAS", "PLN", "PLOGIT", "PRAW"))
      res <- ""
  }
  else{
    if (sm=="OR" )
      res <- "Log Odds Ratio"
    ##
    else if (sm=="RR" )
      res <- "Log Risk Ratio"
    ##
    else if (sm=="HR" )
      res <- "Log Hazard Ratio"
    ##
    else if (sm=="IRR" )
      res <- "Log Incidence Rate Ratio"
    ##
    else if (sm=="ZCOR")
      res <- "Fisher's z transformed correlation"
    ##
    else if (sm=="PFT")
      res <- "Freeman-Tukey Double Arcsine Transformed Proportion"
    ##
    else if (sm=="PAS")
      res <- "Arcsine Transformed Proportion"
    ##
    else if (sm=="PLN")
      res <- "Log Transformed Proportion"
    ##
    else if (sm=="PLOGIT")
      res <- "Logit Transformed Proportion"
    ##
    else if (sm=="PRAW")
      res <- "Untransformed Proportion"
  }
  
  
  if (is.null(res))
    res <- sm

  
  res
}
