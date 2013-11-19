xlab <- function(sm, pscale=1){
  if      (sm=="OR" ) res <- "Odds Ratio"
  else if (sm=="RR" ) res <- "Risk Ratio"
  else if (sm=="RD" ) res <- "Risk Difference"
  else if (sm=="AS" ) res <- "Arcus Sinus Transformation"
  ##
  else if (sm=="SMD") res <- "Standardised mean difference"
  else if (sm=="WMD"|sm=="MD") res <- "Mean difference"
  ##
  else if (sm=="HR" ) res <- "Hazard Ratio"
  ##
  else if (sm=="IRR" ) res <- "Incidence Rate Ratio"
  else if (sm=="IRD" ) res <- "Incidence Rate Difference"
  ##
  else if (sm=="ZCOR") res <- "Correlation (based on Fisher's z transformation)"
  else if (sm=="COR") res <- "Correlation"
  ##else if (sm=="COR") res <- "Correlation (untransformed)"
  ##
  else if (sm %in% c("PFT", "PAS", "PRAW", "PLN", "PLOGIT") & pscale==100)
    res <- "Proportion (in %)"
  else if (sm %in% c("PFT", "PAS", "PRAW", "PLN", "PLOGIT"))
    res <- "Proportion"
  ##
  else res <- sm
  ##
  res
}
