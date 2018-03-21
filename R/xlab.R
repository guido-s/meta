xlab <- function(sm, backtransf,
                 pscale = 1, irscale = 1, irunit = "person-years",
                 newline = FALSE, revman5 = FALSE,
                 big.mark = big.mark) {
  
  res <- NULL
  
  
  newline <- if (newline) "\n" else " "
  
  
  if (sm == "SMD")
    res <- paste(if (revman5) "Std. Mean" else "Standardised Mean",
                 newline,
                 "Difference", sep = "")
  ##
  else if (sm == "WMD" | sm == "MD")
    res <- "Mean Difference"
  ##
  else if (sm == "COR")
    res <- "Correlation"
  ##
  else if (sm == "RD")
    if (pscale == 1)
      res <- "Risk Difference"
    else
      res <- paste("Risk Difference\n(events per ",
                   format(pscale, scientific = FALSE, big.mark = big.mark),
                   " obs.)",
                   sep = "")
  ##
  else if (sm == "ASD")
    res <- paste("Arcus Sinus", newline,
                 "Difference", sep = "")
  ##
  else if (sm == "IRD")
    if (irscale == 1)
      res <- paste("Incidence Rate", newline,
                   "Difference", sep = "")
    else
      res <- paste("Incidence Rate Diff.\n(events per ",
                   format(irscale, scientific = FALSE, big.mark = big.mark),
                   newline,
                   irunit, sep = "")
  ##
  else if (sm == "IR")
    res <- "Incidence Rate"
  ##
  else if (sm == "MRAW")
    res <- "Mean"
  ##
  else if (backtransf) {
    if (sm == "OR")
      res <- "Odds Ratio"
    ##
    else if (sm == "RR")
      res <- "Risk Ratio"
    ##
    else if (sm == "ROM")
      res <- "Ratio of Means"
    ##
    else if (sm == "ZCOR")
      res <- "Correlation"
    ##
    else if (sm == "HR")
      res <- "Hazard Ratio"
    ##
    else if (sm == "IRR")
      res <- paste("Incidence Rate", newline,
                   "Ratio", sep = "")
    ##
    else if (is.prop(sm)) {
      if (pscale == 1)
        res <- ""
      else
        res <- paste("Events per ",
                     format(pscale, scientific = FALSE, big.mark = big.mark),
                     newline,
                     "observations", sep = "")
    }
    ##
    else if (is.rate(sm)) {
      if (irscale == 1)
        res <- "Incidence Rate"
      else
        res <- paste("Events per ",
                     format(irscale, scientific = FALSE, big.mark = big.mark),
                     newline,
                     irunit, sep = "")
    }
    ##
    else if (sm == "MLN")
      res <- "Mean"
  }
  else {
    if (sm == "OR")
      res <- "Log Odds Ratio"
    ##
    else if (sm == "RR")
      res <- "Log Risk Ratio"
    ##
    else if (sm == "ROM")
      res <- paste("Log Ratio of", newline,
                   "Means", sep = "")
    ##
    else if (sm == "HR")
      res <- paste("Log Hazard", newline,
                   "Ratio", sep = "")
    ##
    else if (sm == "IRR")
      res <- paste("Log Incidence Rate", newline,
                   "Ratio", sep = "")
    ##
    else if (sm == "ZCOR")
      res <- paste("Fisher's z transformed", newline,
                   "correlation", sep = "")
    ##
    else if (sm == "PFT")
      res <- paste("Freeman-Tukey Double Arcsine", newline,
                   "Transformed Proportion", sep = "")
    ##
    else if (sm == "PAS")
      res <- paste("Arcsine Transformed", newline,
                   "Proportion", sep = "")
    ##
    else if (sm == "PLN")
      res <- paste("Log Transformed", newline,
                   "Proportion", sep = "")
    ##
    else if (sm == "PLOGIT")
      res <- paste("Logit Transformed", newline,
                   "Proportion", sep = "")
    ##
    else if (sm == "PRAW")
      res <- paste("Untransformed", newline,
                   "Proportion", sep = "")
    ##
    else if (sm == "IR")
      res <- "Incidence Rate"
    ##
    else if (sm == "IRLN")
      res <- paste("Log Incidence", newline,
                   "Rate", sep = "")
    ##
    else if (sm == "IRS")
      res <- paste("Square Root of", newline,
                   "Incidence Rate", sep = "")
    ##
    else if (sm == "IRFT")
      res <- paste("Freeman-Tukey Double Arcsine", newline,
                   "Transformed Rate", sep = "")
    ##
    else if (sm == "MLN")
      res <- "Log Mean"
  }
  
  
  if (is.null(res))
    res <- sm

  
  res
}
