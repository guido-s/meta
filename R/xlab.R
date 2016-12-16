xlab <- function(sm, backtransf,
                 pscale = 1, irscale = 1, irunit = "person-years",
                 newline = FALSE) {
  
  res <- NULL
  
  
  newline <- if (newline) "\n" else " "
  
  
  if (sm == "SMD")
    res <- paste("Standardised mean", newline,
                 "difference", sep = "")
  ##
  else if (sm == "WMD" | sm == "MD")
    res <- "Mean difference"
  ##
  else if (sm == "COR")
    res <- "Correlation"
  ##
  else if (sm == "RD")
    res <- "Risk Difference"
  ##
  else if (sm == "ASD")
    res <- paste("Arcus Sinus", newline,
                 "Difference", sep = "")
  ##
  else if (sm == "IRD")
    res <- paste("Incidence Rate", newline,
                 "Difference", sep = "")
  ##
  else if (sm == "IR")
    res <- "Incidence Rate"
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
    else if (sm %in% c("PFT", "PAS", "PLN", "PLOGIT", "PRAW")) {
      if (pscale == 1)
        res <- ""
      else
        res <- paste("Events per ",
                     format(pscale, scientific = FALSE),
                     newline,
                     "observations", sep = "")
    }
    ##
    else if (sm %in% c("IR", "IRLN", "IRS", "IRFT")) {
      if (irscale == 1)
        res <- "Incidence Rate"
      else
        res <- paste("Events per ",
                     format(irscale, scientific = FALSE),
                     newline,
                     irunit, sep = "")
    }
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
  }
  
  
  if (is.null(res))
    res <- sm

  
  res
}
