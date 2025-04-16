xlab_meta <- function(sm, backtransf,
                      pscale = 1, irscale = 1, irunit = "person-years",
                      newline = FALSE, revman5 = FALSE,
                      big.mark = gs("big.mark"),
                      func.transf = NULL, func.backtransf = NULL) {
  
  res <- sm
  
  
  newline <- if (newline) "\n" else " "
  
  
  if (sm == "SMD")
    res <- paste0(if (revman5) "Std. Mean" else "Standardised Mean",
                  newline, "Difference")
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
      res <- paste0("Risk Difference\n(events per ",
                    format(pscale, scientific = FALSE, big.mark = big.mark),
                    " obs.)")
  ##
  else if (sm == "ASD")
    res <- paste0("Arcus Sinus", newline, "Difference")
  ##
  else if (sm == "IRD")
    if (irscale == 1)
      res <- paste0("Incidence Rate", newline, "Difference")
    else
      res <- paste0("Incidence Rate Diff.\n(events per ",
                    format(irscale, scientific = FALSE, big.mark = big.mark),
                    newline, irunit)
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
    else if (sm == "DOR")
      res <- "Diagnostic Odds Ratio"
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
      res <- paste0("Incidence Rate", newline, "Ratio")
    ##
    else if (sm == "VE")
      res <- "Vaccine Eff."
    ##
    else if (is_prop(sm)) {
      if (pscale == 1)
        res <- ""
      else
        res <- paste0("Events per ",
                      format(pscale, scientific = FALSE, big.mark = big.mark),
                      newline, "observations")
    }
    ##
    else if (is_rate(sm)) {
      if (irscale == 1)
        res <- "Incidence Rate"
      else
        res <- paste0("Events per ",
                      format(irscale, scientific = FALSE, big.mark = big.mark),
                      newline, irunit)
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
    else if (sm == "DOR")
      res <- "Log Diagnostic Odds Ratio"
    ##
    else if (sm == "ROM")
      res <- paste0("Log Ratio of", newline, "Means")
    ##
    else if (sm == "HR")
      res <- paste0("Log Hazard", newline, "Ratio")
    ##
    else if (sm == "IRR")
      res <- paste0("Log Incidence Rate", newline, "Ratio")
    ##
    else if (sm == "ZCOR")
      res <- paste0("Fisher's z transformed", newline, "correlation")
    ##
    else if (sm == "VE")
      res <- "Log Risk Ratio"
    ##
    else if (sm == "PFT")
      res <- paste0("Freeman-Tukey Double Arcsine", newline,
                    "Transformed Proportion")
    ##
    else if (sm == "PAS")
      res <- paste0("Arcsine Transformed", newline, "Proportion")
    ##
    else if (sm == "PLN")
      res <- paste0("Log Transformed", newline, "Proportion")
    ##
    else if (sm == "PLOGIT")
      res <- paste0("Logit Transformed", newline, "Proportion")
    ##
    else if (sm == "PRAW")
      res <- paste0("Untransformed", newline, "Proportion")
    ##
    else if (sm == "IR")
      res <- "Incidence Rate"
    ##
    else if (sm == "IRLN")
      res <- paste0("Log Incidence", newline, "Rate")
    ##
    else if (sm == "IRS")
      res <- paste0("Square Root of", newline, "Incidence Rate")
    ##
    else if (sm == "IRFT")
      res <- paste0("Freeman-Tukey Double Arcsine", newline,
                   "Transformed Rate")
    ##
    else if (sm == "MLN")
      res <- "Log Mean"
    else if (!is.null(func.transf))
      res <- paste0(func.transf, "(", sm, ")")
    else if (!is.null(func.backtransf)) {
      if (func.backtransf == "exp")
        res <- paste0("log(", sm, ")")
      else if (func.backtransf == "z2cor")
        res <-  paste0("Fisher's z transformed", newline, "correlation")
      else if (func.backtransf == "logit2p")
        res <-  paste0("Logit Transformed", newline, "Proportion")
      else if (func.backtransf == "logVR2VE")
        res <-  "Log Vaccine Ratio"
    }
  }
  
  if (is.null(res))
    res <- sm

  
  res
}
