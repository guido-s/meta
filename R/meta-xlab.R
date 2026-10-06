# Auxiliary function
#
# Package: meta
# Author: Guido Schwarzer <guido.schwarzer@uniklinik-freiburg.de>
# License: GPL (>= 2)
#

xlab_meta <- function(sm, backtransf,
                      pscale = 1, irscale = 1, irunit = "person-years",
                      newline = FALSE, revman5 = FALSE,
                      big.mark = gs("big.mark"),
                      func.transf = NULL, func.backtransf = NULL) {
  
  res <- sm
  
  newline <- if (newline) "\n" else " "
  
  #
  # metabin() - gs("sm4bin")
  #
  if (sm == "RD")
    if (pscale == 1)
      res <- "Risk Difference"
  else
    res <- paste0("Risk Difference\n(events per ",
                  format(pscale, scientific = FALSE, big.mark = big.mark),
                  " obs.)")
  #
  else if (sm == "ASD")
    res <- paste0("Arcus Sinus Difference")
  #
  # metacont() - gs("sm4cont")
  #
  else if (sm == "MD" || sm == "WMD")
    res <- "Mean Difference"
  #
  else if (sm == "SMD")
    res <- paste0(if (revman5) "Std. Mean" else "Standardised Mean",
                  newline, "Difference")
  #
  # metacor() - gs("sm4cor")
  #
  else if (sm == "COR")
    res <- "Correlation"
  #
  # metainc() - gs("sm4inc")
  #
  else if (sm == "IRD")
    if (irscale == 1)
      res <- paste0("Incidence Rate Difference")
  else
    res <- paste0("Incidence Rate Difference\n(events per ",
                  format(irscale, scientific = FALSE, big.mark = big.mark),
                  newline, irunit)
  #
  # metamean() - gs("sm4mean")
  #
  else if (sm == "MRAW")
    res <- "Mean"
  #
  else if (backtransf) {
    #
    # metabin() - gs("sm4bin")
    #
    if (sm == "OR")
      res <- "Odds Ratio"
    #
    else if (sm == "RR")
      res <- "Risk Ratio"
    #
    else if (sm == "DOR")
      res <- "Diagnostic Odds Ratio"
    #
    else if (sm == "VE")
      res <- "Vaccine Eff."
    #
    # metacont() - gs("sm4cont")
    #
    else if (sm == "ROM")
      res <- "Ratio of Means"
    #
    # metacor() - gs("sm4cor")
    #
    else if (sm == "ZCOR")
      res <- "Correlation"
    #
    # metamean() - gs("sm4mean")
    #
    else if (is_mean(sm))
      res <- "Mean"
    #
    # metainc() - gs("sm4inc")
    #
    else if (sm == "IRR")
      res <- "Incidence Rate Ratio"
    #
    else if (sm == "IRSD")
      res <- "Incidence Rate Difference"
    #
    # metaprop() - gs("sm4prop")
    #
    else if (is_prop(sm)) {
      if (pscale == 1)
        res <- "Proportion"
      else
        res <- paste0("Events per ",
                      format(pscale, scientific = FALSE, big.mark = big.mark),
                      newline, "observations")
    }
    #
    # metarate() - gs("sm4rate")
    #
    else if (is_rate(sm)) {
      if (irscale == 1)
        res <- "Incidence Rate"
      else
        res <- paste0("Events per ",
                      format(irscale, scientific = FALSE, big.mark = big.mark),
                      newline, irunit)
    }
    #
    # metagen()
    #
    else if (sm == "HR")
      res <- "Hazard Ratio"
  }
  else {
    #
    # metabin() - gs("sm4bin")
    #
    if (sm == "OR")
      res <- "Log Odds Ratio"
    #
    else if (sm == "RR")
      res <- "Log Risk Ratio"
    #
    else if (sm == "DOR")
      res <- "Log Diagnostic Odds Ratio"
    #
    else if (sm == "VE")
      res <- "Log Vaccine Ratio"
    #
    # metacont() - gs("sm4cont")
    #
    else if (sm == "ROM")
      res <- "Log Ratio of Means"
    #
    # metacor() - gs("sm4cor")
    #
    else if (sm == "ZCOR")
      res <- paste0("Fisher's z transformed", newline, "Correlation")
    #
    # metamean() - gs("sm4mean")
    #
    else if (sm == "MLN")
      res <- "Log Mean"
    #
    # metainc() - gs("sm4inc")
    #
    else if (sm == "IRR")
      res <- "Log Incidence Rate Ratio"
    #
    else if (sm == "IRSD")
      res <-  paste0("Square Root of", newline, "Incidence Rate Difference")
    #
    # metaprop() - gs("sm4prop")
    #
    else if (sm == "PLOGIT")
      res <- paste0("Logit Transformed", newline, "Proportion")
    #
    else if (sm == "PLN")
      res <- paste0("Log Transformed", newline, "Proportion")
    #
    else if (sm == "PRAW")
      res <- paste0("Untransformed Proportion")
    #
    else if (sm == "PAS")
      res <- paste0("Arcsine Transformed", newline, "Proportion")
    #
    else if (sm == "PFT")
      res <- paste0("Freeman-Tukey Double Arcsine", newline,
                    "Transformed Proportion")
    #
    # metarate() - gs("sm4rate")
    #
    else if (sm == "IRLN")
      res <- paste0("Log Incidence Rate")
    #
    else if (sm == "IRS")
      res <- paste0("Square Root of", newline, "Incidence Rate")
    #
    else if (sm == "IRFT")
      res <- paste0("Freeman-Tukey Double Arcsine", newline,
                   "Transformed Rate")
    #
    # metagen()
    #
    else if (sm == "HR")
      res <- "Log Hazard Ratio"
    #
    else if (!is.null(func.transf))
      res <- paste0(func.transf, "(", sm, ")")
    else if (!is.null(func.backtransf)) {
      if (func.backtransf == "exp")
        res <- paste0("log(", sm, ")")
      else if (func.backtransf == "z2cor")
        res <-  paste0("Fisher's z transformed", newline, "Correlation")
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
