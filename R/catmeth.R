catmeth <- function(method,
                    method.tau = NULL,
                    sm = "",
                    k.all,
                    hakn = FALSE,
                    class = "",
                    tau.common = FALSE,
                    tau.preset = NULL,
                    sparse = FALSE,
                    incr = NULL,
                    allincr = FALSE,
                    addincr = FALSE,
                    allstudies = FALSE,
                    doublezeros = FALSE,
                    MH.exact = FALSE,
                    method.ci = NULL,
                    pooledvar = FALSE,
                    method.smd,
                    sd.glass,
                    exact.smd = FALSE,
                    model.glmm,
                    pscale = 1,
                    irscale = 1,
                    irunit = "person-years",
                    null.effect = NA,
                    big.mark = "",
                    digits = gs("digits"),
                    digits.tau2 = gs("digits.tau2"),
                    text.tau2 = gs("text.tau2"),
                    method.miss, IMOR.e, IMOR.c
                    ) {
  
  metabin  <- "metabin"  %in% class
  metacont <- "metacont" %in% class
  metainc  <- "metainc"  %in% class
  metaprop <- "metaprop" %in% class
  metarate <- "metarate" %in% class
  trimfill <- "trimfill" %in% class
  metamiss <- "metamiss" %in% class
  
  if (is.null(allstudies)) allstudies <- FALSE
  if (is.null(doublezeros)) doublezeros <- FALSE
  
  if (sm == "ZCOR")
    sm.details <- "\n- Fisher's z transformation of correlations"
  else if (sm == "COR")
    sm.details <- "\n- Untransformed correlations"
  ##
  else if  (sm == "PFT" | sm == "IRFT")
    sm.details <- "\n- Freeman-Tukey double arcsine transformation"
  else if (sm == "PAS")
    sm.details <- "\n- Arcsine transformation"
  else if (is.log.effect(sm))
    sm.details <- "\n- Log transformation"
  else if (sm == "PLOGIT")
    sm.details <- "\n- Logit transformation"
  else if (sm == "PRAW")
    sm.details <- "\n- Untransformed proportions"
  ##
  else if  (sm == "IR")
    sm.details <- "\n- Untransformed rates"
  else if  (sm == "IRS")
    sm.details <- "\n- Square root transformation"
  ##
  else if  (sm == "MRAW")
    sm.details <- "\n- Untransformed (raw) means"
  ##
  else
    sm.details <- ""
  ##
  if (metaprop && !is.null(method.ci)) {
    if  (method.ci == "CP")
      method.ci.details <- "\n- Clopper-Pearson confidence interval for individual studies"
    else if (method.ci == "WS")
      method.ci.details <- "\n- Wilson Score confidence interval for individual studies"
    else if (method.ci == "WSCC")
      method.ci.details <- "\n- Wilson Score confidence interval with continuity correction for individual studies"
    else if (method.ci == "AC")
      method.ci.details <- "\n- Agresti-Coull confidence interval for individual studies"
    else if (method.ci == "SA")
      method.ci.details <- "\n- Simple approximation confidence interval for individual studies"
    else if (method.ci == "SACC")
      method.ci.details <- "\n- Simple approximation confidence interval with continuity correction for individual studies"
    else if (method.ci == "NAsm")
      method.ci.details <- "\n- Normal approximation confidence interval for individual studies"
    ##
    sm.details <- paste(sm.details, method.ci.details, sep = "")
  }
  ##
  if (metacont && sm == "SMD" && !is.null(method.smd)) {
    if  (method.smd == "Hedges")
      if (exact.smd)
        method.details <- "\n- Hedges' g (bias corrected standardised mean difference; using exact formulae)"
      else
        method.details <- "\n- Hedges' g (bias corrected standardised mean difference)"
    else if  (method.smd == "Cohen")
      if (exact.smd)
        method.details <- "\n- Cohen's d (standardised mean difference; using exact formulae)"
      else
        method.details <- "\n- Cohen's d (standardised mean difference)"
    else if  (method.smd == "Glass") {
      if (!is.null(sd.glass) && sd.glass == "control")
        method.details <- "\n- Glass' delta (standardised mean difference; based on control group)"
      else if (!is.null(sd.glass) && sd.glass == "experimental")
        method.details <- "\n- Glass' delta (standardised mean difference; based on experimental group)"
    }
    ##
    sm.details <- paste(sm.details, method.details, sep = "")
  }
  ##
  if (metabin | metainc | metaprop | metarate) {
    if (!(sm == "ASD" | method == "Peto")) {
      if (addincr) {
        if (all(incr == "TACC"))
          sm.details <- paste(sm.details,
                              "\n- Treatment arm continuity correction in all studies",
                              if (method == "GLMM") "\n  (only used to calculate individual study results)",
                              sep = "")
        else if (all(incr != 0))
          sm.details <- paste(sm.details,
                              "\n- Continuity correction",
                              if (length(unique(incr)) == 1)
                                paste(" of ", round(incr, 4), sep = ""),
                              " in all studies",
                              if (method == "GLMM") "\n  (only used to calculate individual study results)",
                              sep = "")
      }
      else if (sparse) {
        if (allincr == FALSE) {
          if (all(incr == "TACC"))
            sm.details <- paste(sm.details,
                                "\n- Treatment arm continuity correction in studies with zero cell frequencies",
                                if (method == "GLMM") "\n  (only used to calculate individual study results)",
                                sep = "")
          else if (any(incr != 0))
            sm.details <- paste(sm.details,
                                "\n- Continuity correction",
                                if (length(unique(incr)) == 1)
                                  paste(" of ", round(incr, 4), sep = ""),
                                " in studies with zero cell frequencies",
                                if (method == "GLMM") "\n  (only used to calculate individual study results)",
                                sep = "")
        }
        else {
          if (all(incr == "TACC"))
            sm.details <- paste(sm.details,
                                "\n- Treatment arm continuity correction in all studies",
                                if (method == "GLMM") "\n  (only used to calculate individual study results)",
                               sep = "")
          else if (any(incr != 0))
            sm.details <- paste(sm.details,
                                "\n- Continuity correction",
                                if (length(unique(incr)) == 1)
                                  paste(" of ", round(incr, 4), sep = ""),
                                " in all studies",
                                if (method == "GLMM") "\n  (only used to calculate individual study results)",
                                sep = "")
        }
        ##
        if (allstudies & doublezeros)
          sm.details <- paste(sm.details,
                              "\n- Studies with double zeros included in meta-analysis")
      }
    }
  }
  
  
  if (pscale != 1)
    sm.details <- paste(sm.details,
                        "\n- Events per ",
                        format(pscale, scientific = FALSE,
                               big.mark = big.mark),
                        " observations", sep = "")
  ##
  if (irscale != 1)
    sm.details <- paste(sm.details,
                        "\n- Events per ",
                        format(irscale, scientific = FALSE,
                               big.mark = big.mark),
                        " ", irunit, sep = "")

  if (!is.na(null.effect) && null.effect != 0) {
    if (pscale != 1)
      sm.details <- paste(sm.details,
                          "\n- Null hypothesis: effect is equal to ",
                          format(round(null.effect * pscale, digits),
                                 scientific = FALSE, big.mark = big.mark),
                          " events per ",
                          format(pscale, scientific = FALSE,
                                 big.mark = big.mark),
                          " observations", sep = "")
    else if (irscale != 1)
      sm.details <- paste(sm.details,
                          "\n- Null hypothesis: effect is equal to ",
                          format(round(null.effect * irscale, digits),
                                 scientific = FALSE, big.mark = big.mark),
                          " events per ",
                          format(irscale, scientific = FALSE,
                                 big.mark = big.mark),
                          " ", irunit, sep = "")
    else
      sm.details <- paste(sm.details,
                          "\n- Null hypothesis: effect is equal to ",
                          format(null.effect, scientific = FALSE,
                                 big.mark = big.mark),
                          sep = "")
  }
  
  
  lab.method.details <- ""
  ##
  if (is.null(method.tau))
    lab.method.tau <- ""
  else {
    if (!is.null(tau.preset)) {
      tau2 <- tau.preset^2
      tau2 <- formatPT(tau2, lab = TRUE, labval = text.tau2,
                       digits = digits.tau2,
                       lab.NA = "NA",
                       big.mark = big.mark)
      ##
      lab.method.tau <- paste("\n- Preset between-study variance: ",
                              tau2, sep = "")
      ##
      lab.method.details <- lab.method.tau
    }
    else {
      i.lab.method.tau <- charmatch(method.tau,
                                    c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB"),
                                    nomatch = NA)
      ##
      lab.method.tau <- c("\n- DerSimonian-Laird estimator for tau^2",
                          "\n- Paule-Mandel estimator for tau^2",
                          "\n- Restricted maximum-likelihood estimator for tau^2",
                          "\n- Maximum-likelihood estimator for tau^2",
                          "\n- Hunter-Schmidt estimator for tau^2",
                          "\n- Sidik-Jonkman estimator for tau^2",
                          "\n- Hedges estimator for tau^2",
                          "\n- Empirical Bayes estimator for tau^2")[i.lab.method.tau]
      ##
      if (tau.common)
        lab.method.tau <- paste(lab.method.tau, " (assuming common tau^2 in subgroups)", sep = "")
      ##
      if (hakn)
        lab.hakn <- "\n- Hartung-Knapp adjustment for random effects model"
      else
        lab.hakn <- ""
      ##      
      lab.method.details <- paste(lab.method.tau, lab.hakn, sep = "")
    }
  }
  ##
  imeth <- charmatch(method,
                     c("MH", "Peto", "Inverse", "Cochran", "GLMM"),
                     nomatch = NA)
  ##
  if ((metabin|metainc) & imeth == 1 & (sparse | addincr))
    if (MH.exact | metainc)
      lab.method.details <- paste(" (without continuity correction)",
                                  lab.method.details, sep = "")
  ##
  if (metacont && !is.null(pooledvar) && pooledvar)
    lab.method.details <- paste(" (with pooled variance for individual studies)",
                                lab.method.details, sep = "")
  ##
  method <- c("\n- Mantel-Haenszel method",
              "\n- Peto method",
              "\n- Inverse variance method",
              "\n- Cochran method",
              "GLMM")[imeth]
  ##
  if (method == "GLMM") {
    
    i.model.glmm <- charmatch(model.glmm,
                              c("UM.FS", "UM.RS", "CM.EL", "CM.AL"))
    ##
    if (metabin)
      method <- c("\n- Logistic regression model (fixed study effects)",
                  "\n- Mixed-effects logistic regression model (random study effects)",
                  "\n- Generalised linear mixed model (conditional Hypergeometric-Normal)",
                  "\n- Generalised linear mixed model (conditional Binomial-Normal)")[i.model.glmm]
    else if (metainc)
      method <- c("\n- Poisson regression model (fixed study effects)",
                  "\n- Mixed-effects Poisson regression model (random study effects)",
                  "\n- Generalised linear mixed model (conditional Poisson-Normal)")[i.model.glmm]
    else if (metaprop)
      method <- "\n- Random intercept logistic regression model"
    else if (metarate)
      method <- "\n- Random intercept Poisson regression model"
  }
  ##
  method <- paste(method, lab.method.details, sep = "")
  ##
  if (k.all > 1) {
    cat(paste("\nDetails on meta-analytical method:", method, sep = ""))
    ##
    if (trimfill)
      cat("\n- Trim-and-fill method to adjust for funnel plot asymmetry")
    ##
    if (metamiss) {
      if (method.miss == "IMOR")
        mmiss <- paste0("IMOR.e = ", IMOR.e, ", IMOR.c = ", IMOR.c)
      else {
        meths <- c("Gamble-Hollis analysis",
                   "impute no events (ICA-0)",
                   "impute events (ICA-1)",
                   "observed risk in control group (ICA-pc)",
                   "observed risk in experimental group (ICA-pe)",
                   "observed group-specific risks (ICA-p)",
                   "best-case scenario (ICA-b)",
                   "worst-case scenario (ICA-w)")
        ##
        mm <- c("GH", "0", "1", "pc", "pe", "p", "b", "w")
        ##
        idx <- charmatch(method.miss, mm)
        ##
        mmiss <- meths[idx]
      }
      ##
      cat(paste0("\n- Imputation method: ", mmiss))
    }
  }
  else
    cat(paste("\nDetails:", method, sep = ""))
  ##
  if (sm.details != "")
    cat(sm.details,
        "\n",
        sep = "")
  else
    cat("\n")
  
  invisible(NULL)
}
