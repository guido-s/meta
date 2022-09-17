catmeth <- function(method,
                    method.tau = NULL,
                    sm = "",
                    k.all,
                    method.random.ci = "",
                    df.random = NA,
                    adhoc.hakn.ci = "",
                    method.predict = "",
                    adhoc.hakn.pi = "",
                    df.predict = NA,
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
                    RR.Cochrane = FALSE,
                    Q.Cochrane = FALSE,
                    method.ci = NULL,
                    print.tau.ci = FALSE,
                    method.tau.ci = "",
                    pooledvar = FALSE,
                    method.smd,
                    sd.glass,
                    exact.smd = TRUE,
                    model.glmm,
                    pscale = 1,
                    irscale = 1,
                    irunit = "person-years",
                    null.effect = NA,
                    big.mark = "",
                    digits = gs("digits"),
                    digits.tau = gs("digits.tau"),
                    text.tau = gs("text.tau"),
                    text.tau2 = gs("text.tau2"),
                    method.miss, IMOR.e, IMOR.c,
                    three.level = FALSE
                    ) {
  
  
  metabin  <- "metabin"  %in% class
  metacont <- "metacont" %in% class
  metainc  <- "metainc"  %in% class
  metaprop <- "metaprop" %in% class
  metarate <- "metarate" %in% class
  trimfill <- "trimfill" %in% class
  metamiss <- "metamiss" %in% class
  ##
  allstudies <- replaceNULL(allstudies, FALSE)
  doublezeros <- replaceNULL(doublezeros, FALSE)
  ##
  method.ci <- replaceNULL(method.ci, "")
  ##
  method.random.ci <- replaceNULL(method.random.ci, "")
  adhoc.hakn.ci <- replaceNULL(adhoc.hakn.ci, "")
  method.predict <- replaceNULL(method.predict, "")
  adhoc.hakn.pi <- replaceNULL(adhoc.hakn.pi, "")
  
  
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
  if (method.ci == "CP")
    method.ci.details <-
      "\n- Clopper-Pearson confidence interval for individual studies"
  else if (method.ci == "WS")
    method.ci.details <-
      "\n- Wilson Score confidence interval for individual studies"
  else if (method.ci == "WSCC")
    method.ci.details <-
      paste0("\n- Wilson Score confidence interval with ",
             "continuity correction\n", "  for individual studies")
  else if (method.ci == "AC")
    method.ci.details <-
      "\n- Agresti-Coull confidence interval for individual studies"
  else if (method.ci == "SA")
    method.ci.details <-
      "\n- Simple approximation confidence interval for individual studies"
  else if (method.ci == "SACC")
    method.ci.details <-
      paste0("\n- Simple approximation confidence interval with ",
             "continuity correction for individual studies")
  else if (method.ci == "NAsm")
    method.ci.details <-
      "\n- Normal approximation confidence interval for individual studies"
  else if (method.ci == "Poisson")
    method.ci.details <-
      "\n- Exact Poisson confidence interval for individual studies"
  else if (method.ci == "t")
    method.ci.details <-
      "\n- Confidence interval for individual studies based on t-distribution"
  else
    method.ci.details <- ""
  ##
  sm.details <- paste0(sm.details, method.ci.details)
  ##
  if (metacont && sm == "SMD" && !is.null(method.smd)) {
    if (method.smd == "Hedges")
      if (exact.smd)
        method.details <-
          paste0("\n- Hedges' g (bias corrected standardised mean difference; ",
                 "using exact formulae)")
      else
        method.details <-
          "\n- Hedges' g (bias corrected standardised mean difference)"
    else if (method.smd == "Cohen")
      if (exact.smd)
        method.details <-
          "\n- Cohen's d (standardised mean difference; using exact formulae)"
      else
        method.details <-
          "\n- Cohen's d (standardised mean difference)"
    else if  (method.smd == "Glass") {
      if (!is.null(sd.glass) && sd.glass == "control")
        method.details <-
          paste0("\n- Glass' delta (standardised mean difference; ",
                 "based on control group)")
      else if (!is.null(sd.glass) && sd.glass == "experimental")
        method.details <-
          paste0("\n- Glass' delta (standardised mean difference; ",
                 "based on experimental group)")
    }
    ##
    sm.details <- paste0(sm.details, method.details)
  }
  ##
  if (metabin | metainc | metaprop | metarate) {
    txtCC <- !(any(method == "MH") & MH.exact & all(k.all == 1))
    txtCC.ind <- (any(method == "MH") & MH.exact) | any(method == "GLMM")
    ##
    if (!(sm == "ASD" | any(method == "Peto"))) {
      if (addincr) {
        if (all(incr == "TACC") && txtCC)
          sm.details <-
            paste0(sm.details,
                   "\n- Treatment arm continuity correction in all studies",
                   if (txtCC.ind)
                     "\n  (only used to calculate individual study results)")
        else if (all(incr != 0) && txtCC)
          sm.details <-
            paste0(sm.details,
                   "\n- Continuity correction",
                   if (length(unique(incr)) == 1)
                     paste0(" of ", round(incr, 4)),
                   " in all studies",
                   if (txtCC.ind)
                     "\n  (only used to calculate individual study results)")
      }
      else if (sparse) {
        if (allincr == FALSE) {
          if (all(incr == "TACC") && txtCC)
            sm.details <-
              paste0(sm.details,
                     paste0("\n- Treatment arm continuity correction in ",
                            "studies with",
                            if (options()$width > 70) " " else "\n  ",
                            "zero cell frequencies"),
                     if (txtCC.ind)
                       "\n  (only used to calculate individual study results)")
          else if (any(incr != 0) && txtCC)
            sm.details <-
              paste0(sm.details,
                     "\n- Continuity correction",
                     if (length(unique(incr)) == 1)
                       paste0(" of ", round(incr, 4)),
                     " in studies with",
                     if (options()$width > 70) " " else "\n  ",
                     "zero cell frequencies",
                     if (txtCC.ind)
                       "\n  (only used to calculate individual study results)")
        }
        else {
          if (all(incr == "TACC")) {
            if (txtCC)
              sm.details <-
                paste0(sm.details,
                       "\n- Treatment arm continuity correction in all studies",
                       if (txtCC.ind)
                         "\n  (only used to calculate individual study results)")
          }
          else if (any(incr != 0) && txtCC)
            sm.details <-
              paste0(sm.details,
                     "\n- Continuity correction",
                     if (length(unique(incr)) == 1)
                       paste0(" of ", round(incr, 4)),
                     " in all studies",
                     if (txtCC.ind)
                       "\n  (only used to calculate individual study results)")
        }
        ##
        if (allstudies & doublezeros)
          sm.details <-
            paste0(sm.details,
                  "\n- Studies with double zeros included in meta-analysis")
      }
    }
  }
  
  
  if (pscale != 1)
    sm.details <- paste0(sm.details,
                         "\n- Events per ",
                         format(pscale, scientific = FALSE,
                                big.mark = big.mark),
                         " observations")
  ##
  if (irscale != 1)
    sm.details <- paste0(sm.details,
                        "\n- Events per ",
                        format(irscale, scientific = FALSE,
                               big.mark = big.mark),
                        " ", irunit)

  if (!is.na(null.effect) && null.effect != 0) {
    if (pscale != 1)
      sm.details <- paste0(sm.details,
                           "\n- Null hypothesis: effect is equal to ",
                           format(round(null.effect * pscale, digits),
                                  scientific = FALSE, big.mark = big.mark),
                           " events per ",
                           format(pscale, scientific = FALSE,
                                  big.mark = big.mark),
                           " observations")
    else if (irscale != 1)
      sm.details <- paste0(sm.details,
                           "\n- Null hypothesis: effect is equal to ",
                           format(round(null.effect * irscale, digits),
                                  scientific = FALSE, big.mark = big.mark),
                           " events per ",
                           format(irscale, scientific = FALSE,
                                  big.mark = big.mark),
                           " ", irunit)
    else
      sm.details <- paste0(sm.details,
                           "\n- Null hypothesis: effect is equal to ",
                           format(null.effect, scientific = FALSE,
                                  big.mark = big.mark))
  }
  
  
  ##
  ## Methods to estimate tau2
  ##
  lab.method.tau <- ""
  ##
  if (!is.null(tau.preset)) {
    tau.preset <- formatPT(tau.preset, lab = TRUE, labval = text.tau,
                           digits = digits.tau,
                           lab.NA = "NA",
                           big.mark = big.mark)
    ##
    lab.method.tau <-
      paste0("\n- Preset square root of between-study variance: ", tau.preset)
  }
  else {
    if (any(method.tau == "DL"))
      lab.method.tau <- "\n- DerSimonian-Laird estimator"
    ##
    if (any(method.tau == "PM"))
      lab.method.tau <-
        paste0(lab.method.tau, "\n- Paule-Mandel estimator")
    ##
    if (any(method.tau == "REML"))
      lab.method.tau <-
        paste0(lab.method.tau, "\n- Restricted maximum-likelihood estimator")
    ##
    if (any(method.tau == "ML"))
      lab.method.tau <-
        paste0(lab.method.tau, "\n- Maximum-likelihood estimator")
    ##
    if (any(method.tau == "HS"))
      lab.method.tau <-
        paste0(lab.method.tau, "\n- Hunter-Schmidt estimator")
    ##
    if (any(method.tau == "SJ"))
      lab.method.tau <-
        paste0(lab.method.tau, "\n- Sidik-Jonkman estimator")
    ##
    if (any(method.tau == "HE"))
      lab.method.tau <-
        paste0(lab.method.tau, "\n- Hedges estimator")
    ##
    if (any(method.tau == "EB"))
      lab.method.tau <-
        paste0(lab.method.tau, "\n- Empirical Bayes estimator")
    ##
    if (any(lab.method.tau != "")) {
      lab.method.tau <- paste(lab.method.tau, "for", text.tau2)
      if (tau.common)
        lab.method.tau <-
          paste0(lab.method.tau,
                 "\n  (assuming common ", text.tau2, " in subgroups)")
    }
    ##
    if (any(method.tau == "DL") & metabin & Q.Cochrane)
      lab.method.tau <-
        paste0(lab.method.tau,
               "\n- Mantel-Haenszel estimator used in ",
               "calculation of Q and ", text.tau2,
               if (options()$width > 70)
                 " (like RevMan 5)"
               else
                 "\n  (like RevMan 5)")
  }

  
  ##
  ## Methods to calculate confidence interval for tau2
  ##
  lab.tau.ci <- ""
  ##
  if (print.tau.ci) {
    if (any(method.tau.ci == "QP"))
      lab.tau.ci <- "\n- Q-Profile method"
    ##
    if (any(method.tau.ci == "BJ"))
      lab.tau.ci <-
        paste0(lab.tau.ci, "\n- Biggerstaff and Jackson method")
    ##
    if (any(method.tau.ci == "J"))
      lab.tau.ci <-
        paste0(lab.tau.ci, "\n- Jackson method")
    ##
    if (any(method.tau.ci == "PL"))
      lab.tau.ci <-
        paste0(lab.tau.ci, "\n- Profile-Likelihood method")
    ##
    lab.tau.ci <-
      paste(lab.tau.ci,
            "for confidence interval of", text.tau2, "and", text.tau)
  }
  
  
  ##
  ## Method to calculate random effects confidence interval
  ##
  lab.random.ci <- ""
  ##
  if (any(method.random.ci == "HK")) {
    lab.random.ci <-
      paste0("\n- Hartung-Knapp (HK) adjustment for ",
             "random effects model (df = ",
             rmSpace(unique(df.random[method.random.ci == "HK"])),
             ")")
    ##
    if (any(adhoc.hakn.ci != ""))
      lab.random.ci <-
        paste0(lab.random.ci,
               "\n  (with ",
               if (any(method.random.ci == "HK" & adhoc.hakn.ci == ""))
                 "and without ",
               "ad hoc correction)")
  }
  ##
  if (any(method.random.ci == "classic-KR")) {
    lab.random.ci <-
      paste0(lab.random.ci,
             "\n- Classic method instead of Kenward-Roger adjustment ",
             "(classic-KR) used for random effects model")
  }
  ##
  if (any(method.random.ci == "KR")) {
    lab.random.ci <-
      paste0(lab.random.ci,
             "\n- Kenward-Roger (KR) adjustment for ",
             "random effects model (df = ",
             rmSpace(unique(df.random[method.random.ci == "KR"])),
             ")")
  }
  
  
  ##
  ## Method to calculate prediction interval
  ##
  lab.predict <- ""
  ##
  if (any(method.predict == "HTS")) {
    lab.predict <-
      paste0("\n- Prediction interval based on t-distribution (HTS) ",
             "(df = ",
             rmSpace(unique(df.predict[method.predict == "HTS"])),
             ")")
  }
  ##
  if (any(method.predict == "HK")) {
    lab.predict <-
      paste0(lab.predict,
             paste0("\n- Hartung-Knapp (HK) prediction interval (df = ",
                    rmSpace(unique(df.predict[method.predict == "HK"])),
                    ")"))
    if (any(adhoc.hakn.pi != ""))
      lab.predict <-
        paste0(lab.predict,
               "\n  (with ",
               if (any(method.predict == "HK" &
                       adhoc.hakn.pi == ""))
                 "and without ",
               "ad hoc correction)")
  }
  ##
  if (any(method.predict == "HTS-KR")) {
    lab.predict <-
      paste0("\n- Prediction interval based on t-distribution (df = ",
             rmSpace(unique(df.predict[method.predict == "HTS"])),
             ") instead of ",
             "Kenward-Roger adjustment")
  }
  ##
  if (any(method.predict == "KR")) {
    lab.predict <-
      paste0(lab.predict,
             "\n- Kenward-Roger (KR) prediction interval (df = ",
             rmSpace(unique(df.predict[method.predict == "KR"])),
             ")")
  }
  ##
  if (any(method.predict == "HTS-KR")) {
    lab.predict <-
      paste0(lab.predict,
             "\n- Kenward-Roger (KR) prediction interval (df = ",
             rmSpace(unique(df.predict[method.predict == "KR"])),
             ")")
  }
  ##
  if (any(method.predict == "NNF")) {
    lab.predict <-
      paste0(lab.predict,
             "\n- Boot-strap prediction interval (NNF) (df = ",
             rmSpace(unique(df.predict[method.predict == "NNF"])),
             ")")
  }
  ##
  if (any(method.predict == "S")) {
    lab.predict <-
      paste0(lab.predict,
             "\n- Prediction interval based on ",
                 "standard normal distribution (S)")
  }
  
  
  ##
  ## Meta-analysis method
  ##
  lab.method <- ""
  ##
  if (any(method == "MH")) {
    lab.method <- "\n- Mantel-Haenszel method"
    if ((metabin | metainc) & (sparse | addincr) & MH.exact)
      lab.method <-
        paste(lab.method, "(without continuity correction)")
  }
  ##
  if (any(method == "Peto"))
    lab.method <- paste0(lab.method, "\n- Peto method")
  ##
  if (any(method == "Inverse")) {
    lab.method <-
      paste0(lab.method, "\n- Inverse variance method")
    if (three.level)
      lab.method <-
        paste(lab.method, "(three-level model)")
    if (metacont && !is.null(pooledvar) && pooledvar)
      lab.method <-
        paste(lab.method,
              "(with pooled variance for individual studies)")
  }
  ##
  if (any(method == "Cochran"))
    lab.method <- paste0(lab.method, "\n- Cochran method")
  ##
  if (any(method == "SSW"))
    lab.method <- paste0(lab.method, "\n- Sample size method")
  ##
  if (any(method == "GLMM")) {
    if (metabin)
      lab.method <-
        paste0(lab.method,
               if (any(model.glmm == "UM.FS"))
                 "\n- Logistic regression model (fixed study effects)",
               if (any(model.glmm == "UM.RS"))
                 paste("\n- Mixed-effects logistic regression model",
                       "(random study effects)"),
               if (any(model.glmm == "CM.EL"))
                 paste("\n- Generalised linear mixed model",
                       "(conditional Hypergeometric-Normal)"),
               if (any(model.glmm == "CM.AL"))
                 paste("\n- Generalised linear mixed model",
                       "(conditional Binomial-Normal)"))
    else if (metainc)
      lab.method <-
        paste0(lab.method,
               if (any(model.glmm == "UM.FS"))
                 "\n- Poisson regression model (fixed study effects)",
               if (any(model.glmm == "UM.RS"))
                 paste("\n- Mixed-effects Poisson regression model",
                       "(random study effects)"),
               if (any(model.glmm == "CM.EL"))
                 paste("\n- Generalised linear mixed model",
                       "(conditional Poisson-Normal)"))
    else if (metaprop)
      lab.method <-
        paste0(lab.method,
               "\n- Random intercept logistic regression model")
    else if (metarate)
      lab.method <-
        paste0(lab.method,
               "\n- Random intercept Poisson regression model")
  }
  
  
  details <-
    paste0(lab.method, lab.method.tau, lab.tau.ci, lab.random.ci, lab.predict)
  ##
  if (any(k.all > 1)) {
    if (details != "")
      cat(paste0("\nDetails on meta-analytical method:", details))
    ##
    if (trimfill)
      cat("\n- Trim-and-fill method to adjust for funnel plot asymmetry")
  }
  else {
    if (details != "" | sm.details != "")
      cat("\nDetails:")
    if (details != "")
      cat(details)
  }
  ##
  if (metamiss) {
    if (method.miss == "IMOR") {
      mmiss <- "Informative Missingness Odds Ratio"
      if (length(unique(IMOR.e)) == 1 & length(unique(IMOR.c)) == 1)
        mmiss <-
          paste0("\n  ", mmiss,
                 " (IMOR.e = ", round(unique(IMOR.e), 4),
                 ", IMOR.c = ", round(unique(IMOR.c), 4), ")")
      else
        mmiss <- paste(mmiss, "(IMOR)")
    }
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
  ##
  if (sm.details != "")
    cat(paste0(sm.details, "\n"))
  else if (details != "")
    cat("\n")
  
  
  invisible(NULL)
}


NULL
