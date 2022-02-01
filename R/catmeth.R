catmeth <- function(method,
                    method.tau = NULL,
                    sm = "",
                    k.all,
                    hakn = FALSE,
                    adhoc.hakn = FALSE,
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
  if (!is.null(method.ci)) {
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
    else if (method.ci == "t")
      method.ci.details <-
        "\n- Confidence interval for individual studies based on t-distribution"
    else
      method.ci.details <- ""
    ##
    sm.details <- paste0(sm.details, method.ci.details)
  }
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
    txtCC <- !(method == "MH" & MH.exact & k.all == 1)
    txtCC.ind <- (method == "MH" & MH.exact) | method == "GLMM"
    ##
    if (!(sm == "ASD" | method == "Peto")) {
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
            paste(sm.details,
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
  
  
  lab.method.details <- ""
  ##
  if (is.null(method.tau))
    lab.method.tau <- ""
  else {
    if (!is.null(tau.preset)) {
      tau.preset <- formatPT(tau.preset, lab = TRUE, labval = text.tau,
                             digits = digits.tau,
                             lab.NA = "NA",
                             big.mark = big.mark)
      ##
      lab.method.tau <-
        paste0("\n- Preset square root of between-study variance: ", tau.preset)
      ##
      lab.method.details <- lab.method.tau
    }
    else {
      i.lab.method.tau <-
        charmatch(method.tau, c(gs("meth4tau"), ""), nomatch = NA)
      ##
      lab.method.tau <-
        c("\n- DerSimonian-Laird estimator",
          "\n- Paule-Mandel estimator",
          "\n- Restricted maximum-likelihood estimator",
          "\n- Maximum-likelihood estimator",
          "\n- Hunter-Schmidt estimator",
          "\n- Sidik-Jonkman estimator",
          "\n- Hedges estimator",
          "\n- Empirical Bayes estimator",
          "")[i.lab.method.tau]
      if (lab.method.tau != "")
        lab.method.tau <- paste(lab.method.tau, "for", text.tau2)
      ##
      if (lab.method.tau != "" & tau.common)
        lab.method.tau <-
          paste0(lab.method.tau,
                 if (options()$width <= 70 || i.lab.method.tau == 3)
                   "\n  " else " ",
                 "(assuming common ", text.tau2,
                 " in subgroups)")
      ##
      if (lab.method.tau != "" & metabin & Q.Cochrane)
        lab.method.tau <-
          paste0(lab.method.tau,
                 "\n- Mantel-Haenszel estimator used in ",
                 "calculation of Q and ", text.tau2,
                 if (options()$width > 70)
                   " (like RevMan 5)"
                 else
                   "\n  (like RevMan 5)")
      ##
      if (print.tau.ci & method.tau.ci %in% c("QP", "BJ", "J")) {
        i.lab.tau.ci <-
          charmatch(method.tau.ci, c("QP", "BJ", "J"), nomatch = NA)
        ##
        lab.tau.ci <-
          paste("\n-",
                 c("Q-profile method",
                   "Biggerstaff and Jackson method",
                   "Jackson method")[i.lab.tau.ci],
                 "for confidence interval of",
                 text.tau2, "and", text.tau)
        ##
        lab.method.tau <- paste0(lab.method.tau, lab.tau.ci)
      }
      ##
      if (hakn) {
        lab.hakn <- "\n- Hartung-Knapp adjustment for random effects model"
        if (adhoc.hakn)
          lab.hakn <- paste0(lab.hakn,
                             "\n  (with ad hoc variance correction)")
      }
      else
        lab.hakn <- ""
      ##      
      lab.method.details <- paste0(lab.method.tau, lab.hakn)
    }
  }
  ##
  imeth <-
    charmatch(method,
              c("MH", "Peto", "Inverse", "Cochran", "SSW", "GLMM",
                "NoMA", ""),
              nomatch = NA)
  ##
  if ((metabin|metainc) & imeth == 1 & (sparse | addincr))
    if (MH.exact | metainc)
      lab.method.details <-
        paste0(" (without continuity correction)", lab.method.details)
  if (three.level)
    lab.method.details <-
      paste0(" (three-level model)", lab.method.details)
  ##
  if (metacont && !is.null(pooledvar) && pooledvar)
    lab.method.details <-
      paste0(" (with pooled variance for individual studies)",
             lab.method.details)
  ##
  method <- c("\n- Mantel-Haenszel method",
              "\n- Peto method",
              "\n- Inverse variance method",
              "\n- Cochran method",
              "\n- Sample size method",
              "GLMM",
              "",
              "")[imeth]
  ##
  if (method == "GLMM") {
    
    i.model.glmm <- charmatch(model.glmm,
                              c("UM.FS", "UM.RS", "CM.EL", "CM.AL"))
    ##
    if (metabin)
      method <-
        c("\n- Logistic regression model (fixed study effects)",
          "\n- Mixed-effects logistic regression model (random study effects)",
          paste("\n- Generalised linear mixed model",
                "(conditional Hypergeometric-Normal)"),
          paste("\n- Generalised linear mixed model",
                "(conditional Binomial-Normal)")
          )[i.model.glmm]
    else if (metainc)
      method <-
        c("\n- Poisson regression model (fixed study effects)",
          "\n- Mixed-effects Poisson regression model (random study effects)",
          "\n- Generalised linear mixed model (conditional Poisson-Normal)"
          )[i.model.glmm]
    else if (metaprop)
      method <- "\n- Random intercept logistic regression model"
    else if (metarate)
      method <- "\n- Random intercept Poisson regression model"
  }
  ##
  method <- paste0(method, lab.method.details)
  ##
  if (k.all > 1) {
    if (method != "")
      cat(paste0("\nDetails on meta-analytical method:", method))
    ##
    if (trimfill)
      cat("\n- Trim-and-fill method to adjust for funnel plot asymmetry")
  }
  else {
    if (method != "" | sm.details != "")
      cat("\nDetails:")
    if (method != "")
      cat(method)
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
  else if (method != "")
    cat("\n")
  
  
  invisible(NULL)
}


NULL
