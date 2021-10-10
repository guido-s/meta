.onAttach <- function(libname, pkgname) {
  msg <-
    paste0("Loading 'meta' package (version ",
           utils::packageDescription("meta")$Version,
           ").",
           "\nType 'help(meta)' for a brief overview.",
           "\nReaders of 'Use-R! Meta-Analysis with R' should install",
           "\nolder version of 'meta' package: ",
           "https://tinyurl.com/dt4y5drs")
  packageStartupMessage(msg)
}


updateversion <- function(x) {
  ##
  ## Update older meta objects
  ##
  if (!is.null(x$version))
    meta.version <- as.numeric(unlist(strsplit(x$version, "-"))[1])
  if (is.null(x$version) || meta.version < gs("version.update"))
    x <- update(x, warn = FALSE, warn.deprecated = FALSE)
  ##
  x
}


cathet <- function(k,
                   tau2, lower.tau2, upper.tau2,
                   print.tau2, print.tau2.ci, text.tau2, digits.tau2,
                   tau, lower.tau, upper.tau,
                   print.tau, print.tau.ci, text.tau, digits.tau,
                   sign.lower.tau, sign.upper.tau,
                   I2, lowI2, uppI2, 
                   print.I2, print.I2.ci, text.I2, digits.I2,
                   H, lowH, uppH,
                   print.H, digits.H,
                   Rb, lowRb, uppRb,
                   print.Rb, text.Rb,
                   big.mark,
                   detail.tau = "") {
  
  
  if (is.null(lower.tau2))
    lower.tau2 <- NA
  if (is.null(upper.tau2))
    upper.tau2 <- NA
  if (is.null(lower.tau))
    lower.tau <- NA
  if (is.null(upper.tau))
    upper.tau <- NA
  ##
  if (all(is.na(lower.tau2)) && all(is.na(upper.tau2)))
    print.tau2.ci <- FALSE
  if (all(is.na(lower.tau)) && all(is.na(upper.tau)))
    print.tau.ci <- FALSE
  
  
  stau <- length(tau) == 1
  ##
  if (!stau) {
    text.tau2 <- paste(text.tau2, seq_along(tau), sep = ".")
    text.tau <- paste(text.tau, seq_along(tau), sep = ".")
  }
  ##
  detail.tau <- ifelse(detail.tau != "", paste0(" (", detail.tau, ")"), "")
  
  
  cat(
    paste(
      if (print.tau2 | print.tau | print.I2 | print.H | print.Rb)
        " ",
      if (print.tau2)
        paste0(formatPT(tau^2,
                        lab = TRUE, labval = text.tau2,
                        digits = digits.tau2,
                        lab.NA = "NA",
                        big.mark = big.mark),
               if (print.tau2.ci)
                 pasteCI(lower.tau2, upper.tau2, digits.tau2, big.mark,
                         sign.lower.tau, sign.upper.tau),
               if (!print.tau) detail.tau),
      ##
      if (print.tau)
        paste0(
          if (print.tau2) "; " else "",
          formatPT(tau,
                   lab = TRUE, labval = text.tau,
                   digits = digits.tau,
                   lab.NA = "NA",
                   big.mark = big.mark),
          if (print.tau.ci)
            pasteCI(lower.tau, upper.tau, digits.tau, big.mark,
                    sign.lower.tau, sign.upper.tau),
          detail.tau),
      sep = "", collapse = "\n")
  )
  ##
  cat(
    paste0(
      if (print.I2)
        paste0(
          ifelse(
            print.tau2 | print.tau,
          ifelse(!stau | print.tau2.ci | print.tau.ci |
                 (options()$width < 70 & print.I2.ci),
                 "\n", ";"),
          ""),
          if (print.tau2 | print.tau)
            " ",
          text.I2, " = ",
          if (is.na(I2))
            "NA"
          else
            paste0(formatN(I2, digits.I2), "%"),
          if (print.I2.ci)
            pasteCI(lowI2, uppI2, digits.I2, big.mark, unit = "%")
        ),
      ##
      if (print.H)
        paste0(
          if (print.tau2 | print.tau | print.I2)
            "; ",
          "H = ",
          if (is.na(H))
            "NA"
          else
            formatN(H, digits.H, "NA", big.mark = big.mark),
          if (!(is.na(lowH) | is.na(uppH)))
            pasteCI(lowH, uppH, digits.H, big.mark)
        ),
      ##
      if (print.Rb)
        paste0(
          if (print.tau2 | print.tau | print.I2 | print.H)
            ";\n",
          text.Rb, " = ",
          if (is.na(Rb))
            "NA"
          else
            paste0(formatN(Rb, digits.I2, big.mark = big.mark), "%"),
          if (!(is.na(lowRb) | is.na(uppRb)))
            pasteCI(lowRb, uppRb, digits.I2, big.mark, unit = "%")
        ),
      ##
      if (print.tau2 | print.tau | print.I2 | print.H | print.Rb)
        "\n"
    )
  )
  
  
  invisible(NULL)
}


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
                    exact.smd = FALSE,
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
        charmatch(method.tau, c(.settings$meth4tau, ""), nomatch = NA)
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
  else {
    if (method != "" | sm.details != "")
      cat("\nDetails:")
    if (method != "")
      cat(method)
  }
  ##
  if (sm.details != "")
    cat(paste0(sm.details, "\n"))
  else if (method != "")
    cat("\n")
  
  
  invisible(NULL)
}


## The following R code is based on the file snowfall-internal.R from
## R package snowfall (Maintainer: Jochen Knaus <jo@imbi.uni-freiburg.de>)
##
## Helpers for managing the internal variables in the package namespace without
## awake the R CMD check for later R versions (which basically blaims many
## global assignments).
##
## The given solution has an advantage: only writing is affected. Reading of the
## objects can remain the same (thanks to Uwe Ligges for the tipp):
##   reading:  gs("CIbracket")
##   writing:  setOption("CIbracket", "(")
##

##
## Set an option in the meta option list.
## (Basically this is the setting of a list entry).
## key - character: object name
## val - object (everything is allowed, even NULL)
##
setOption <- function(key = NULL, val = NULL) {
  if(!is.null(key) && is.character(key)) {
    option <- getVar(".settings") # Get from NS
    option[[key]] <- val
    setVar(".settings", option) # Write to NS
    ##
    return(invisible(TRUE))
  }
  ##
  stop("Argument 'key' or 'val' is NULL or 'key' is no string.")
}


##
## Get a specific variable from the meta namespace.
## var - character: object name
##
getVar <- function(var = NULL) {
  if(!is.null(var) && is.character(var)) {
    tmp <- try(getFromNamespace(var, "meta"))
    ##
    if(inherits(tmp, "try-error"))
      stop("Object", var, "not found in meta package.")
    ##
    return(tmp)
  }

  stop("Argument 'var' is NULL or not a string.")
}


##
## Write a specific variable to the meta namespace.
## var - character: object name
## arg - object (NULL allowed)
##
setVar <- function(var = NULL, arg = NULL) {
  if(!is.null(var) && is.character(var)) {
    assignInNamespace(var, arg, "meta")

    return(invisible(TRUE))
  }

  stop("var is NULL or no character");
}






.settings <- list()
##
## List of internal settings
##
argslist.internal <-
  c("comb.fixed", "comb.random", "level.comb", "digits.zval",
    "print.byvar", "byseparator",
    "Wan2014.Table1", "Wan2014.Table2",
    "sm4bin", "sm4cont", "sm4cor", "sm4inc", "sm4mean", "sm4prop", "sm4rate",
    "ci4cont", "ci4prop",
    "meth4bin", "meth4inc", "meth4prop", "meth4rate",
    "meth4tau", "meth4tau.ci",
    "adhoc4hakn",
    "meth4bias", "meth4bias.old",
    "version.update")
##
setOption("argslist.internal", argslist.internal)
##
## Set defaults (for internal options)
##
setOption("sm4bin", c("OR", "RD", "RR", "ASD", "DOR"))
setOption("sm4cont", c("MD", "SMD", "ROM"))
setOption("sm4cor", c("ZCOR", "COR"))
setOption("sm4inc", c("IRR", "IRD", "IRSD"))
setOption("sm4mean", c("MRAW", "MLN"))
setOption("sm4prop", c("PLOGIT", "PLN", "PRAW", "PAS", "PFT"))
setOption("sm4rate", c("IR", "IRLN", "IRS", "IRFT"))
##
setOption("ci4cont", c("z", "t"))
setOption("ci4prop", c("CP", "WS", "WSCC", "AC", "SA", "SACC", "NAsm"))
##
setOption("meth4bin", c("Inverse", "MH", "Peto", "GLMM", "SSW"))
setOption("meth4inc", c("Inverse", "MH", "Cochran", "GLMM"))
setOption("meth4prop", c("Inverse", "GLMM"))
setOption("meth4rate", c("Inverse", "GLMM"))
##
setOption("meth4tau", c("DL", "PM", "REML", "ML", "HS", "SJ", "HE", "EB"))
setOption("meth4tau.ci", c("QP", "BJ", "J", "PL", ""))
setOption("adhoc4hakn", c("", "se", "ci", "iqwig6"))
##
setOption("meth4bias.old", c("rank", "linreg", "mm", "count", "score"))
setOption("meth4bias", c("Begg", "Egger", "Thompson", "Schwarzer",
                         "Harbord", "Peters", "Deeks",
                         "Pustejovsky", "Macaskill"))
##
setOption("version.update", 5.0)
##
## List of arguments that can be changed by user
##
argslist <-
  c("level", "level.ma", "fixed", "random",
    "hakn", "adhoc.hakn", "method.tau", "method.tau.ci", "tau.common",
    "prediction", "level.predict",
    "method.bias",
    "text.fixed", "text.random", "text.predict",
    "text.w.fixed", "text.w.random",
    "title", "complab", "CIbracket", "CIseparator",
    "print.subgroup.name", "sep.subgroup",
    "keepdata", "warn", "warn.deprecated",
    "backtransf",
    "smbin", "smcont", "smcor", "sminc", "smmean", "smprop", "smrate",
    "incr", "allincr", "addincr",
    "method", "allstudies", "MH.exact",
    "RR.Cochrane", "Q.Cochrane", "model.glmm", "print.CMH",
    "pooledvar", "method.smd", "sd.glass", "exact.smd",
    "method.ci.cont", "method.ci.prop",
    "label.e", "label.c", "label.left", "label.right",
    "layout",
    "test.overall", "test.subgroup", "test.effect.subgroup",
    "digits", "digits.se", "digits.zval", "digits.stat",
    "digits.Q", "digits.tau2", "digits.tau", "digits.H", "digits.I2",
    "digits.prop", "digits.weight",
    "digits.pval", "digits.pval.Q", "digits.forest",
    "scientific.pval", "big.mark", "zero.pval", "JAMA.pval",
    "print.I2", "print.H", "print.Rb",
    "text.tau2", "text.tau", "text.I2", "text.Rb"
    )
##
setOption("argslist", argslist)
##
## General settings
##
setOption("level", 0.95)
setOption("level.ma", 0.95)
setOption("level.comb", 0.95)
setOption("fixed", TRUE)
setOption("comb.fixed", TRUE)
setOption("random", TRUE)
setOption("comb.random", TRUE)
setOption("hakn", FALSE)
setOption("adhoc.hakn", "")
setOption("method.tau", "REML")
setOption("method.tau.ci", NULL)
setOption("tau.common", FALSE)
setOption("prediction", FALSE)
setOption("level.predict", 0.95)
setOption("method.bias", "Egger")
setOption("text.fixed", "Common effect model")
setOption("text.random", "Random effects model")
setOption("text.predict", "Prediction interval")
setOption("text.w.fixed", "common")
setOption("text.w.random", "random")
setOption("title", "")
setOption("complab", "")
setOption("CIbracket", "[")
setOption("CIseparator", "; ")
setOption("print.subgroup.name", TRUE)
setOption("print.byvar", TRUE)
setOption("sep.subgroup", " = ")
setOption("byseparator", " = ")
setOption("keepdata", TRUE)
setOption("warn", TRUE)
setOption("warn.deprecated", TRUE)
setOption("backtransf", TRUE)
setOption("digits", 4)
setOption("digits.se", 4)
setOption("digits.zval", 2)
setOption("digits.stat", 2)
setOption("digits.Q", 2)
setOption("digits.tau2", 4)
setOption("digits.tau", 4)
setOption("digits.H", 2)
setOption("digits.I2", 1)
setOption("digits.prop", 4)
setOption("digits.weight", 1)
setOption("digits.pval", 4)
setOption("digits.pval.Q", 4)
setOption("scientific.pval", FALSE)
setOption("big.mark", "")
setOption("zero.pval", TRUE)
setOption("JAMA.pval", FALSE)
setOption("print.I2", TRUE)
setOption("print.H", TRUE)
setOption("print.Rb", FALSE)
setOption("text.tau2", "tau^2")
setOption("text.tau", "tau")
setOption("text.I2", "I^2")
setOption("text.Rb", "Rb")
##
## Default summary measure
##
setOption("smbin", "RR")
setOption("smcont", "MD")
setOption("smcor", "ZCOR")
setOption("sminc", "IRR")
setOption("smmean", "MRAW")
setOption("smprop", "PLOGIT")
setOption("smrate", "IRLN")
##
## Settings for R functions metabin, metainc, metaprop
##
setOption("incr", 0.5)
setOption("allincr", FALSE)
setOption("addincr", FALSE)
##
## Additional settings for R function metabin
##
setOption("method", "MH")
setOption("allstudies", FALSE)
setOption("MH.exact", FALSE)
setOption("RR.Cochrane", FALSE)
setOption("Q.Cochrane", FALSE)
setOption("model.glmm", "UM.FS")
setOption("print.CMH", FALSE)
##
## Additional setting for R function metacont
##
setOption("pooledvar", FALSE)
setOption("method.smd", "Hedges")
setOption("sd.glass", "control")
setOption("exact.smd", FALSE)
setOption("method.ci.cont", "z")
##
## Additional setting for R function metaprop
##
setOption("method.ci.prop", "CP")
##
## Settings for R functions comparing two treatments
##
setOption("label.e", "Experimental")
setOption("label.c", "Control")
setOption("label.left", "")
setOption("label.right", "")
##
## Settings for R function forest.meta
##
setOption("layout", "meta")
setOption("test.overall", FALSE)
setOption("test.subgroup", TRUE)
setOption("test.effect.subgroup", FALSE)
setOption("digits.forest", 2)


setOption("Wan2014.Table1",
          c(0.000, 1.128, 1.693, 2.059, 2.326,
            2.534, 2.704, 2.847, 2.970, 3.078,
            3.173, 3.259, 3.336, 3.407, 3.472,
            3.532, 3.588, 3.640, 3.689, 3.735,
            3.778, 3.819, 3.858, 3.895, 3.931,
            3.964, 3.997, 4.027, 4.057, 4.086,
            4.113, 4.139, 4.165, 4.189, 4.213,
            4.236, 4.259, 4.280, 4.301, 4.322,
            4.341, 4.361, 4.379, 4.398, 4.415,
            4.433, 4.450, 4.466, 4.482, 4.498))
##
setOption("Wan2014.Table2",
          c(0.990, 1.144, 1.206, 1.239, 1.260,
            1.274, 1.284, 1.292, 1.298, 1.303,
            1.307, 1.311, 1.313, 1.316, 1.318,
            1.320, 1.322, 1.323, 1.324, 1.326,
            1.327, 1.328, 1.329, 1.330, 1.330,
            1.331, 1.332, 1.332, 1.333, 1.333,
            1.334, 1.334, 1.335, 1.335, 1.336,
            1.336, 1.336, 1.337, 1.337, 1.337,
            1.338, 1.338, 1.338, 1.338, 1.339,
            1.339, 1.339, 1.339, 1.339, 1.340))
