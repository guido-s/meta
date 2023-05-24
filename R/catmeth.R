catmeth <- function(x,
                    common, random, prediction, overall, overall.hetstat,
                    func.transf, backtransf, func.backtransf,
                    big.mark, digits, digits.tau, text.tau, text.tau2,
                    print.tau.ci = FALSE
                    ) {

  ##
  ##
  ## (1) Some settings
  ##
  ##
  
  metabin  <- inherits(x, "metabin")
  metacont <- inherits(x, "metacont")
  metainc  <- inherits(x, "metainc")
  metaprop <- inherits(x, "metaprop")
  metarate <- inherits(x, "metarate")
  trimfill <- inherits(x, "trimfill")
  metamiss <- inherits(x, "metamiss")
  metabind <- inherits(x, "metabind")
  ##
  by <- !is.null(x$subgroup)
  ##
  gm <- gm(x)
  meth <- gm$meth
  pred <- gm$pred
  ##
  selmod <-
    if (!common & !random)
      "common"
    else if (common & !random)
      "common"
    else if (!common & random)
      "random"
    else
      c("common", "random")
  ##
  details <- ""
  
  
  ##
  ##
  ## (2) Meta-analysis methods
  ##
  ##

  if (common & (overall | metabind | by)) {
    
    ##
    ## Mantel-Haenszel method
    ##
    if ((metabin | metainc)) {
      if (metabin)
        vars <- c("incr", "method.incr", "sparse", "MH.exact")
      else
        vars <- c("incr", "method.incr", "sparse")
      ##
      dat.mh <-
        unique(meth[meth$model == "common" & meth$method == "MH", vars])
      ##
      for (i in seq_len(nrow(dat.mh))) {
        details <- paste0(details, "\n- Mantel-Haenszel method")
        ##
        if ((dat.mh$sparse[i] | dat.mh$method.incr[i] == "all") &
            (!is.null(dat.mh$MH.exact) && dat.mh$MH.exact[i]))
          details <-
            paste0(details,
                   if (random)
                     ", without continuity correction)"
                   else
                     " (without continuity correction)")
      }
    }
    else if (any(meth$method == "MH"))
      details <- paste0(details,"\n- Mantel-Haenszel method")
    
    ##
    ## Peto method
    ##
    if (any(meth$method == "Peto"))
      details <- paste0(details,"\n- Peto method")
    
    ##
    ## Cochran method
    ##
    if (any(meth$method == "Cochran"))
      details <- paste0(details,"\n- Cochran method")
        
  }

  if (overall | metabind | by) {
    
    ##
    ## Inverse variance method
    ##
    vars <-
      if (metacont)
        c("three.level", "pooledvar")
      else
        "three.level"
    ##      
    dat.iv <-
      unique(meth[meth$model %in% selmod & meth$method == "Inverse",
                  vars, drop = FALSE])
    ##
    for (i in seq_len(nrow(dat.iv))) {
      details <- paste0(
        details,
        "\n- Inverse variance method",
        if (dat.iv$three.level[i]) " (three-level model)" else "",
      if (metacont && !is.na(dat.iv$pooledvar[i]) && dat.iv$pooledvar[i])
        " (with pooled variance for individual studies)" else "")
    }
    
    ##
    ## SSW method
    ##
    if (any(meth$method[meth$model %in% selmod] == "SSW"))
      details <- paste0(details,"\n- Sample size method")
    
    ##
    ## GLMM
    ##
    
    if (metabin) {
      mgbin <-
        unique(meth$model.glmm[meth$model %in% selmod & meth$method == "GLMM"])
      for (mgbin.i in mgbin) {
        details <- paste0(
          details,
          if (mgbin.i == "UM.FS")
            "\n- Logistic regression model (fixed study effects)"
          else if (mgbin.i == "UM.RS")
            paste0(
              "\n- Mixed-effects logistic regression model ",
              "(random study effects)")
          else if (mgbin.i == "CM.EL")
            paste0(
              "\n- Generalised linear mixed model ",
              "(conditional Hypergeometric-Normal)")
          else if (mgbin.i == "CM.AL")
            "\n- Generalised linear mixed model (conditional Binomial-Normal)")
      }
    }
    
    if (metainc) {
      mginc <-
        unique(meth$model.glmm[meth$model %in% selmod & meth$method == "GLMM"])
      for (mginc.i in mginc) {
        details <- paste0(
          details,
          if (mginc.i == "UM.FS")
            "\n- Poisson regression model (fixed study effects)"
          else if (mginc.i == "UM.RS")
            paste0(
              "\n- Mixed-effects Poisson regression model ",
              "(random study effects)")
          else if (mginc.i == "CM.EL")
            paste0(
              "\n- Generalised linear mixed model ",
              "(conditional Poisson-Normal)"))
      }
    }
    
    if (metaprop & any(meth$method == "GLMM"))
      details <-
        paste0(details, "\n- Random intercept logistic regression model")
    
    if (metarate & any(meth$method == "GLMM"))
      details <-
        paste0(details, "\n- Random intercept Poisson regression model")
  }
  
  
  ##
  ##
  ## (3) Estimation of between-study variance
  ##
  ##
  
  if ((random & (overall | metabind | by)) | overall.hetstat) {
    dat.mt <-
      unique(meth[meth$model %in% "random", c("tau.preset", "method.tau")])
    ##
    tau.preset <- dat.mt$tau.preset[!is.na(dat.mt$tau.preset)]
    ##
    if (length(tau.preset) >= 1) {
      tau.preset <-
        formatPT(tau.preset, lab = TRUE, labval = text.tau,
                 digits = digits.tau, lab.NA = "NA", big.mark = big.mark)
      ##
      details <-
        paste0(details,
               "\n- Preset square root of between-study variance: ",
               cond(tau.preset, only.finite = FALSE))
    }
    ##
    method.tau <- unique(dat.mt$method.tau[is.na(dat.mt$tau.preset)])
    ##
    if (length(method.tau) >= 1) {
      for (mti in method.tau) {
        details <-
          paste0(details,
                 if (mti == "DL")
                   "\n- DerSimonian-Laird estimator"
                 ##
                 else if (mti == "PM")
                   "\n- Paule-Mandel estimator"
                 ##
                 else if (mti == "REML")
                   "\n- Restricted maximum-likelihood estimator"
                 ##
                 else if (mti == "ML")
                   "\n- Maximum-likelihood estimator"
                 ##
                 else if (mti == "HS")
                   "\n- Hunter-Schmidt estimator"
                 ##
                 else if (mti == "SJ")
                   "\n- Sidik-Jonkman estimator"
                 ##
                 else if (mti == "HE")
                   "\n- Hedges estimator"
                 ##
                 else if (mti == "EB")
                   "\n- Empirical Bayes estimator")
        ##
        if (mti != "")
          details <- paste(details, "for", text.tau2)
        ##
        if (replaceNULL(x$tau.common, FALSE))
          details <-
            paste0(details,
                   "\n  (assuming common ", text.tau2,
                   " in subgroups)")
      }
    }
    ##
    if (metabin) {
      dat.qc <-
        unique(meth[meth$model %in% "random" & is.na(meth$tau.preset),
                    c("method.tau", "Q.Cochrane")])
      if (any(dat.qc$method.tau == "DL" & dat.qc$Q.Cochrane))
        details <-
          paste0(details,
                 "\n- Mantel-Haenszel estimator used in ",
                 "calculation of Q and ", text.tau2,
                 if (options()$width > 70) " (like RevMan 5)"
                 else "\n  (like RevMan 5)")
    }
  }
  
  
  ##
  ##
  ## (4) Confidence interval for between-study variance
  ##
  ##

  if (print.tau.ci) {
    method.tau.ci <- unique(meth$method.tau.ci[meth$model %in% "random"])
    ##
    if (length(method.tau) >= 1) {
      for (mtci in method.tau.ci) {
        details <-
          paste0(details,
                 if (mtci == "QP")
                   "\n- Q-Profile method"
                 ##
                 else if (mtci == "BJ")
                   "\n- Biggerstaff and Jackson method"
                 ##
                 else if (mtci == "J")
                   "\n- Jackson method"
                 ##
                 else if (mtci == "PL")
                   "\n- Profile-Likelihood method")
        ##
        if (mtci != "")
          details <-
            paste(details,
                  "for confidence interval of", text.tau2, "and", text.tau)
      }
    }
  }
  
  
  ##
  ##
  ## (5) Confidence interval of random effects estimate
  ##
  ##

  if (random) {
    vars <- c("method", "method.random.ci", "df.random", "three.level")
    dat.rc <- unique(meth[meth$model == "random", vars])
    ##
    dat.rc.hk1 <-
      subset(dat.rc,
             dat.rc$method.random.ci == "HK" &
             !(dat.rc$method == "GLMM" | dat.rc$three.level))
    ##
    dat.rc.hk2 <-
      subset(dat.rc,
             dat.rc$method.random.ci == "HK" &
             (dat.rc$method == "GLMM" | dat.rc$three.level))
    ##
    dat.rc.ckr <- subset(dat.rc, dat.rc$method.random.ci == "classic-KR")
    dat.rc.kr <- subset(dat.rc, dat.rc$method.random.ci == "KR")
    ##
    more.ci <- sum(1L * (nrow(dat.rc.hk1) > 0) +
                   1L * (nrow(dat.rc.hk2) > 0) +
                   1L * (nrow(dat.rc.ckr) > 0) +
                   1L * (nrow(dat.rc.kr) > 0)) > 1
    ##
    if (nrow(dat.rc.hk1) > 0) {
      details <-
        paste0(
          details,
          "\n- Hartung-Knapp ",
          if (more.ci) "(HK) ",
          "adjustment for random effects model (df = ",
          cond(dat.rc.hk1$df.random),
          ")")
      ##
      if (any(dat.rc.hk1$adhoc.hakn.ci != ""))
        details <- paste0(
          details,
          "\n  (with ",
          if (any(dat.rc.hk1$adhoc.hakn.ci == ""))
            "and without ",
          "ad hoc correction)")
    }
    ##
    if (nrow(dat.rc.hk2) > 0) {
      details <-
        paste0(
          details,
          "\n- Random effects confidence interval based on t-distribution ",
          if (more.ci) "(T) ",
          "(df = ",
          cond(dat.rc.hk2$df.random),
          ")")
    }
    ##
    if (nrow(dat.rc.ckr) > 0)
      details <-
        paste0(
          details,
          "\n- Classic method instead of Kenward-Roger adjustment ",
          if (more.ci) "(classic-KR) ",
          "used for random effects model")
    ##
    if (nrow(dat.rc.kr) > 0)
      details <-
        paste0(
          details,
          "\n- Kenward-Roger ",
          if (more.ci) "(KR) ",
          "adjustment for random effects model (df = ",
          cond(dat.rc.kr$df.random),
          ")")
  }
  
  
  ##
  ##
  ## (6) Prediction interval
  ##
  ##

  if (prediction) {
    dat.pr <- unique(pred)
    ##
    dat.pr.hts1 <- subset(dat.pr, dat.pr$method.predict == "HTS")
    dat.pr.hk <- subset(dat.pr, dat.pr$method.predict == "HK")
    dat.pr.hts2 <- subset(dat.pr, dat.pr$method.predict == "HTS-KR")
    dat.pr.kr <- subset(dat.pr, dat.pr$ethod.predict == "KR")
    dat.pr.nnf <- subset(dat.pr, dat.pr$method.predict == "NNF")
    dat.pr.s <- subset(dat.pr, dat.pr$method.predict == "S")
    ##
    more.pi <- sum(1L * (nrow(dat.pr.hts1) > 0) +
                   1L * (nrow(dat.pr.hk) > 0) +
                   1L * (nrow(dat.pr.hts2) > 0) +
                   1L * (nrow(dat.pr.kr) > 0) +
                   1L * (nrow(dat.pr.nnf) > 0) +
                   1L * (nrow(dat.pr.s) > 0)) > 1
    ##
    if (nrow(dat.pr.hts1) > 0)
      details <-
        paste0(
          details,
          "\n- Prediction interval based on t-distribution ",
          if (more.pi) "(HTS) ",
          "(df = ",
          cond(dat.pr.hts1$df.predict),
          ")")
    ##
    if (nrow(dat.pr.hk) > 0) {
      details <-
        paste0(
          details,
          "\n- Hartung-Knapp ",
          if (more.pi) "(HK) ",
          "prediction interval (df = ",
          cond(dat.pr.hk$df.random),
          ")")
      ##
      if (any(dat.pr.hk$adhoc.hakn.ci != ""))
        details <- paste0(
          details,
          "\n  (with ",
          if (any(dat.pr.hk$adhoc.hakn.ci == ""))
            "and without ",
          "ad hoc correction)")
    }
    ##
    if (nrow(dat.pr.hts2) > 0)
      details <-
        paste0(
          details,
          "\n- Prediction interval based on t-distribution ",
          if (more.pi) "(HTS) ",
          "(df = ",
          cond(dat.pr.hts2$df.predict),
          ") instead of ",
          "Kenward-Roger adjustment")
    ##
    if (nrow(dat.pr.kr) > 0)
      details <-
        paste0(
          details,
          "\n- Kenward-Roger ",
          if (more.pi) "(KR) ",
          "prediction interval (df = ",
          cond(dat.pr.kr$df.random),
          ")")
    ##
    if (nrow(dat.pr.nnf) > 0)
      details <-
        paste0(
          details,
          "\n- Boot-strap prediction interval ",
          if (more.pi) "(NNF) ",
          "(df = ",
          cond(dat.pr.nnf$df.random),
          ")")
    ##
    if (nrow(dat.pr.s) > 0)
      details <-
        paste0(
          details,
          "\n- Prediction interval based on standard normal distribution",
          if (more.pi) "(S)")
  }
  
  
  ##
  ##
  ## (7) Trim-and-fill method
  ##
  ##

  if (overall & trimfill) {
    type <- unique(meth$type[meth$model %in% selmod])
    type <- type[type != ""]
    ##
    if (length(type) > 0) {
      details <-
        paste0(details,
               "\n- Trim-and-fill method to adjust for funnel plot asymmetry",
               if (length(type) > 1)
                 "\n  ("
               else
                 " (",
               paste0(type, "-estimator", collapse = ", "), ")")
    }
  }
  
  
  ##
  ##
  ## (8) metamiss
  ##
  ##

  if (metamiss) {
    method.miss <- x$method.miss
    IMOR.e <- x$IMOR.e
    IMOR.c <- x$IMOR.c
    ##
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
    details <- paste0(details, "\n- Imputation method: ", mmiss)
  }
  
  
  ##
  ##
  ## (9) Information on effect measure
  ##
  ##
  
  sm <- x$sm
  ##
  details <-
    paste0(details,
           if (sm == "ZCOR")
             "\n- Fisher's z transformation of correlations"
           else if (sm == "COR")
             "\n- Untransformed correlations"
           else if (sm == "PFT" | sm == "IRFT")
             "\n- Freeman-Tukey double arcsine transformation"
           else if (sm == "PAS")
             "\n- Arcsine transformation"
           else if (is_log_effect(sm))
             "\n- Log transformation"
           else if (sm == "PLOGIT")
             "\n- Logit transformation"
           else if (sm == "PRAW")
             "\n- Untransformed proportions"
           else if (sm == "IR")
             "\n- Untransformed rates"
           else if (sm == "IRS")
             "\n- Square root transformation"
           else if (sm == "MRAW")
             "\n- Untransformed (raw) means")
  ##
  if (!is.null(func.transf))
    details <-
      paste0(details, "\n- User-specified transformation: ", func.transf)
  ##
  if (backtransf & !is.null(func.backtransf))
    details <-
      paste0(details,
             "\n- User-specified back-transformation: ", func.backtransf)
  ##
  ## Add information on method for standardised mean difference
  ##
  if (!is.null(meth$method.smd)) {
    dat.ms <-
      unique(meth[meth$model %in% selmod,
                  c("method.smd", "exact.smd", "sd.glass")])
    ##
    for (i in seq_len(nrow(dat.ms))) {
      method.smd.i <- dat.ms$method.smd[i]
      exact.smd.i <- dat.ms$exact.smd[i]
      ##
      if (method.smd.i == "Hedges")
        details <-
          paste0(details,
                 "\n- Hedges' g (bias corrected standardised mean ",
                 "difference",
                 if (exact.smd.i)
                   "; using exact formulae)"
                 else
                   ")")
      else if (method.smd.i == "Cohen")
        details <-
          paste0(details,
                 "\n- Cohen's d (standardised mean difference",
                 if (exact.smd.i)
                   "; using exact formulae)"
                 else
                   ")")
      else if (method.smd.i == "Glass")
        details <-
          paste0(details,
                 "\n- Glass' delta (standardised mean difference; ",
                 if (dat.ms$sd.glass[i] == "control")
                   "based on control group)"
                 else
                   "based on experimental group)")
    }
  }
  
  
  ##
  ##
  ## (10) Information on confidence interval for individual studies
  ##
  ##
  
  info.ci <- attr(x, ".print.study.results.")
  details.ci <- ""
  ##
  if ((!is.null(info.ci) && info.ci) & !is.null(x$method.ci)) {
    if (any(x$k.all > 1)) {
      if (x$method.ci == "WSCC")
        fis <- "\n  for individual studies"
      else
        fis <- " for individual studies"
    }
    else
      fis <- ""
    ##
    details.ci <-
      if (x$method.ci == "CP")
        paste0("\n- Clopper-Pearson confidence interval", fis)
      else if (x$method.ci == "WS")
        paste0("\n- Wilson Score confidence interval", fis)
      else if (x$method.ci == "WSCC")
        paste0("\n- Wilson Score confidence interval with ",
               "continuity correction", fis)
      else if (x$method.ci == "AC")
        paste0("\n- Agresti-Coull confidence interval", fis)
      else if (x$method.ci == "SA")
        paste0("\n- Simple approximation confidence interval", fis)
      else if (x$method.ci == "SACC")
        paste0("\n- Simple approximation confidence interval with ",
               "continuity correction", fis)
      else if (x$method.ci == "NAsm")
        paste0("\n- Normal approximation confidence interval", fis)
      else if (x$method.ci == "Poisson")
        paste0("\n- Exact Poisson confidence interval", fis)
      else if (x$method.ci == "t")
        paste0("\n- Confidence interval", fis, " based on t-distribution")
  }
  ##
  details <- paste0(details, details.ci)
  
  
  ##
  ##
  ## (11) Information on continuity correction
  ##
  ##

  if (metabin | metainc | metaprop | metarate) {
    vars <- c("method", "incr", "method.incr", "sparse", "k.all")
    ##
    if (metabin)
      vars <- c(vars, "allstudies", "doublezeros", "MH.exact", "RR.Cochrane")
    ##
    dat.cc <- unique(meth[meth$model %in% selmod, vars])
    ##
    if (!(metabin & sm == "ASD")) {
      details.cc <- ""
      ##
      for (i in seq_len(nrow(dat.cc))) {
        method.incr.i <- dat.cc$method.incr[i]
        incr.i <- dat.cc$incr[i]
        sparse.i <- dat.cc$sparse[i]
        k.all.i <- dat.cc$k.all[i]
        ##
        details.rr <- NULL
        ##
        if (metabin) {
          txtCC.ind.i <-
            (dat.cc$method[i] == "MH" & dat.cc$MH.exact[i]) |
            dat.cc$method[i] == "GLMM"
          ##
          if (dat.cc$RR.Cochrane[i] &
              (method.incr.i == "all" | (sparse.i & incr.i > 0))) {
            details.rr <-
              if (options()$width <= 70)
                " (applied twice to sample sizes, like RevMan 5)"
              else
                "\n  (applied twice to sample sizes, like RevMan 5)"
          }
        }
        else
          txtCC.ind.i <- dat.cc$method[i] == "GLMM"
        ##
        if (method.incr.i == "all") {
          if (incr.i == "TACC") {
            details.cc <- c(
              details.cc,
              if (k.all.i > 1)
                "\n- Treatment arm continuity correction in all studies"
              else
                "\n- Treatment arm continuity correction")
          }
          else if (as.numeric(incr.i) > 0)
            details.cc <- c(
              details.cc,
              paste0("\n- Continuity correction of ",
                     round(as.numeric(incr.i), 4),
                     if (k.all.i > 1)
                       " in all studies",
                     details.rr))
        }
        else if (sparse.i) {
          if (incr.i == "TACC") {
            details.cc <- c(
              details.cc,
              paste0(
                "\n- Treatment arm continuity correction in ",
                if (k.all.i > 1)
                  "studies "
                else
                  "study ",
                "with ",
                if (options()$width > 70)
                  " "
                else
                  "\n  ",
                "zero cell frequencies",
                details.rr
              ))
          }
          else if (as.numeric(incr.i) > 0) {
            details.cc <- c(
              details.cc,
              paste0(
                "\n- Continuity correction of ",
                round(as.numeric(incr.i), 4),
                if (k.all.i > 1)
                  " in studies with",
                if (k.all.i > 1 & options()$width > 70)
                  " "
                else if (k.all.i > 1)
                  "\n  ",
                if (k.all.i > 1)
                  "zero cell frequencies",
                details.rr))
          }
          ##
          if ((incr.i == "TACC" || as.numeric(incr.i) > 0) && txtCC.ind.i)
            details.cc <- c(
              details.cc,
              "\n  (only used to calculate individual study results)")
        }
        ##
        if (metabin) {
          if (dat.cc$allstudies[i] & dat.cc$doublezeros[i])
            details.cc <- c(
              details.cc,
              if (k.all.i > 1)
                "\n- Studies with double zeros included in meta-analysis"
              else
                "\n- Study with double zeros considered")
        }
      }
      ##
      details.cc <- paste0(unique(details.cc), collapse = "")
      details <- paste0(details, details.cc)
    }
  }
  
  
  ##
  ##
  ## (12) Information on number of events and null hypothesis
  ##
  ##
  
  pscale <- x$pscale
  irscale <- x$irscale
  irunit <- x$irunit
  ##
  if (!is.null(pscale) && pscale != 1)
    details <-
      paste0(details,
             "\n- Events per ",
             format(pscale, scientific = FALSE, big.mark = big.mark),
             " observations")
  ##
  if (!is.null(irscale) && irscale != 1)
    details <-
      paste0(details,
             "\n- Events per ",
             format(irscale, scientific = FALSE, big.mark = big.mark),
             " ", irunit)
  ##
  null.effect <- x$null.effect
  ##
  if (!is.na(null.effect) && null.effect != 0) {
    details <- paste0(details, "\n- Null hypothesis: effect is equal to ")
    ##
    if (!is.null(pscale) && pscale != 1)
      details <-
        paste0(details,
               format(round(null.effect * pscale, digits), scientific = FALSE,
                      big.mark = big.mark),
               " events per ",
               format(pscale, scientific = FALSE, big.mark = big.mark),
               " observations")
    else if (!is.null(irscale) && irscale != 1)
      details <-
        paste0(details,
               format(round(null.effect * irscale, digits),
                      scientific = FALSE, big.mark = big.mark),
               " events per ",
               format(irscale, scientific = FALSE, big.mark = big.mark),
               " ", irunit)
    else
      details <- paste0(details,
                        format(null.effect, scientific = FALSE,
                               big.mark = big.mark))
  }
  
  
  if (details != "")
    cat(paste0("\nDetails",
               if (any(x$k.all > 1)) " on meta-analytical method",
               ":", details, "\n"))
  
  invisible(NULL)
}


NULL
