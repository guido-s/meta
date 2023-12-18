catmeth <- function(x,
                    common, random, prediction, overall, overall.hetstat,
                    func.transf, backtransf, func.backtransf,
                    big.mark, digits, digits.tau, text.tau, text.tau2,
                    print.tau2 = TRUE, print.tau2.ci = FALSE,
                    print.tau = FALSE, print.tau.ci = FALSE, 
                    forest = FALSE,
                    print.df = TRUE
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
  bothmod <- length(selmod) > 1
  ##
  details <- NULL
  ##
  width <- options()$width
  ##
  method <-
    if (metacont)
      "metacont"
    else if (metabin)
      "metabin"
    else if (metainc)
      "metainc"
    else if (metaprop)
      "metaprop"
    else if (metarate)
      "metarate"
  ##
  if (forest) {
    text.tau2 <- "tau^2"
    text.tau <- "tau"
  }
  ##
  text.t <- ""
  ##
  if (print.tau2 | metabind)
    text.t <- text.tau2
  else if (print.tau)
    text.t <- text.tau
  
  
  ##
  ##
  ## (2) Meta-analysis methods
  ##
  ##

  if (overall | metabind | by) {
    
    meth.ma <- meth[meth$model %in% selmod, , drop = FALSE]
    ##
    vars.ma <- c("model", "method")
    ##
    if ((metabin | metainc) & any(meth.ma$method == "GLMM"))
      vars.ma <- c(vars.ma, "model.glmm")
    ##
    vars.ma <- c(vars.ma, "three.level", "rho")
    ##
    if (metacont)
      vars.ma <- c(vars.ma, "pooledvar")
    else if (metabin)
      vars.ma <- c(vars.ma, "incr", "method.incr", "sparse", "MH.exact")
    else if (metainc)
      vars.ma <- c(vars.ma, "incr", "method.incr", "sparse")
    ##
    meth.ma <-
      unique(meth.ma[, vars.ma, drop = FALSE])
    ##
    details.i <- vector("character", length = nrow(meth.ma))
    for (i in seq_len(nrow(meth.ma)))
      details.i[i] <- methtxt(meth.ma, i, random, method)
    ##
    details <- paste(c(details, unique(details.i)), collapse = "")
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
    if (length(tau.preset) >= 1)
      details <-
        paste0(details,
               "\n- Preset square root of between-study variance: ",
               text.tau, " = ",
               cond(tau.preset, only.finite = FALSE, digits = digits.tau,
                    big.mark = big.mark))
    ##
    method.tau <- unique(dat.mt$method.tau[is.na(dat.mt$tau.preset)])
    ##
    if ((print.tau2 | print.tau | metabind) & length(method.tau) >= 1) {
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
          details <- paste(details, "for", text.t)
        ##
        if (replaceNULL(x$tau.common, FALSE))
          details <-
            paste0(details,
                   if (forest) " " else "\n  ",
                   "(assuming common ", text.t, " in subgroups)")
      }
    }
    ##
    if (metabin) {
      dat.qc <-
        unique(meth[meth$model %in% "random" & is.na(meth$tau.preset),
                    c("method.tau", "Q.Cochrane")])
      dat.qc$Q.Cochrane <- replaceNA(dat.qc$Q.Cochrane, FALSE)
      ##
      if (any(dat.qc$method.tau == "DL" & dat.qc$Q.Cochrane))
        details <-
          paste0(details,
                 "\n- Mantel-Haenszel estimator used in ",
                 "calculation of Q and ", text.t,
                 if (width > 70 | forest) " (like RevMan 5)"
                 else "\n  (like RevMan 5)")
    }
  }
  
  
  ##
  ##
  ## (4) Confidence interval for between-study variance
  ##
  ##

  if (print.tau.ci | print.tau2.ci) {
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
            paste0(details,
                  " for confidence interval of ",
                  if (print.tau2 & print.tau)
                    paste(text.tau2, "and", text.tau)
                  else if (print.tau2)
                    text.tau2
                  else
                    text.tau)
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
          "adjustment for random effects model",
          if (print.df)
            paste0(" (df = ", cond(dat.rc.hk1$df.random, digits = 0), ")")
        )
      ##
      if (any(dat.rc.hk1$adhoc.hakn.ci != ""))
        details <- paste0(
          details,
          if (forest) " " else "\n  ", "(with ",
          if (any(dat.rc.hk1$adhoc.hakn.ci == ""))
            "and without ",
          "ad hoc correction)")
    }
    ##
    if (nrow(dat.rc.hk2) > 0) {
      details <-
        paste0(
          details,
          "\n- Random effects confidence interval based on t-distribution",
          if (more.ci) " (T)",
          if (print.df)
            paste0(" (df = ", cond(dat.rc.hk2$df.random, digits = 0), ")")
        )
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
          "adjustment for random effects model",
          if (print.df)
            paste0(" (df = ", cond(dat.rc.kr$df.random), ")")
        )
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
    dat.pr.kr <- subset(dat.pr, dat.pr$method.predict == "KR")
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
          cond(dat.pr.hts1$df.predict, digits = 0),
          ")")
    ##
    if (nrow(dat.pr.hk) > 0) {
      details <-
        paste0(
          details,
          "\n- Hartung-Knapp ",
          if (more.pi) "(HK) ",
          "prediction interval (df = ",
          cond(dat.pr.hk$df.predict, digits = 0),
          ")")
      ##
      if (any(dat.pr.hk$adhoc.hakn.ci != ""))
        details <- paste0(
          details,
          if (forest) " " else "\n  ", "(with ",
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
          cond(dat.pr.hts2$df.predict, digits = 0),
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
          cond(dat.pr.kr$df.predict),
          ")")
    ##
    if (nrow(dat.pr.nnf) > 0)
      details <-
        paste0(
          details,
          "\n- Boot-strap prediction interval ",
          if (more.pi) "(NNF) ",
          "(df = ",
          cond(dat.pr.nnf$df.predict, digits = 0),
          ")")
    ##
    if (nrow(dat.pr.s) > 0)
      details <-
        paste0(
          details,
          "\n- Prediction interval based on standard normal distribution",
          if (more.pi) " (S)")
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
               if (length(type) == 1 | forest)
                 " ("
               else
                 "\n  (",
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
          paste0(if (forest) " " else "\n  ", mmiss,
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
                   paste0(if (width >= 80 | forest) "; " else ";\n  ",
                          "using exact formulae)")
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
                   "based on second group)"
                 else
                   "based on first group)")
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
        fis <- paste0(if (forest) " " else "\n  ",
                      "for individual studies")
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
        if (!is.na(dat.cc$incr[i]))
          incr.i <- dat.cc$incr[i]
        else
          incr.i <- 0
        sparse.i <- !is.na(dat.cc$sparse[i]) && dat.cc$sparse[i]
        k.all.i <- dat.cc$k.all[i]
        ##
        details.rr <- NULL
        ##
        if (metabin) {
          txtCC.ind.i <-
            (dat.cc$method[i] == "MH" & dat.cc$MH.exact[i]) |
            dat.cc$method[i] == "GLMM"
          ##
          if (!is.na(dat.cc$RR.Cochrane[i]) &&
              dat.cc$RR.Cochrane[i] &
              (method.incr.i == "all" | (sparse.i & incr.i > 0))) {
            details.rr <-
              if (width >= 70 | forest)
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
                "with",
                if (width > 70 | forest)
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
                if (k.all.i > 1 & width > 70)
                  " "
                else if (k.all.i > 1 & !forest)
                  "\n  ",
                if (k.all.i > 1)
                  "zero cell frequencies",
                details.rr))
          }
          ##
          if ((incr.i == "TACC" || as.numeric(incr.i) > 0) && txtCC.ind.i)
            details.cc <- c(
              details.cc,
              if (forest) " " else "\n  ",
              "(only used to calculate individual study results)")
        }
        ##
        if (metabin) {
          if ((!is.na(dat.cc$allstudies[i]) && dat.cc$allstudies[i]) &
              (!is.na(dat.cc$doublezeros[i]) && dat.cc$doublezeros[i]))
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
    if (is_prop(sm) || sm == "RD")
      details <-
        paste0(details,
               "\n- Events per ",
               format(pscale, scientific = FALSE, big.mark = big.mark),
               " observations")
    else
      details <-
        paste0(details,
               "\n- Scaling factor for results: ",
               format(pscale, scientific = FALSE, big.mark = big.mark))   
  ##
  if (!is.null(irscale) && irscale != 1)
    if (is_rate(sm) || sm == "IRD")
      details <-
        paste0(details,
               "\n- Events per ",
               format(irscale, scientific = FALSE, big.mark = big.mark),
               " ", irunit)
    else
      details <-
        paste0(details,
               "\n- Scaling factor for results: ",
               format(irscale, scientific = FALSE, big.mark = big.mark))   
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
  
  
  if (!is.null(details) && length(details) > 0 && details != "") {
    details <-
      paste0("\nDetails",
             if ((common | random | prediction) && any(x$k.all > 1))
               " on meta-analytical method",
             ":", details)
    ##
    if (!forest)
      cat(paste0(details, "\n"))
  }
  
  invisible(details)
}


NULL
