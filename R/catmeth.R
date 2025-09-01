catmeth <- function(x,
                    common, random, prediction, overall, overall.hetstat,
                    prediction.subgroup = FALSE,
                    #
                    func.transf, backtransf, func.backtransf,
                    #
                    big.mark, digits, digits.tau, text.tau, text.tau2,
                    #
                    print.tau2 = TRUE, print.tau2.ci = FALSE,
                    print.tau = FALSE, print.tau.ci = FALSE,
                    #
                    print.I2 = FALSE, text.I2,
                    #
                    print.df = TRUE,
                    #
                    forest = FALSE) {

  ##
  ##
  ## (1) Some settings
  ##
  ##
  
  metacum.metainf <- inherits(x, "metacum") | inherits(x, "metainf")
  metabind <- inherits(x, "metabind")
  #
  metabin  <- inherits(x, "metabin") & !metacum.metainf
  metacont <- inherits(x, "metacont") & !metacum.metainf
  metacor  <- inherits(x, "metacor") & !metacum.metainf
  metagen  <- inherits(x, "metagen") & !metacum.metainf
  metainc  <- inherits(x, "metainc") & !metacum.metainf
  metamean <- inherits(x, "metamean") & !metacum.metainf
  metaprop <- inherits(x, "metaprop") & !metacum.metainf
  metarate <- inherits(x, "metarate") & !metacum.metainf
  #
  trimfill <- inherits(x, "trimfill") & !metacum.metainf
  metamiss <- inherits(x, "metamiss") & !metacum.metainf
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
  #
  # Print information on meta-analysis of n-of-1-trials first
  #
  if (!is.null(x$cycles)) {
    details <-
      paste0("\n- ",
             if (any(x$k.all > 1)) "Meta-a" else "A",
             "nalysis of n-of-1 trials (pooled SD=",
             x$sd.n_of_1,
             ")")
  }
  else
    details <- NULL
  #
  width <- options()$width
  ##
  if (metabin)
    method <- "metabin"
  else if (metacont)
    method <- "metacont"
  else if (metacor)
    method <- "metacor"
  else if (metagen)
    method <- "metagen"
  else if (metainc)
    method <- "metainc"
  else if (metamean)
    method <- "metamean"
  else if (metaprop)
    method <- "metaprop"
  else if (metarate)
    method <- "metarate"
  else if (metacum.metainf | metabind)
    method <- x$classes
  ##
  if (forest) {
    text.tau2 <- "tau^2"
    text.tau <- "tau"
    #
    text.I2 <- "I^2"
  }
  ##
  text.t <- ""
  ##
  if (print.tau2 | metabind | metacum.metainf)
    text.t <- text.tau2
  else if (print.tau)
    text.t <- text.tau
  #
  method.I2 <- replaceNULL(x$method.I2, "Q")
  
  
  ##
  ##
  ## (2) Meta-analysis methods
  ##
  ##
    
  if (overall | metabind | metacum.metainf | by) {
    
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
    else if (metabin) {
      vars.ma <- c(vars.ma, "incr", "method.incr", "sparse", "MH.exact")
      if (any(meth.ma$method == "LRP"))
        vars.ma <- c(vars.ma, "phi")
    }
    else if (metainc)
      vars.ma <- c(vars.ma, "incr", "method.incr", "sparse")
    #
    meth.ma <- unique(meth.ma[, vars.ma, drop = FALSE])
    ##
    details.i <- vector("character", length = nrow(meth.ma))
    for (i in seq_len(nrow(meth.ma)))
      details.i[i] <- text_meth(meth.ma, i, random, method)
    ##
    details <- paste(c(details, unique(details.i)), collapse = "")
    #
    # Information on user-specified weights
    #
    if (!is.null(x$weights.common) | !is.null(x$weights.random)) {
      if (!is.null(x$weights.common) & !is.null(x$weights.random))
        details.usw <- ""
      else if  (!is.null(x$weights.common))
        details.usw <- " (common effect model)"
      else
        details.usw <- " (random effects model)"
      #
      details <-
        paste0(details, "\n- User-specified weights", details.usw)
    }
  }
  
  
  ##
  ##
  ## (3) Estimation of between-study variance
  ##
  ##
  if (((random | (metacum.metainf & (print.tau2 | print.tau))) & 
       (overall | metabind | by)) | overall.hetstat) {
    if (metacum.metainf)
      sel.mt <- meth$model %in% c("common", "random")
    else
      sel.mt <- meth$model == "random"
    #
    dat.mt <- unique(meth[sel.mt, c("tau.preset", "method.tau")])
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
    if (metabin | metacum.metainf) {
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
  ## (5) Method to estimate I2
  ##
  ##
    
  if (print.I2) {
    method.I2 <- unique(meth$method.I2)
    ##
    if (length(method.I2) >= 1) {
      for (mi2i in method.I2) {
        details <-
          paste0(details,
                 if (mi2i == "Q")
                   paste0("\n- Calculation of ", text.I2,
                          " based on Q")
                 else if (mi2i == "tau2")
                   paste0("\n- Calculation of ", text.I2,
                          " based on ", text.tau2))
      }
    }
  }
  
  
  #
  #
  # (6) Confidence interval of common effect estimate
  #
  #
  
  if (common &
      any(meth$model == "common" & meth$method.common.ci == "IVhet"))
    details <- paste0(details,
                      "\n- Inverse variance heterogeneity model")
  
  
  ##
  ##
  ## (7) Confidence interval of random effects estimate
  ##
  ##

  if (random) {
    vars <- c("method", "method.random.ci", "df.random", "three.level")
    #
    dat.rc <- unique(meth[meth$model == "random", vars])
    dat.rc.hk <-
      unique(meth[meth$model == "random", c(vars, "adhoc.hakn.ci")])
    ##
    dat.rc.hk <-
      subset(dat.rc.hk,
             dat.rc.hk$method.random.ci == "HK" &
             !(dat.rc.hk$method %in% c("GLMM", "LRP") | dat.rc.hk$three.level))
    ##
    dat.rc.hk.tdist <-
      subset(dat.rc,
             dat.rc$method.random.ci == "HK" &
             (dat.rc$method %in% c("GLMM", "LRP") | dat.rc$three.level))
    ##
    dat.rc.ckr <- subset(dat.rc, dat.rc$method.random.ci == "classic-KR")
    dat.rc.kr <- subset(dat.rc, dat.rc$method.random.ci == "KR")
    ##
    more.ci <- sum(1L * (nrow(dat.rc.hk) > 0) +
                   1L * (nrow(dat.rc.hk.tdist) > 0) +
                   1L * (nrow(dat.rc.ckr) > 0) +
                   1L * (nrow(dat.rc.kr) > 0)) > 1
    ##
    if (nrow(dat.rc.hk) > 0) {
      details <-
        paste0(
          details,
          "\n- Hartung-Knapp ",
          if (more.ci) "(HK) ",
          "adjustment for random effects model",
          if (print.df)
            paste0(" (df = ", cond(dat.rc.hk$df.random, digits = 0), ")")
        )
      ##
      if (any(dat.rc.hk$adhoc.hakn.ci != ""))
        details <- paste0(
          details,
          if (forest) " " else "\n  ", "(with ",
          if (any(dat.rc.hk$adhoc.hakn.ci == ""))
            "and without ",
          "ad hoc correction)")
    }
    ##
    if (nrow(dat.rc.hk.tdist) > 0) {
      details <-
        paste0(
          details,
          "\n- Random effects confidence interval based on t-distribution",
          if (more.ci) " (T)",
          if (print.df)
            paste0(" (df = ", cond(dat.rc.hk.tdist$df.random, digits = 0), ")")
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
  ## (8) Prediction interval
  ##
  ##

  more.pi <- NULL
  #
  if (!prediction & replaceNULL(prediction.subgroup, FALSE)) {
    prediction <- TRUE
    print.df <- FALSE
    more.pi <- TRUE
  }
  #
  if (prediction) {
    dat.pr <- unique(pred)
    ##
    dat.pr.v <- subset(dat.pr, dat.pr$method.predict == "V")
    dat.pr.v.kr <- subset(dat.pr, dat.pr$method.predict == "V-KR")
    dat.pr.hts1 <- subset(dat.pr, dat.pr$method.predict == "HTS")
    dat.pr.hk <- subset(dat.pr, dat.pr$method.predict == "HK")
    dat.pr.hts2 <- subset(dat.pr, dat.pr$method.predict == "HTS-KR")
    dat.pr.kr <- subset(dat.pr, dat.pr$method.predict == "KR")
    dat.pr.nnf <- subset(dat.pr, dat.pr$method.predict == "NNF")
    dat.pr.s <- subset(dat.pr, dat.pr$method.predict == "S")
    ##
    if (is.null(more.pi))
      more.pi <- sum(1L * (nrow(dat.pr.v) > 0) +
                       1L * (nrow(dat.pr.v.kr) > 0) +
                       1L * (nrow(dat.pr.hts1) > 0) +
                       1L * (nrow(dat.pr.hk) > 0) +
                       1L * (nrow(dat.pr.hts2) > 0) +
                       1L * (nrow(dat.pr.kr) > 0) +
                       1L * (nrow(dat.pr.nnf) > 0) +
                       1L * (nrow(dat.pr.s) > 0)) > 1
    #
    if (nrow(dat.pr.v) > 0)
      details <-
        paste0(
          details,
          "\n- Prediction interval based on t-distribution",
          if (more.pi) " (V)",
          if (print.df)
            paste0(" (df = ",
                   cond(dat.pr.v$df.predict, digits = 0),
                   ")")
          )
    #
    if (nrow(dat.pr.v.kr) > 0)
      details <-
      paste0(
        details,
        "\n- Prediction interval based on t-distribution",
        if (more.pi) " (V)",
        if (print.df)
          paste0(" (df = ",
                 cond(dat.pr.v.kr$df.predict, digits = 0),
                 ")"),
        " instead of Kenward-Roger adjustment")
    #
    if (nrow(dat.pr.hts1) > 0)
      details <-
      paste0(
        details,
        "\n- Prediction interval based on t-distribution",
        if (more.pi) " (HTS)",
        if (print.df)
          paste0(" (df = ",
                 cond(dat.pr.hts1$df.predict, digits = 0),
                 ")")
        )
    ##
    if (nrow(dat.pr.hk) > 0) {
      details <-
        paste0(
          details,
          "\n- Hartung-Knapp",
          if (more.pi) " (HK)",
          " prediction interval",
          if (print.df)
            paste0(" (df = ",
                   cond(dat.pr.hk$df.predict, digits = 0),
                   ")")
        )
      ##
      if (any(dat.pr.hk$adhoc.hakn.pi != ""))
        details <- paste0(
          details,
          if (forest) " " else "\n  ", "(with ",
          if (any(dat.pr.hk$adhoc.hakn.pi == ""))
            "and without ",
          "ad hoc correction)")
    }
    ##
    if (nrow(dat.pr.hts2) > 0)
      details <-
        paste0(
          details,
          "\n- Prediction interval based on t-distribution",
          if (more.pi) " (HTS) ",
          if (print.df)
            paste0(" (df = ",
                   cond(dat.pr.hts2$df.predict, digits = 0),
                   ")"),
          " instead of Kenward-Roger adjustment")
    ##
    if (nrow(dat.pr.kr) > 0)
      details <-
        paste0(
          details,
          "\n- Kenward-Roger",
          if (more.pi) " (KR)",
          " prediction interval",
          if (print.df)
            paste0(" (df = ",
                   cond(dat.pr.kr$df.predict),
                   ")")
        )
    ##
    if (nrow(dat.pr.nnf) > 0)
      details <-
        paste0(
          details,
          "\n- Boot-strap prediction interval",
          if (more.pi) " (NNF)",
          if (print.df)
            paste0(" (df = ",
                   cond(dat.pr.nnf$df.predict, digits = 0),
                   ")")
        )
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
  ## (9) Trim-and-fill method
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
  ## (10) R function metamiss()
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
  ## (11) Information on effect measure
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
  ## (12) Information on confidence interval for individual studies
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
  ## (13) Information on continuity correction
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
        #
        sparse.i <- !is.na(dat.cc$sparse[i]) && dat.cc$sparse[i]
        k.all.i <- dat.cc$k.all[i]
        ##
        details.rr <- NULL
        ##
        if (metabin) {
          txtCC.ind.i <-
            (dat.cc$method[i] == "MH" & dat.cc$MH.exact[i]) |
            dat.cc$method[i] %in% c("GLMM", "LRP")
          ##
          if (!is.na(dat.cc$RR.Cochrane[i]) &&
              dat.cc$RR.Cochrane[i] &
              (method.incr.i == "all" | (sparse.i & incr.i != 0))) {
            details.rr <-
              if (width >= 70 | forest)
                " (applied twice to sample sizes, like RevMan 5)"
              else
                "\n  (applied twice to sample sizes, like RevMan 5)"
          }
        }
        else
          txtCC.ind.i <- dat.cc$method[i] %in% c("GLMM", "LRP")
        #
        # No continuity correction for generalized linear mixed model and
        # argument 'method.ci != "NAsm"'
        #
        if (!(txtCC.ind.i &
              (!is.null(x$method.ci) && x$method.ci != "NAsm"))) {
          if (method.incr.i == "user") {
            if (incr.i != 0)
              details.cc <- "\n- User-defined continuity correction"
          }
          else  if (method.incr.i == "all") {
            if (incr.i == "TACC") {
              details.cc <- c(
                details.cc,
                if (k.all.i > 1)
                  "\n- Treatment arm continuity correction in all studies"
                else
                  "\n- Treatment arm continuity correction")
            }
            else if (incr.i != 0)
              details.cc <- c(
                details.cc,
                paste0("\n- Continuity correction of ",
                       incr.i,
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
            else if (incr.i != 0) {
              details.cc <- c(
                details.cc,
                paste0(
                  "\n- Continuity correction of ",
                  incr.i,
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
            if ((incr.i == "TACC" || incr.i != 0) && txtCC.ind.i)
              details.cc <- c(
                details.cc,
                if (forest) " " else "\n  ",
                "(only used to calculate individual study results)")
          }
        }
        #
        if (metabin & method.incr.i != "user") {
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
  ## (14) Information on number of events and null hypothesis
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
  if (!is.na(null.effect) && (null.effect != 0 | metaprop | metarate)) {
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
               " of meta-analysis methods",
             ":", details)
    ##
    if (!forest)
      cat(paste0(details, "\n"))
  }
  
  invisible(details)
}


NULL
