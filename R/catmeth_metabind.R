catmeth_metabind <- function(x,
                             #
                             common, random, prediction,
                             #
                             print.tau2,
                             print.tau2.ci,
                             print.tau,
                             print.tau.ci,
                             print.I2,
                             print.pval,
                             #
                             backtransf,
                             #
                             digits,
                             digits.tau2,
                             digits.tau,
                             #
                             pscale,
                             irscale, irunit,
                             #
                             big.mark,
                             #
                             forest) {
  # Get rid of warning 'no visible binding for global variable'
  type <- method <- df <- id <- name <- outclab <- title <- label.right <-
    Q <- upper.Rb <- TE <- upper <- level.predict <- text.t <-
    method.interval <- method.incr <- NULL
  
  #
  #
  # (1) Check and set
  #
  #
  
  chkclass(x, "metabind")
  #
  print.df <- TRUE
  #
  if (!(common || random || prediction))
    return(invisible(""))
  #
  width <- options()$width
  #
  if (forest) {
    text.tau2 <- "tau^2"
    text.tau <- "tau"
    #
    text.I2 <- "I^2"
  }
  else {
    text.tau2 <- gs("text.tau2")
    text.tau <- gs("text.tau")
    #
    text.I2 <- gs("text.I2")
  }
  #
  results <- NULL
  if (common || random) {
    results <- as.data.frame(x, common = common, random = random,
                             prediction = FALSE, subgroup = FALSE)
    #
    results %<>%
      mutate(
        type =
          case_match(
            type,
            "Common effect model" ~ "common",
            "Random effects model" ~ "random",
            .default = type))
    #
    results$n <- NULL
    results$n.e <- NULL
    results$n.c <- NULL
    #
    results$event <- NULL
    results$event.e <- NULL
    results$event.c <- NULL
    #
    results$time <- NULL
    results$time.e <- NULL
    results$time.c <- NULL
    #
    results$tau.common <- replaceNULL(results$tau.common)
    results$Q.Cochrane <- replaceNULL(results$Q.Cochrane)
    #
    common <- any(results$type == "common")
    random <- any(results$type == "random")
  }
  common_or_random <- common || random
  #
  print.tau2 <- print.tau2 && random
  print.tau2.ci <- print.tau2.ci && print.tau2
  print.tau <- print.tau && random
  print.tau.ci <- print.tau.ci && print.tau
  #
  predict <- if (prediction) x$predict else NULL
  #
  ids <- unique(c(results$id, predict$id))
  overall <- x$overall %>%
    filter(id %in% ids) %>%
    select(-name, -outclab, -(title:label.right), -(Q:upper.Rb))
  #
  overall$n <- NULL
  overall$n.e <- NULL
  overall$n.c <- NULL
  #
  overall$event <- NULL
  overall$event.e <- NULL
  overall$event.c <- NULL
  #
  overall$time <- NULL
  overall$time.e <- NULL
  overall$time.c <- NULL
  #
  overall$tau.common <- replaceNULL(overall$tau.common)
  overall$Q.Cochrane <- replaceNULL(overall$Q.Cochrane)
  #
  if (nrow(overall) == 0)
    return(invisible(""))
  #
  if (common_or_random) {
    meth <- results %>%
      select(-id, -text, -(TE:upper), -outclab,
             -name, -(title:label.right), -(Q:upper.Rb),
             -common, -random, -prediction, -level.predict) %>%
      unique()
    #
    meth$statistic <- NULL
    meth$pval <- NULL
    #
    meth$tau.preset <- replaceNULL(meth$tau.preset)
    #
    sel.random <- meth$type == "random"
  }
  else {
    meth <- left_join(predict, overall, by = "id")
    sel.random <- rep(FALSE, nrow(meth))
  }
  
  #
  #
  # (2) Analysis types
  #
  #
  
  details <- character()
  #
  for (cl in unique(overall$class)) {
    description <- switch(
      cl,
      metabin = "comparing two groups with binary outcomes",
      metacont = "comparing two groups with continuous outcomes",
      metainc = "comparing two groups with incidence rates",
      metacor = "of single correlations",
      metamean = "of single means",
      metaprop = "of single proportions",
      metarate = "of single rates",
      metagen = "using the generic inverse variance method",
      NULL)
    #
    if (!is.null(description))
      details <- c(details, paste("\n- Meta-analysis", description))
  }
  #
  # Estimation methods, retaining component-specific metadata
  #
  if (common_or_random) {
    any.random <- any(sel.random)
    #
    meth.cr <- meth %>% mutate(model = type)
    #
    for (i in seq_len(nrow(meth.cr)))
      details <- c(details, text_meth(meth.cr, i, any.random, meth.cr$class[i]))
  }
  
  #
  #
  # (3) Estimation of between-study variance
  #
  #
  
  if (print.tau2 || print.tau) {
    dat.preset <- unique(meth[sel.random, c("tau.preset", "method.tau")])
    #
    # Information on tau.preset
    #
    tau.preset <- dat.preset$tau.preset[!is.na(dat.preset$tau.preset)]
    #
    if (length(tau.preset) >= 1) {
      details.preset <-
        paste0("\n- Preset square root of between-study variance: ",
               text.tau, " = ",
               cond(tau.preset, only.finite = FALSE, digits = digits.tau,
                    big.mark = big.mark))
    }
    else
      details.preset <- NULL
    #
    details <- c(details, details.preset)
    #
    # Information on method to estimate tau2
    #
    dat.tau <- unique(meth[sel.random, c("method.tau", "tau.common")])
    #
    details.tau <- NULL
    #
    if (nrow(dat.tau) >= 1) {
      for (i in seq_len(nrow(dat.tau))) {
        mti <- dat.tau$method.tau[i]
        #
        if (mti == "DL")
          details.tau.i <- "\n- DerSimonian-Laird estimator"
        #
        else if (mti == "PM")
          details.tau.i <- "\n- Paule-Mandel estimator"
        #
        else if (mti == "REML")
          details.tau.i <- "\n- Restricted maximum-likelihood estimator"
        #
        else if (mti == "ML")
          details.tau.i <- "\n- Maximum-likelihood estimator"
        #
        else if (mti == "HS")
          details.tau.i <- "\n- Hunter-Schmidt estimator"
        #
        else if (mti == "SJ")
          details.tau.i <- "\n- Sidik-Jonkman estimator"
        #
        else if (mti == "HE")
          details.tau.i <- "\n- Hedges estimator"
        #
        else if (mti == "EB")
          details.tau.i <- "\n- Empirical Bayes estimator"
        #
        if (mti != "")
          details.tau.i <- paste(details.tau.i, "for", text.tau2)
        #
        if (replaceNA(meth$tau.common[i], FALSE))
          details.tau.i <-
            paste0(details.tau.i,
                   if (forest) " " else "\n  ",
                   "(assuming common ", text.t, " in subgroups)")
        #
        details.tau <- c(details.tau, details.tau.i)
      }
    }
    else
      details.tau <- NULL
    #
    details <- c(details, details.tau)
    #
    if ("metabin" %in% overall$class) {
      dat.qc <-
        unique(meth[sel.random & is.na(meth$tau.preset),
                    c("method.tau", "Q.Cochrane")])
      dat.qc$Q.Cochrane <- replaceNA(dat.qc$Q.Cochrane, FALSE)
      #
      if (any(dat.qc$method.tau == "DL" & dat.qc$Q.Cochrane))
        details <-
        c(details,
          paste0("\n- Mantel-Haenszel estimator used in ",
                 "calculation of Q and ", text.tau2,
                 if (width > 70 | forest) " (like RevMan 5)"
                 else "\n  (like RevMan 5)"))
    }
  }
  
  #
  #
  # (4) Confidence interval for between-study variance
  #
  #
  
  details.tau.ci <- NULL
  #
  if (print.tau2.ci || print.tau.ci) {
    method.tau.ci <- unique(meth$method.tau.ci[sel.random])
    #
    for (mtci in method.tau.ci) {
      details.tau.ci.i <- switch(
        mtci,
        QP = "\n- Q-Profile method",
        BJ = "\n- Biggerstaff and Jackson method",
        J  = "\n- Jackson method",
        PL = "\n- Profile-Likelihood method")
      #
      if (mtci != "")
        details.tau.ci.i <-
          paste0(details.tau.ci.i,
                 " for confidence interval of ",
                 if (print.tau2 & print.tau)
                   paste(text.tau2, "and", text.tau)
                 else if (print.tau2)
                   text.tau2
                 else
                   text.tau)
    }
    #
    details.tau.ci <- c(details.tau.ci, details.tau.ci.i)
  }
  #
  details <- c(details, details.tau.ci)
  
  #
  #
  # (5) Method to estimate I2
  #
  #
  
  details.I2 <- NULL
  #
  if (print.I2) {
    method.I2 <- unique(meth$method.I2)
    #
    for (mi2i in method.I2) {
      details.I2.i <-
        paste0("\n- Calculation of ", text.I2,
               " based on ", switch(mi2i, Q = "Q", tau2 = text.tau2))
    }
    #
    details.I2 <- c(details.I2, details.I2.i)
  }
  #
  details <- c(details, details.I2)
  
  #
  #
  # (6) Confidence interval of common effect estimate
  #
  #
    
  if (common &&
      any(meth$type == "common" & meth$method.common.ci == "IVhet"))
    details <- c(details, "\n- Inverse variance heterogeneity model")
  
  #
  #
  # (7) Confidence interval of random effects estimate
  #
  #
  
  details.random.ci <- NULL
  #
  if (random) {
    vars <- c("method", "method.interval", "df", "three.level")
    #
    dat.rc <- unique(meth[meth$type == "random", vars])
    dat.rc.hk <-
      unique(meth[meth$type == "random", c(vars, "adhoc.hakn")])
    #
    dat.rc.hk <-
      subset(dat.rc.hk,
             dat.rc.hk$method.interval == "HK" &
               !(dat.rc.hk$method %in% c("GLMM", "LRP") | dat.rc.hk$three.level))
    #
    dat.rc.hk.tdist <-
      subset(dat.rc,
             dat.rc$method.interval == "HK" &
               (dat.rc$method %in% c("GLMM", "LRP") | dat.rc$three.level))
    #
    dat.rc.ckr <- subset(dat.rc, dat.rc$method.interval == "classic-KR")
    dat.rc.kr <- subset(dat.rc, dat.rc$method.interval == "KR")
    dat.rc.cr <-
      subset(dat.rc, dat.rc$method.interval %in% c("CR0", "CR1", "CR2"))
    #
    more.ci <- sum(1L * (nrow(dat.rc.hk) > 0) +
                     1L * (nrow(dat.rc.hk.tdist) > 0) +
                     1L * (nrow(dat.rc.ckr) > 0) +
                     1L * (nrow(dat.rc.kr) > 0) +
                     1L * (nrow(dat.rc.cr) > 0)) > 1
    #
    if (nrow(dat.rc.hk) > 0) {
      details.random.ci <-
        c(details.random.ci,
          paste0(
            "\n- Hartung-Knapp ",
            if (more.ci) "(HK) ",
            "adjustment for random effects model",
            if (print.df)
              paste0(" (df = ", cond(dat.rc.hk$df, digits = 0), ")")
          )
        )
      #
      if (any(dat.rc.hk$adhoc.hakn != "")) {
        details.random.ci <-
          c(details.random.ci,
            paste0(
              if (forest) " " else "\n  ", "(with ",
              if (any(dat.rc.hk$adhoc.hakn == ""))
                "and without ",
              "ad hoc correction)"
              )
            )
      }
    }
    #
    if (nrow(dat.rc.hk.tdist) > 0) {
      details.random.ci <-
        c(details.random.ci,
          paste0(
            "\n- Random effects confidence interval based on t-distribution",
            if (more.ci) " (T)",
            if (print.df)
              paste0(" (df = ", cond(dat.rc.hk.tdist$df, digits = 0), ")")
          )
        )
    }
    #
    if (nrow(dat.rc.ckr) > 0) {
      details.random.ci <-
        c(details.random.ci,
          paste0(
            "\n- Classic method instead of Kenward-Roger adjustment ",
            if (more.ci) "(classic-KR) ",
            "used for random effects model"
          )
        )
    }
    #
    if (nrow(dat.rc.kr) > 0) {
      details.random.ci <-
        c(details.random.ci,
          paste0(
            "\n- Kenward-Roger ",
            if (more.ci) "(KR) ",
            "adjustment for random effects model",
            if (print.df)
              paste0(" (df = ", cond(dat.rc.kr$df), ")")
          )
        )
    }
    #
    if (nrow(dat.rc.cr) > 0) {
      details.cr <- NULL
      #
      for (i in seq_len(nrow(dat.rc.cr))) {
        df.cr <- dat.rc.cr$df[i]
        #
        details.cr.i <-
          paste0(
            "\n- Cluster-robust variance estimator (",
            dat.rc.cr$method.interval[i],
            if (print.df)
              paste0(
                ", df = ",
                cond(df.cr,
                     digits = if (is.finite(df.cr) && df.cr == round(df.cr))
                       0
                     else
                       2)),
            ")"
          )
        #
        details.cr <- c(details.cr, details.cr.i)
      }
      #
      details.random.ci <- c(details.random.ci, details.cr)
    }
  }
  #
  details <- c(details, details.random.ci)
  
  #
  #
  # (8) Prediction interval
  #
  #
  
  details.pi <- NULL
  #
  if (prediction) {
    dat.pr <- unique(predict)
    #
    methods.pi <- c("V", "V-KR", "HTS", "HK", "HTS-KR", "KR", "NNF", "S")
    methods.pi <- methods.pi[methods.pi %in% dat.pr$method.interval]
    #
    more.pi <- length(methods.pi) > 1
    #
    for (mpi in methods.pi) {
      dat.pr.i <- subset(dat.pr, method.interval == mpi)
      #
      details.pi.i <- switch(
        mpi,
        "HTS-KR" =,
        "V-KR" =,
        HTS =,
        V = "\n- Prediction interval based on t-distribution",
        HK = "\n- Hartung-Knapp",
        KR = "\n- Kenward-Roger",
        NNF = "\n- Boot-strap prediction interval",
        S = "\n- Prediction interval based on standard normal distribution"
      )
      #
      if (more.pi) {
        label.pi <- switch(
          mpi,
          `V-KR` = "V",
          `HTS-KR` = "HTS",
          mpi
        )
        details.pi.i <- paste0(details.pi.i, " (", label.pi, ")")
      }
      #
      if (mpi %in% c("HK", "KR"))
        details.pi.i <- paste0(details.pi.i, " prediction interval")
      #
      if (print.df && mpi != "S") {
        details.pi.i <-
          paste0(
            details.pi.i,
            " (df = ",
            if (mpi == "KR")
              cond(dat.pr.i$df)
            else
              cond(dat.pr.i$df, digits = 0),
            ")"
          )
      }
      #
      if (mpi %in% c("V-KR", "HTS-KR")) {
        details.pi.i <-
          paste0(details.pi.i, " instead of Kenward-Roger adjustment")
      }
      #
      if (mpi == "HK" && any(dat.pr.i$adhoc.hakn.pi != "")) {
        details.pi.i <-
          paste0(
            details.pi.i,
            if (forest) " " else "\n  ",
            "(with ",
            if (any(dat.pr.i$adhoc.hakn.pi == "")) "and without ",
            "ad hoc correction)"
          )
      }
      #
      details.pi <- c(details.pi, details.pi.i)
    }
  }
  #
  details <- c(details, details.pi)
  
  #
  #
  # (9) Trim-and-fill method
  #
  #
  
  details.trimfill <- NULL
  #
  if ("trimfill" %in% overall$class) {
    type <- unique(meth$type.tf)
    type <- type[!is.na(type) & type != ""]
    #
    if (length(type) > 0)
      details.trimfill <-
      paste0(
        "\n- Trim-and-fill method to adjust for funnel plot asymmetry",
        if (length(type) == 1 || forest) " (" else "\n  (",
        paste0(type, "-estimator", collapse = ", "),
        ")"
      )
  }
  #
  details <- c(details, details.trimfill)  
  
  #
  #
  # (10) R function metamiss()
  #
  #
  
  details.miss <- NULL
  #
  if ("metamiss" %in% overall$class) {
    dat.miss <- unique(
      meth %>%
        select(method.miss, IMOR.e, IMOR.c) %>%
        filter(!is.na(method.miss))
    )
    #
    meths <- c(
      "Gamble-Hollis analysis",
      "impute no events (ICA-0)",
      "impute events (ICA-1)",
      "observed risk in control group (ICA-pc)",
      "observed risk in experimental group (ICA-pe)",
      "observed group-specific risks (ICA-p)",
      "best-case scenario (ICA-b)",
      "worst-case scenario (ICA-w)"
    )
    #
    mm <- c("GH", "0", "1", "pc", "pe", "p", "b", "w")
    #
    for (i in seq_len(nrow(dat.miss))) {
      method.miss <- dat.miss$method.miss[i]
      IMOR.e <- dat.miss$IMOR.e[i]
      IMOR.c <- dat.miss$IMOR.c[i]
      #
      if (method.miss == "IMOR")
        mmiss <-
        paste0(
          if (forest) " " else "\n  ",
          "Informative Missingness Odds Ratio",
          " (IMOR.e = ", IMOR.e,
          ", IMOR.c = ", IMOR.c, ")"
        )
      else {
        idx <- charmatch(method.miss, mm)
        mmiss <- meths[idx]
      }
      #
      details.miss <- c(
        details.miss,
        paste0("\n- Imputation method: ", mmiss)
      )
    }
  }
  #
  details <- c(details, details.miss)  
  
  #
  #
  # (11) Information on effect measures and (back-)transformations
  #
  #
  
  details.sm <- NULL
  #
  for (sm.i in unique(as.character(meth$sm[!is.na(meth$sm)]))) {
    details.sm.i <- switch(
      sm.i,
      ZCOR = "\n- Fisher's z transformation of correlations",
      COR = "\n- Untransformed correlations",
      PFT = "\n- Freeman-Tukey double arcsine transformation",
      IRFT = "\n- Freeman-Tukey transformation",
      PAS = "\n- Arcsine transformation",
      PLOGIT = "\n- Logit transformation",
      PRAW = "\n- Untransformed proportions",
      IR = "\n- Untransformed rates",
      IRS = "\n- Square root transformation",
      MRAW = "\n- Untransformed (raw) means",
      if (is_log_effect(sm.i)) "\n- Log transformation" else NULL
    )
    #
    details.sm <- c(details.sm, details.sm.i)
  }
  #
  details.func <- NULL
  #
  if (!is.null(meth$func.transf)) {
    func.transf <- unique(meth$func.transf)
    func.transf <- func.transf[!is.na(func.transf) & nzchar(func.transf)]
    #
    if (length(func.transf) > 0)
      details.func <- paste("\n- User-specified transformation:", func.transf)
  }
  #
  details.bfunc <- NULL
  #
  if (backtransf && !is.null(meth$func.backtransf)) {
    func.backtransf <- unique(meth$func.backtransf)
    func.backtransf <-
      func.backtransf[!is.na(func.backtransf) & nzchar(func.backtransf)]
    #
    if (length(func.backtransf) > 0) {
      details.bfunc <-
        paste("\n- User-specified back-transformation:", func.backtransf)
    }
  }
  #
  # Add information on method for standardised mean difference
  #
  details.smd <- NULL
  #
  if (!is.null(meth$method.smd)) {
    dat.ms <- unique(
      meth[meth$method.smd %in% c("Hedges", "Cohen", "Glass"),
           c("method.smd", "exact.smd", "sd.glass"),
           drop = FALSE])
    #
    for (i in seq_len(nrow(dat.ms))) {
      method.smd.i <- dat.ms$method.smd[i]
      exact.smd.i <- dat.ms$exact.smd[i]
      #
      if (method.smd.i == "Hedges") {
        details.smd.i <-
          paste0("\n- Hedges' g (bias corrected standardised mean difference",
                 if (isTRUE(exact.smd.i))
                   paste0(if (width >= 80 | forest) "; " else ";\n  ",
                          "using exact formulae)")
                 else
                   ")")
      }
      else if (method.smd.i == "Cohen") {
        details.smd.i <-
          paste0("\n- Cohen's d (standardised mean difference",
                 if (isTRUE(exact.smd.i))
                   "; using exact formulae)"
                 else
                   ")")
      }
      else if (method.smd.i == "Glass") {
        details.smd.i <-
          paste0("\n- Glass' delta (standardised mean difference; ",
                 if (dat.ms$sd.glass[i] == "control")
                   "based on second group)"
                 else
                   "based on first group)")
      }
      #
      details.smd <- c(details.smd, details.smd.i)
    }
  }
  #
  details <- c(details, details.sm, details.func, details.bfunc, details.smd)
  
  #
  #
  # (12) Information on continuity correction
  #
  #
  
  details.cc <- NULL
  #
  if (common_or_random &&
      any(c("metabin", "metainc", "metaprop", "metarate") %in% meth$class)) {
    vars <- c("class", "method", "sm", "incr", "method.incr", "sparse", "k.all")
    #
    if ("metabin" %in% meth$class)
      vars <- c(vars, "allstudies", "doublezeros", "MH.exact", "RR.Cochrane")
    #
    dat.cc <-
      unique(
        meth[, vars, drop = FALSE] %>%
          filter(class %in% c("metabin", "metainc", "metaprop", "metarate"),
                 !is.na(method.incr)))
    #
    for (i in seq_len(nrow(dat.cc))) {
      metabin.i <- dat.cc$class[i] == "metabin"
      sm.i <- dat.cc$sm[i]
      #
      if (!(metabin.i && sm.i == "ASD")) {
        details.cc.i <- NULL
        method.i <- dat.cc$method[i]
        method.incr.i <- dat.cc$method.incr[i]
        incr.i <- if (!is.na(dat.cc$incr[i])) dat.cc$incr[i] else "0"
        sparse.i <- isTRUE(dat.cc$sparse[i])
        k.all.i <- dat.cc$k.all[i]
        #
        details.rr <- NULL
        #
        txtCC.ind.i <-
          method.i %in% c("GLMM", "LRP") ||
          (metabin.i && method.i == "MH" &&
             isTRUE(dat.cc$MH.exact[i]))
        #
        if (metabin.i &&
            isTRUE(dat.cc$RR.Cochrane[i]) &&
            (method.incr.i == "all" || (sparse.i && incr.i != "0"))) {
          details.rr <-
            if (width >= 70 || forest)
              " (applied twice to sample sizes, like RevMan 5)"
          else
            "\n  (applied twice to sample sizes, like RevMan 5)"
        }
        #
        # Describe corrections affecting the pooled estimates only.
        #
        if (!txtCC.ind.i) {
          if (method.incr.i == "user") {
            if (incr.i != "0")
              details.cc.i <- "\n- User-defined continuity correction"
          }
          else if (method.incr.i == "all") {
            if (incr.i == "TACC") {
              details.cc.i <-
                if (k.all.i > 1)
                  "\n- Treatment arm continuity correction in all studies"
              else
                "\n- Treatment arm continuity correction"
            }
            else if (incr.i != "0") {
              details.cc.i <-
                paste0("\n- Continuity correction of ",
                       incr.i,
                       if (k.all.i > 1) " in all studies",
                       details.rr)
            }
          }
          else if (sparse.i) {
            if (incr.i == "TACC") {
              details.cc.i <-
                paste0(
                  "\n- Treatment arm continuity correction in ",
                  if (k.all.i > 1) "studies " else "study ",
                  "with",
                  if (width > 70 || forest) " " else "\n  ",
                  "zero cell frequencies",
                  details.rr)
            }
            else if (incr.i != "0") {
              details.cc.i <-
                paste0(
                  "\n- Continuity correction of ",
                  incr.i,
                  if (k.all.i > 1) " in studies with",
                  if (k.all.i > 1)
                    if (width > 70 || forest) " " else "\n  ",
                  if (k.all.i > 1) "zero cell frequencies",
                  details.rr)
            }
          }
        }
        #
        if (metabin.i && method.incr.i != "user" &&
            isTRUE(dat.cc$allstudies[i]) &&
            isTRUE(dat.cc$doublezeros[i])) {
          details.cc.i <-
            c(
              details.cc.i,
              if (k.all.i > 1)
                "\n- Studies with double zeros included in meta-analysis"
              else
                "\n- Study with double zeros considered"
            )
        }
        #
        details.cc <- c(details.cc, details.cc.i)
      }
    }
  }
  #
  details <- c(details, details.cc)  
  
  #
  #
  # (13) Information on number of events
  #
  #
  
  if (!is.null(pscale) && pscale != 1 &&
      (any(is_prop(meth$sm)) || any(meth$sm == "RD"))) {
    details <-
      c(details,
        paste0(
          "\n- Events per ",
          format(pscale, scientific = FALSE, big.mark = big.mark),
          " observations"))
  }
  #
  if (!is.null(irscale) && irscale != 1 &&
      (any(is_rate(meth$sm)) || any(meth$sm == "IRD"))) {
    details <-
      c(details,
        paste0(
          "\n- Events per ",
          format(irscale, scientific = FALSE, big.mark = big.mark),
          " ", irunit))
  }
  
  #
  #
  # (14) Information on null hypothesis
  #
  #
  
  if ((common || random) && print.pval) {
    details.null <- NULL
    dat.null <- NULL
    #
    for (i in seq_len(nrow(meth))) {
      null.effect <- meth$null.effect[i]
      #
      if (is.null(null.effect) || is.na(null.effect))
        next
      #
      sm.i <- meth$sm[i]
      null.effect.original <-
        meth$class[i] %in% c("metacor", "metamean", "metaprop", "metarate")
      #
      # Use [[i]] for scalar columns and transformation-argument list-columns
      #
      func.transf.i <- meth$func.transf[[i]]
      func.backtransf.i <- meth$func.backtransf[[i]]
      args.transf.i <- meth$args.transf[[i]]
      args.backtransf.i <- meth$args.backtransf[[i]]
      #
      # bind_rows() can introduce NA for absent transformation names.
      #
      if (is.character(func.transf.i) && anyNA(func.transf.i))
        func.transf.i <- NULL
      #
      if (is.character(func.backtransf.i) && anyNA(func.backtransf.i))
        func.backtransf.i <- NULL
      #
      n.i <- if (sm.i == "PFT") meth$n.harmonic.mean[[i]] else NULL
      time.i <- if (sm.i == "IRFT") meth$t.harmonic.mean[[i]] else NULL
      #
      if (!backtransf && null.effect.original)
        null.effect <-
        transf(
          null.effect, sm.i,
          n = n.i,
          time = time.i,
          func = func.transf.i,
          args = args.transf.i
        )
      else if (backtransf && !null.effect.original)
        null.effect <-
        backtransf(
          null.effect, sm.i,
          n = n.i,
          time = time.i,
          func = func.backtransf.i,
          args = args.backtransf.i
        )
      #
      if (is.na(null.effect))
        next
      #
      if (!is.null(pscale) && pscale != 1) {
        prefix <- "\n- Null hypothesis: "
        null.effect <- null.effect * pscale
        suffix <- paste0(
          " events per ",
          format(pscale, scientific = FALSE, big.mark = big.mark),
          " observations")
      }
      else if (!is.null(irscale) && irscale != 1) {
        prefix <- "\n- Null hypothesis: "
        null.effect <- null.effect * irscale
        suffix <- paste0(
          " events per ",
          format(irscale, scientific = FALSE, big.mark = big.mark),
          " ", irunit)
      }
      else {
        sm.lab <- smlab(
          sm.i,
          backtransf = backtransf,
          pscale = 1,
          irscale = 1,
          forest = forest,
          func.backtransf = func.backtransf.i
        )
        #
        if (sm.lab == "")
          sm.lab <- "effect"
        #
        prefix <- paste0("\n- Null hypothesis: ", sm.lab, " = ")
        suffix <- ""
      }
      #
      dat.null <- bind_rows(
        dat.null,
        data.frame(prefix = prefix, value = null.effect, suffix = suffix)
      )
    }
    #
    if (!is.null(dat.null)) {
      labels.null <- unique(dat.null[, c("prefix", "suffix"), drop = FALSE])
      #
      for (i in seq_len(nrow(labels.null))) {
        sel.null <- dat.null$prefix == labels.null$prefix[i] &
          dat.null$suffix == labels.null$suffix[i]
        details.null <- c(
          details.null,
          paste0(labels.null$prefix[i],
                 cond(dat.null$value[sel.null], only.finite = FALSE,
                      digits = digits, big.mark = big.mark),
                 labels.null$suffix[i])
        )
      }
    }
    #
    details <- c(details, details.null)
  }
  
  #
  #
  # (15) Return details
  #
  #
  
  details <- unique(details)
  details <- details[!is.na(details) & nzchar(details)]
  
  if (length(details) == 0)
    return(invisible(""))
  
  details <- paste0(
    "\nDetails of meta-analysis methods:",
    paste0(details, collapse = "")
  )
  
  if (!forest)
    cat(details, "\n", sep = "")
  
  invisible(details)
}
