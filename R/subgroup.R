subgroup <- function(x, tau.preset = NULL, subgroup.rma, ...) {
  
  
  subgroup <- x$subgroup
  ##
  levs <- bylevs(subgroup)
  n.levs <- length(levs)
  methci <- paste(x$method.random.ci,
                  toupper(substring(x$adhoc.hakn.ci, 1, 2)),
                  sep = "-")
  methci <- gsub("-$", "", methci)
  ##
  methpi <- paste(x$method.predict,
                  toupper(substring(x$adhoc.hakn.pi, 1, 2)),
                  sep = "-")
  methpi <- gsub("-$", "", methpi)
  ##
  if (!(length(subgroup) > 0)) {
    warning("Argument 'subgroup' is missing.")
    return(NULL)
  }
  ##
  bin  <- inherits(x, "metabin")
  cont <- inherits(x, "metacont")
  cor  <- inherits(x, "metacor")
  gen  <- inherits(x, "metagen")
  mean <- inherits(x, "metamean")
  inc  <- inherits(x, "metainc")
  prop <- inherits(x, "metaprop")
  rate <- inherits(x, "metarate")
  ##
  bin.cont.gen <- bin | cont | gen
  bin.inc <- bin | inc
  cor.prop.mean <- cor | prop | mean
  ##
  three.level <- !is.null(x$three.level) && x$three.level
  ##
  sumNA <- function(x)
    if (all(is.na(x)))
      NA
    else
      sum(x, na.rm = TRUE)
  
  
  ##
  ##
  ## (1) Subgroup analysis without common tau2
  ##
  ##
  res.i <- vector(mode = "list")
  ##
  res.w <- matrix(NA, ncol = 63, nrow = n.levs)
  add.w <- matrix("", ncol =  2, nrow = n.levs)
  j <- 0
  ##
  for (i in levs) {
    j <- j + 1
    sel <- subgroup == i
    ##
    if (all(is.na(x$studlab[sel])))
      stop("No data available for subgroup = ", i)
    ##
    ##
    if (bin)
      meta1 <- metabin(x$event.e[sel], x$n.e[sel],
                       x$event.c[sel], x$n.c[sel],
                       studlab = x$studlab[sel],
                       exclude = x$exclude[sel],
                       cluster =
                         if (!is.null(x$cluster)) x$cluster[sel] else NULL,
                       ##
                       method = x$method,
                       sm = x$sm,
                       incr = if (length(x$incr) == 1) x$incr else x$incr[sel],
                       method.incr = x$method.incr,
                       allstudies = x$allstudies,
                       MH.exact = x$MH.exact,
                       RR.Cochrane = x$RR.Cochrane,
                       Q.Cochrane = x$Q.Cochrane,
                       model.glmm = x$model.glmm,
                       ##
                       level.ma = x$level.ma,
                       method.random.ci = x$method.random.ci,
                       adhoc.hakn.ci = x$adhoc.hakn.ci,
                       ##
                       level.predict = x$level.predict,
                       method.predict = x$method.predict,
                       adhoc.hakn.pi = x$adhoc.hakn.pi,
                       ##
                       method.tau = x$method.tau,
                       method.tau.ci = x$method.tau.ci,
                       tau.preset = tau.preset,
                       TE.tau = x$TE.tau,
                       ##
                       keepdata = FALSE,
                       warn = x$warn,
                       control = x$control)
    ##
    else if (cont)
      meta1 <- metacont(x$n.e[sel], x$mean.e[sel],
                        x$sd.e[sel],
                        x$n.c[sel], x$mean.c[sel],
                        x$sd.c[sel],
                        studlab = x$studlab[sel],
                        exclude = x$exclude[sel],
                        cluster =
                          if (!is.null(x$cluster)) x$cluster[sel] else NULL,
                        ##
                        sm = x$sm,
                        ##
                        pooledvar = x$pooledvar,
                        method.smd = x$method.smd,
                        sd.glass = x$sd.glass,
                        exact.smd = x$exact.smd,
                        ##
                        level.ma = x$level.ma,
                        method.random.ci = x$method.random.ci,
                        adhoc.hakn.ci = x$adhoc.hakn.ci,
                        ##
                        level.predict = x$level.predict,
                        method.predict = x$method.predict,
                        adhoc.hakn.pi = x$adhoc.hakn.pi,
                        ##
                        method.tau = x$method.tau,
                        method.tau.ci = x$method.tau.ci,
                        tau.preset = tau.preset,
                        TE.tau = x$TE.tau,
                        ##
                        keepdata = FALSE,
                        warn = x$warn,
                        control = x$control)
    ##
    else if (cor)
      meta1 <- metacor(x$cor[sel], x$n[sel],
                       studlab = x$studlab[sel],
                       exclude = x$exclude[sel],
                       cluster =
                         if (!is.null(x$cluster)) x$cluster[sel] else NULL,
                       ##
                       sm = x$sm,
                       ##
                       level.ma = x$level.ma,
                       method.random.ci = x$method.random.ci,
                       adhoc.hakn.ci = x$adhoc.hakn.ci,
                       ##
                       level.predict = x$level.predict,
                       method.predict = x$method.predict,
                       adhoc.hakn.pi = x$adhoc.hakn.pi,
                       ##
                       method.tau = x$method.tau,
                       method.tau.ci = x$method.tau.ci,
                       tau.preset = tau.preset,
                       TE.tau = x$TE.tau,
                       ##
                       null.effect = x$null.effect,
                       ##
                       keepdata = FALSE,
                       control = x$control)
    ##
    else if (gen)
      meta1 <- metagen(x$TE[sel], x$seTE[sel],
                       studlab = x$studlab[sel],
                       exclude = x$exclude[sel],
                       cluster =
                         if (!is.null(x$cluster)) x$cluster[sel] else NULL,
                       ##
                       sm = x$sm,
                       ##
                       level.ma = x$level.ma,
                       method.random.ci = x$method.random.ci,
                       adhoc.hakn.ci = x$adhoc.hakn.ci,
                       ##
                       level.predict = x$level.predict,
                       method.predict = x$method.predict,
                       adhoc.hakn.pi = x$adhoc.hakn.pi,
                       ##
                       method.tau = x$method.tau,
                       method.tau.ci = x$method.tau.ci,
                       tau.preset = tau.preset,
                       TE.tau = x$TE.tau,
                       ##
                       null.effect = x$null.effect,
                       n.e = x$n.e[sel], n.c = x$n.c[sel],
                       ##
                       keepdata = FALSE,
                       warn = x$warn,
                       control = x$control)
    ##
    else if (inc)
      meta1 <- metainc(x$event.e[sel], x$time.e[sel],
                       x$event.c[sel], x$time.c[sel],
                       studlab = x$studlab[sel],
                       exclude = x$exclude[sel],
                       cluster =
                         if (!is.null(x$cluster)) x$cluster[sel] else NULL,
                       ##
                       method = x$method,
                       sm = x$sm,
                       incr = if (length(x$incr) == 1) x$incr else x$incr[sel],
                       method.incr = x$method.incr,
                       model.glmm = x$model.glmm,
                       ##
                       level.ma = x$level.ma,
                       method.random.ci = x$method.random.ci,
                       adhoc.hakn.ci = x$adhoc.hakn.ci,
                       ##
                       level.predict = x$level.predict,
                       method.predict = x$method.predict,
                       adhoc.hakn.pi = x$adhoc.hakn.pi,
                       ##
                       method.tau = x$method.tau,
                       method.tau.ci = x$method.tau.ci,
                       tau.preset = tau.preset,
                       TE.tau = x$TE.tau,
                       ##
                       n.e = x$n.e[sel], n.c = x$n.c[sel],
                       ##
                       keepdata = FALSE,
                       warn = x$warn,
                       control = x$control)
    ##
    else if (mean)
      meta1 <- metamean(x$n[sel], x$mean[sel], x$sd[sel],
                        studlab = x$studlab[sel],
                        exclude = x$exclude[sel],
                        cluster =
                          if (!is.null(x$cluster)) x$cluster[sel] else NULL,
                       ##
                       sm = x$sm,
                       ##
                       level.ma = x$level.ma,
                       method.random.ci = x$method.random.ci,
                       adhoc.hakn.ci = x$adhoc.hakn.ci,
                       ##
                       level.predict = x$level.predict,
                       method.predict = x$method.predict,
                       adhoc.hakn.pi = x$adhoc.hakn.pi,
                       ##
                       method.tau = x$method.tau,
                       method.tau.ci = x$method.tau.ci,
                       tau.preset = tau.preset,
                       TE.tau = x$TE.tau,
                       ##
                       null.effect = x$null.effect,
                       ##
                       keepdata = FALSE,
                       warn = x$warn,
                       control = x$control)
    ##
    else if (prop)
      meta1 <- metaprop(x$event[sel], x$n[sel],
                        studlab = x$studlab[sel],
                        exclude = x$exclude[sel],
                        cluster =
                          if (!is.null(x$cluster)) x$cluster[sel] else NULL,
                        ##
                        method = x$method,
                        sm = x$sm,
                        incr = if (length(x$incr) == 1) x$incr else x$incr[sel],
                        method.incr = x$method.incr,
                        ##
                        level.ma = x$level.ma,
                        method.random.ci = x$method.random.ci,
                        adhoc.hakn.ci = x$adhoc.hakn.ci,
                        ##
                        level.predict = x$level.predict,
                        method.predict = x$method.predict,
                        adhoc.hakn.pi = x$adhoc.hakn.pi,
                        ##
                        method.tau = x$method.tau,
                        method.tau.ci = x$method.tau.ci,
                        tau.preset = tau.preset,
                        TE.tau = x$TE.tau,
                        ##
                        null.effect = x$null.effect,
                        ##
                        keepdata = FALSE,
                        warn = x$warn,
                        control = x$control)
    ##
    else if (rate)
      meta1 <- metarate(x$event[sel], x$time[sel],
                        studlab = x$studlab[sel],
                        exclude = x$exclude[sel],
                        cluster =
                          if (!is.null(x$cluster)) x$cluster[sel] else NULL,
                        ##
                        method = x$method,
                        sm = x$sm,
                        incr = if (length(x$incr) == 1) x$incr else x$incr[sel],
                        method.incr = x$method.incr,
                        ##
                        level.ma = x$level.ma,
                        method.random.ci = x$method.random.ci,
                        adhoc.hakn.ci = x$adhoc.hakn.ci,
                        ##
                        level.predict = x$level.predict,
                        method.predict = x$method.predict,
                        adhoc.hakn.pi = x$adhoc.hakn.pi,
                        ##
                        method.tau = x$method.tau,
                        method.tau.ci = x$method.tau.ci,
                        tau.preset = tau.preset,
                        TE.tau = x$TE.tau,
                        ##
                        null.effect = x$null.effect,
                        ##
                        keepdata = FALSE,
                        warn = x$warn,
                        control = x$control)
    ##
    else
      stop("No meta-analysis object used for subgroup analysis.")
    ##
    if (length(meta1$TE.random) == 1 &&
        length(meta1$TE.random) != length(meta1$seTE.random))
      meta1$TE.random <- rep_len(meta1$TE.random, length(meta1$seTE.random))
    ##
    res.i[[j]] <- list(k = meta1$k,
                       k.study = meta1$k.study,
                       k.all = meta1$k.all,
                       k.TE = meta1$k.TE,
                       ##
                       TE.common = meta1$TE.common,
                       seTE.common = meta1$seTE.common,
                       statistic.common = meta1$statistic.common,
                       pval.common = meta1$pval.common,
                       lower.common = meta1$lower.common,
                       upper.common = meta1$upper.common,
                       w.common = sum(x$w.common[sel]),
                       ##
                       TE.random = meta1$TE.random,
                       seTE.random = meta1$seTE.random,
                       statistic.random = meta1$statistic.random,
                       pval.random = meta1$pval.random,
                       df.random = meta1$df.random,
                       lower.random = meta1$lower.random,
                       upper.random = meta1$upper.random,
                       w.random = sum(x$w.random[sel]),
                       ##
                       seTE.classic = meta1$seTE.classic,
                       ##
                       df.hakn.ci = meta1$df.hakn.ci,
                       seTE.hakn.ci = meta1$seTE.hakn.ci,
                       seTE.hakn.adhoc.ci = meta1$seTE.hakn.adhoc.ci,
                       ##
                       df.kero = meta1$df.kero,
                       seTE.kero = meta1$seTE.kero,
                       ##
                       seTE.predict = meta1$seTE.predict,
                       df.predict = meta1$df.predict,
                       lower.predict = meta1$lower.predict,
                       upper.predict = meta1$upper.predict,
                       ##
                       seTE.hakn.pi = meta1$seTE.hakn.pi,
                       seTE.hakn.adhoc.pi = meta1$seTE.hakn.adhoc.pi,
                       ##
                       Q = meta1$Q,
                       ##
                       tau2 = sum(meta1$tau2),
                       tau = sum(meta1$tau),
                       ##
                       H = meta1$H,
                       lower.H = meta1$lower.H,
                       upper.H = meta1$upper.H,
                       ##
                       I2 = meta1$I2,
                       lower.I2 = meta1$lower.I2,
                       upper.I2 = meta1$upper.I2,
                       ##
                       Rb = meta1$Rb,
                       lower.Rb = meta1$lower.Rb,
                       upper.Rb = meta1$upper.Rb,
                       ##
                       event = if (prop) sumNA(meta1$event) else NA,
                       n = if (cor.prop.mean) sumNA(meta1$n) else NA,
                       ##
                       event.e = if (bin.inc) sumNA(meta1$event.e) else NA,
                       n.e = if (bin.cont.gen) sumNA(meta1$n.e) else NA,
                       event.c = if (bin.inc) sumNA(meta1$event.c) else NA,
                       n.c = if (bin.cont.gen) sumNA(meta1$n.c) else NA,
                       ##
                       time.e = if (inc) sumNA(meta1$time.e) else NA,
                       time.c = if (inc) sumNA(meta1$time.c) else NA,
                       ##
                       n.harmonic.mean = 1 / mean(1 / x$n[sel]),
                       t.harmonic.mean = 1 / mean(1 / x$time[sel]))
  }
  ##
  k.w <- extrVec(res.i, "k", levs)
  k.study.w <- extrVec(res.i, "k.study", levs)
  k.all.w <- extrVec(res.i, "k.all", levs)
  k.TE.w <- extrVec(res.i, "k.TE", levs)
  ##
  TE.common.w <- extrVec(res.i, "TE.common", levs)
  seTE.common.w <- extrVec(res.i, "seTE.common", levs)
  statistic.common.w <- extrVec(res.i, "statistic.common", levs)
  pval.common.w <- extrVec(res.i, "pval.common", levs)
  lower.common.w <- extrVec(res.i, "lower.common", levs)
  upper.common.w <- extrVec(res.i, "upper.common", levs)
  w.common.w <- extrVec(res.i, "w.common", levs)
  ##
  TE.random.w <- extrMat(res.i, "TE.random", levs, methci)    
  seTE.random.w <- extrMat(res.i, "seTE.random", levs, methci)
  if (is.matrix(TE.random.w)) {
    TE.random.w <- as.vector(TE.random.w[, 1])
    names(TE.random.w) <- rownames(seTE.random.w)
  }
  statistic.random.w <- extrMat(res.i, "statistic.random", levs, methci)
  pval.random.w <- extrMat(res.i, "pval.random", levs, methci)
  df.random.w <- extrMat(res.i, "df.random", levs, methci)
  lower.random.w <- extrMat(res.i, "lower.random", levs, methci)
  upper.random.w <- extrMat(res.i, "upper.random", levs, methci)
  w.random.w <- extrVec(res.i, "w.random", levs)
  ##
  seTE.classic.w <- extrVec(res.i, "seTE.classic", levs)
  ##
  df.hakn.ci.w <- extrMat(res.i, "df.hakn.ci", levs, methci)
  seTE.hakn.ci.w <- extrVec(res.i, "seTE.hakn.ci", levs)
  seTE.hakn.adhoc.ci.w <- extrMat(res.i, "seTE.hakn.adhoc.ci", levs, methci)
  ##
  df.kero.w <- extrVec(res.i, "df.kero", levs)
  seTE.kero.w <- extrVec(res.i, "seTE.kero", levs)
  ##
  seTE.predict.w <- extrMat(res.i, "seTE.predict", levs, methpi)
  df.predict.w <- extrMat(res.i, "df.predict", levs, methpi)
  lower.predict.w <- extrMat(res.i, "lower.predict", levs, methpi)
  upper.predict.w <- extrMat(res.i, "upper.predict", levs, methpi)
  ##
  seTE.hakn.pi.w <- extrVec(res.i, "seTE.hakn.pi", levs)
  seTE.hakn.adhoc.pi.w <- extrMat(res.i, "seTE.hakn.adhoc.pi", levs, methpi)
  ##
  Q.w <- extrVec(res.i, "Q", levs)
  ##
  tau2.w <- extrVec(res.i, "tau2", levs)
  tau.w <- extrVec(res.i, "tau", levs)
  ##
  H.w <- extrVec(res.i, "H", levs)
  lower.H.w <- extrVec(res.i, "lower.H", levs)
  upper.H.w <- extrVec(res.i, "upper.H", levs)
  ##
  I2.w <- extrVec(res.i, "I2", levs)
  lower.I2.w <- extrVec(res.i, "lower.I2", levs)
  upper.I2.w <- extrVec(res.i, "upper.I2", levs)
  ##
  Rb.w <- extrVec(res.i, "Rb", levs)
  lower.Rb.w <- extrVec(res.i, "lower.Rb", levs)
  upper.Rb.w <- extrVec(res.i, "upper.Rb", levs)
  ##
  event.w <- extrVec(res.i, "event", levs)
  n.w <- extrVec(res.i, "n", levs)
  ##
  event.e.w <- extrVec(res.i, "event.e", levs)
  n.e.w <- extrVec(res.i, "n.e", levs)
  event.c.w <- extrVec(res.i, "event.c", levs)
  n.c.w <- extrVec(res.i, "n.c", levs)
  ##
  time.e.w <- extrVec(res.i, "time.e", levs)
  time.c.w <- extrVec(res.i, "time.c", levs)
  ##
  n.harmonic.mean.w <- extrVec(res.i, "n.harmonic.mean", levs)
  t.harmonic.mean.w <- extrVec(res.i, "t.harmonic.mean", levs)
  ##
  ## Tests for subgroup differences
  ##
  Q.w.common <- sum(Q.w, na.rm = TRUE)
  df.Q.w <- sum((k.w - 1)[!is.na(Q.w)])
  pval.Q.w.common <- pvalQ(Q.w.common, df.Q.w)
  ##
  Q.b.common <- metagen(TE.common.w, seTE.common.w, method.tau = "DL")$Q
  ##
  if (is.matrix(seTE.random.w)) {
    n.random <- ncol(seTE.random.w)
    Q.b.random <- rep(NA, n.random)
    names(Q.b.random) <- colnames(seTE.random.w)
    for (i in seq_len(n.random))
      Q.b.random[i] <-
        metagen(TE.random.w, seTE.random.w[, i], method.tau = "DL")$Q
  }
  else {
    n.random <- 1
    Q.b.random <-
      metagen(TE.random.w, seTE.random.w, method.tau = "DL")$Q
  }
  ##
  df.Q.b <- if (x$k == 0) 0 else x$k - 1 - sum((k.w - 1)[!is.na(Q.w)])
  df.Q.b.common <- df.Q.b
  df.Q.b.random <- rep(df.Q.b, n.random)
  ##
  pval.Q.b.common <- pvalQ(Q.b.common, df.Q.b.common)
  pval.Q.b.random <- pvalQ(Q.b.random, df.Q.b.random)
  
  
  ##
  ##
  ## (2) Subgroup analysis with common tau-squared (three-level model)
  ##
  ##
  if (three.level && !missing(subgroup.rma)) {
    mod <- as.call(~ subgroup.rma - 1)
    mod.Q <- as.call(~ subgroup.rma)
    ##
    cluster <- x$cluster
    runID <- seq_along(cluster)
    ##
    mv.random <-
      runNN(rma.mv,
            list(yi = x$TE, V = x$seTE^2,
                 mods = mod,
                 random = as.call(~ 1 | cluster / runID),
                 method = x$method.tau,
                 test = ifelse(x$method.random.ci == "HK", "t", "z"),
                 level = 100 * x$level.ma,
                 data = data.frame(subgroup.rma, cluster, runID)))
    ##
    mv.random.Q <-
      suppressWarnings(
        runNN(rma.mv,
              list(yi = x$TE, V = x$seTE^2,
                   mods = mod.Q,
                   random = as.call(~ 1 | cluster / runID),
                   method = x$method.tau,
                   test = ifelse(x$method.random.ci == "HK", "t", "z"),
                   level = 100 * x$level.ma,
                   data = data.frame(subgroup.rma, cluster, runID))))
    ##
    if (length(TE.random.w) != length(as.numeric(mv.random$b))) {
      TE.random.w[!is.na(TE.random.w)] <- as.numeric(mv.random$b)
      seTE.random.w[!is.na(seTE.random.w)] <- as.numeric(mv.random$se)
    }
    else {
      TE.random.w   <- as.numeric(mv.random$b)
      seTE.random.w <- as.numeric(mv.random$se)
    }
    ##
    tau2.w <- rep_len(sum(mv.random$sigma2), n.levs)
    tau.w <- sqrt(tau2.w)
    ##
    Rb.w     <- rep_len(NA, n.levs)
    Rb.w.low <- rep_len(NA, n.levs)
    Rb.w.upp <- rep_len(NA, n.levs)
    ##
    ci.common.w  <- ci(TE.common.w, seTE.common.w, x$level.ma)
    ##
    if (!is.null(x$method.random.ci == "HK") && x$method.random.ci == "HK")
      ci.random.w <- ci(TE.random.w, seTE.random.w, x$level.ma, df = k.w - 1)
    else
      ci.random.w <- ci(TE.random.w, seTE.random.w, x$level.ma)
    ##
    lower.common.w <- ci.common.w$lower
    upper.common.w <- ci.common.w$upper
    statistic.common.w <- ci.common.w$statistic
    pval.common.w <- ci.common.w$p
    ##
    lower.random.w <- ci.random.w$lower
    upper.random.w <- ci.random.w$upper
    statistic.random.w <- ci.random.w$statistic
    pval.random.w <- ci.random.w$p
    ##
    ## Tests for subgroup differences
    ##
    Q.w.common <- NA
    df.Q.w <- mv.random.Q$k.eff - mv.random.Q$p.eff
    pval.Q.w.common <- NA
    ##
    Q.b.common <- NA
    Q.b.random <- mv.random.Q$QM
    ##
    df.Q.b <- mv.random.Q$QMdf
    df.Q.b <- df.Q.b[!is.na(df.Q.b)]
    df.Q.b.common <- NA
    df.Q.b.random <- df.Q.b
    ##
    pval.Q.b.common <- NA
    pval.Q.b.random <- mv.random.Q$QMp
    ##
    ## Prediction interval
    ##
    seTE.predict.w <- sqrt(seTE.random.w^2 + tau2.w)
    df.predict.w <- k.w - 2
    ci.p.w <- ci(TE.random.w, seTE.predict.w, x$level.predict, df.predict.w)
    ##
    lower.predict.w <- ci.p.w$lower
    upper.predict.w <- ci.p.w$upper
    ##
    lower.predict.w[k.w < 3] <- NA
    upper.predict.w[k.w < 3] <- NA
    ##
    ## Degrees of freedom of Hartung-Knapp method
    ##
    df.random.w <- df.hakn.ci.w <- k.w - 1
    if (!x$method.random.ci == "HK") {
      df.random.w[!is.na(df.random.w)] <- NA
      df.hakn.ci.w[!is.na(df.hakn.ci.w)] <- NA
    }
  }
  
  
  ##
  ##
  ## (3) Subgroup analysis with common tau-squared (GLMM)
  ##
  ##
  if (x$method == "GLMM" & !missing(subgroup.rma)) {
    mod <- as.call(~ subgroup.rma - 1)
    mod.Q <- as.call(~ subgroup.rma)
    ##
    if (prop) {
      glmm.common <-
        runNN(rma.glmm,
              list(xi = x$event, ni = x$n,
                   mods = mod,
                   method = "FE",
                   test = ifelse(x$method.random.ci == "HK", "t", "z"),
                   level = 100 * x$level.ma,
                   measure = "PLO",
                   data = data.frame(subgroup.rma),
                   ...))
      ##
      glmm.random <-
        runNN(rma.glmm,
              list(xi = x$event, ni = x$n,
                   mods = mod,
                   method = x$method.tau,
                   test = ifelse(x$method.random.ci == "HK", "t", "z"),
                   level = 100 * x$level.ma,
                   measure = "PLO",
                   data = data.frame(subgroup.rma),
                   ...))
      ##
      ## Test for subgroup differences
      ##
      glmm.common.Q <-
        suppressWarnings(
          runNN(rma.glmm,
                list(xi = x$event, ni = x$n,
                     mods = mod.Q,
                     method = "FE",
                     test = ifelse(x$method.random.ci == "HK", "t", "z"),
                     level = 100 * x$level.ma,
                     measure = "PLO",
                     data = data.frame(subgroup.rma),
                     ...)))
      ##
      glmm.random.Q <-
        suppressWarnings(
          runNN(rma.glmm,
                list(xi = x$event, ni = x$n,
                     mods = mod.Q,
                     method = x$method.tau,
                     test = ifelse(x$method.random.ci == "HK", "t", "z"),
                     level = 100 * x$level.ma,
                     measure = "PLO",
                     data = data.frame(subgroup.rma),
                     ...)))
    }
    ##
    else if (rate) {
      glmm.common <-
        runNN(rma.glmm,
              list(xi = x$event, ti = x$time,
                   mods = mod,
                   method = "FE",
                   test = ifelse(x$method.random.ci == "HK", "t", "z"),
                   level = 100 * x$level.ma,
                   measure = "IRLN", control = x$control,
                   data = data.frame(subgroup.rma),
                   ...))
      ##
      glmm.random <-
        runNN(rma.glmm,
              list(xi = x$event, ti = x$time,
                   mods = mod,
                   method = x$method.tau,
                   test = ifelse(x$method.random.ci == "HK", "t", "z"),
                   level = 100 * x$level.ma,
                   measure = "IRLN", control = x$control,
                   data = data.frame(subgroup.rma),
                   ...))
      ##
      ## Test for subgroup differences
      ##
      glmm.common.Q <-
        suppressWarnings(
          runNN(rma.glmm,
                list(xi = x$event, ti = x$time,
                     mods = mod.Q,
                     method = "FE",
                     test = ifelse(x$method.random.ci == "HK", "t", "z"),
                     level = 100 * x$level.ma,
                     measure = "IRLN", control = x$control,
                     data = data.frame(subgroup.rma),
                     ...)))
      ##
      glmm.random.Q <-
        suppressWarnings(
          runNN(rma.glmm,
                list(xi = x$event, ti = x$time,
                     mods = mod.Q,
                     method = x$method.tau,
                     test = ifelse(x$method.random.ci == "HK", "t", "z"),
                     level = 100 * x$level.ma,
                     measure = "IRLN", control = x$control,
                     data = data.frame(subgroup.rma),
                     ...)))
    }
    ##
    else if (inc) {
      glmm.common <-
        runNN(rma.glmm,
              list(x1i = x$event.e, t1i = x$time.e,
                   x2i = x$event.c, t2i = x$time.c,
                   mods = mod,
                   method = "FE",
                   test = ifelse(x$method.random.ci == "HK", "t", "z"),
                   model = x$model.glmm,
                   level = 100 * x$level.ma,
                   measure = "IRR", control = x$control,
                   data = data.frame(subgroup.rma),
                   ...))
      ##
      glmm.random <-
        runNN(rma.glmm,
              list(x1i = x$event.e, t1i = x$time.e,
                   x2i = x$event.c, t2i = x$time.c,
                   mods = mod,
                   method = x$method.tau,
                   model = x$model.glmm,
                   test = ifelse(x$method.random.ci == "HK", "t", "z"),
                   level = 100 * x$level.ma,
                   measure = "IRR", control = x$control,
                   data = data.frame(subgroup.rma),
                   ...))
      ##
      ## Test for subgroup differences
      ##
      glmm.common.Q <-
        suppressWarnings(
          runNN(rma.glmm,
                list(x1i = x$event.e, t1i = x$time.e,
                     x2i = x$event.c, t2i = x$time.c,
                     mods = mod.Q,
                     method = "FE",
                     test = ifelse(x$method.random.ci == "HK", "t", "z"),
                     model = x$model.glmm,
                     level = 100 * x$level.ma,
                     measure = "IRR", control = x$control,
                     data = data.frame(subgroup.rma),
                     ...)))
      ##
      glmm.random.Q <-
        suppressWarnings(
          runNN(rma.glmm,
                list(x1i = x$event.e, t1i = x$time.e,
                     x2i = x$event.c, t2i = x$time.c,
                     mods = mod.Q,
                     method = x$method.tau,
                     model = x$model.glmm,
                     test = ifelse(x$method.random.ci == "HK", "t", "z"),
                     level = 100 * x$level.ma,
                     measure = "IRR", control = x$control,
                     data = data.frame(subgroup.rma),
                     ...)))
    }
    ##
    else if (bin) {
      glmm.common <-
        runNN(rma.glmm,
              list(ai = x$event.e, n1i = x$n.e,
                   ci = x$event.c, n2i = x$n.c,
                   mods = mod,
                   method = "FE",
                   test = ifelse(x$method.random.ci == "HK", "t", "z"),
                   model = x$model.glmm,
                   level = 100 * x$level.ma,
                   measure = "OR", control = x$control,
                   data = data.frame(subgroup.rma),
                   ...))
      ##
      glmm.random <-
        runNN(rma.glmm,
              list(ai = x$event.e, n1i = x$n.e,
                   ci = x$event.c, n2i = x$n.c,
                   mods = mod,
                   method = x$method.tau,
                   model = x$model.glmm,
                   test = ifelse(x$method.random.ci == "HK", "t", "z"),
                   level = 100 * x$level.ma,
                   measure = "OR", control = x$control,
                   data = data.frame(subgroup.rma),
                   ...))
      ##
      ## Test for subgroup differences
      ##
      glmm.common.Q <-
        suppressWarnings(
          runNN(rma.glmm,
                list(ai = x$event.e, n1i = x$n.e,
                     ci = x$event.c, n2i = x$n.c,
                     mods = mod.Q,
                     method = "FE",
                     test = ifelse(x$method.random.ci == "HK", "t", "z"),
                     model = x$model.glmm,
                     level = 100 * x$level.ma,
                     measure = "OR", control = x$control,
                     data = data.frame(subgroup.rma),
                     ...)))
      ##
      glmm.random.Q <-
        suppressWarnings(
          runNN(rma.glmm,
                list(ai = x$event.e, n1i = x$n.e,
                     ci = x$event.c, n2i = x$n.c,
                     mods = mod.Q,
                     method = x$method.tau,
                     model = x$model.glmm,
                     test = ifelse(x$method.random.ci == "HK", "t", "z"),
                     level = 100 * x$level.ma,
                     measure = "OR", control = x$control,
                     data = data.frame(subgroup.rma),
                     ...)))
    }
    ##
    if (length(TE.common.w) != length(as.numeric(glmm.common$b))) {
      TE.common.w[!is.na(TE.common.w)] <- as.numeric(glmm.common$b)
      seTE.common.w[!is.na(seTE.common.w)] <- as.numeric(glmm.common$se)
      TE.random.w[!is.na(TE.random.w)] <- as.numeric(glmm.random$b)
      seTE.random.w[!is.na(seTE.random.w)] <- as.numeric(glmm.random$se)
    }
    else {
      TE.common.w   <- as.numeric(glmm.common$b)
      seTE.common.w <- as.numeric(glmm.common$se)
      TE.random.w   <- as.numeric(glmm.random$b)
      seTE.random.w <- as.numeric(glmm.random$se)
    }
    ##
    tau2.w <- rep_len(glmm.random$tau2, n.levs)
    tau.w <- sqrt(tau2.w)
    ##
    Rb.w     <- rep_len(NA, n.levs)
    Rb.w.low <- rep_len(NA, n.levs)
    Rb.w.upp <- rep_len(NA, n.levs)
    ##
    ci.common.w  <- ci(TE.common.w, seTE.common.w, x$level.ma)
    ##
    if (!is.null(x$method.random.ci == "HK") && x$method.random.ci == "HK")
      ci.random.w <- ci(TE.random.w, seTE.random.w, x$level.ma, df = k.w - 1)
    else
      ci.random.w <- ci(TE.random.w, seTE.random.w, x$level.ma)
    ##
    lower.common.w <- ci.common.w$lower
    upper.common.w <- ci.common.w$upper
    statistic.common.w <- ci.common.w$statistic
    pval.common.w <- ci.common.w$p
    ##
    lower.random.w <- ci.random.w$lower
    upper.random.w <- ci.random.w$upper
    statistic.random.w <- ci.random.w$statistic
    pval.random.w <- ci.random.w$p
    ##
    ## Tests for subgroup differences
    ##
    Q.w.common <- glmm.common.Q$QE.Wld
    df.Q.w <- glmm.common.Q$QE.df
    pval.Q.w.common <- glmm.common.Q$QEp.Wld
    ##
    Q.b.common <- glmm.common.Q$QM
    Q.b.random <- glmm.random.Q$QM
    ##
    df.Q.b <- glmm.common.Q$QMdf
    df.Q.b <- df.Q.b[!is.na(df.Q.b)]
    df.Q.b.common <- df.Q.b
    df.Q.b.random <- df.Q.b
    ##
    pval.Q.b.common <- glmm.common.Q$QMp
    pval.Q.b.random <- glmm.random.Q$QMp
  }
  
  
  ##
  ##
  ## (4) Return list with results of subgroup analysis
  ##
  ##
  res <- list(subgroup.levels = levs,
              ##
              k.w = k.w,
              k.study.w = k.study.w,
              k.all.w = k.all.w,
              k.TE.w = k.TE.w,
              n.w = n.w,
              event.w = event.w,
              ##
              TE.common.w = TE.common.w,
              seTE.common.w = seTE.common.w,
              statistic.common.w = statistic.common.w,
              pval.common.w = pval.common.w,
              lower.common.w = lower.common.w,
              upper.common.w = upper.common.w,
              w.common.w = w.common.w,
              ##
              TE.random.w = TE.random.w,
              seTE.random.w = seTE.random.w,
              statistic.random.w = statistic.random.w,
              pval.random.w = pval.random.w,
              df.random.w = df.random.w,
              lower.random.w = lower.random.w,
              upper.random.w = upper.random.w,
              w.random.w = w.random.w,
              ##
              seTE.classic.w = seTE.classic.w,
              ##
              df.hakn.ci.w = df.hakn.ci.w,
              seTE.hakn.ci.w = seTE.hakn.ci.w,
              seTE.hakn.adhoc.ci.w = seTE.hakn.adhoc.ci.w,
              ##
              df.kero.w = df.kero.w,
              seTE.kero.w = seTE.kero.w,
              ##
              seTE.predict.w = seTE.predict.w,
              df.predict.w = df.predict.w,
              lower.predict.w = lower.predict.w,
              upper.predict.w = upper.predict.w,
              seTE.hakn.pi.w = seTE.hakn.pi.w,
              seTE.hakn.adhoc.pi.w = seTE.hakn.adhoc.pi.w,
              ##
              Q.w = Q.w,
              pval.Q.w = pvalQ(Q.w, k.w - 1),
              ##
              tau2.w = tau2.w,
              tau.w = tau.w,
              ##
              H.w = H.w,
              lower.H.w = lower.H.w,
              upper.H.w = upper.H.w,
              ##
              I2.w = I2.w,
              lower.I2.w = lower.I2.w,
              upper.I2.w = upper.I2.w,
              ##
              Rb.w = Rb.w,
              lower.Rb.w = lower.Rb.w,
              upper.Rb.w = upper.Rb.w,
              ##
              Q.w.common = Q.w.common,
              Q.w.random = NA,
              df.Q.w = df.Q.w,
              pval.Q.w.common = pval.Q.w.common,
              pval.Q.w.random = NA,
              ##
              Q.b.common = Q.b.common,
              Q.b.random = Q.b.random,
              df.Q.b = df.Q.b,
              df.Q.b.common = df.Q.b.common,
              df.Q.b.random = df.Q.b.random,
              pval.Q.b.common = pval.Q.b.common,
              pval.Q.b.random = pval.Q.b.random
              )
  ##
  ## No general list elements
  ##
  res$n.e.w <- n.e.w
  res$n.c.w <- n.c.w
  res$n.harmonic.mean.w <- n.harmonic.mean.w
  ##
  res$event.e.w <- event.e.w
  res$event.c.w <- event.c.w
  ##
  res$time.e.w <- time.e.w
  res$time.c.w <- time.c.w
  res$t.harmonic.mean.w <- t.harmonic.mean.w
  ##
  ## Deprecated list elements
  ##
  res$zval.common.w <- statistic.common.w
  res$zval.random.w <- statistic.random.w
  ##
  ## Backward compatibility
  ##
  res$TE.fixed.w <- res$TE.common.w
  res$seTE.fixed.w <- res$seTE.common.w
  res$lower.fixed.w <- res$lower.common.w
  res$upper.fixed.w <- res$upper.common.w
  res$statistic.fixed.w <- res$statistic.common.w
  res$pval.fixed.w <- res$pval.common.w
  res$zval.fixed.w <- res$zval.common.w
  res$w.fixed.w <- res$w.common.w
  ##
  res$Q.w.fixed <- res$Q.w.common
  res$pval.Q.w.fixed <- res$pval.Q.w.common
  res$Q.b.fixed <- res$Q.b.common
  res$pval.Q.b.fixed <- res$pval.Q.b.common
  ##
  res$bylevs <- res$subgroup.levels
  ##
  res$seTE.hakn.w <- res$seTE.hakn.ci.w
  res$seTE.hakn.adhoc.w <- res$seTE.hakn.adhoc.ci.w
  
  
  res
}
