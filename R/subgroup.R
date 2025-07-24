subgroup <- function(x, tau.preset = NULL, subgroup.rma,
                     seed = NULL, ...) {
  
  subgroup <- x$subgroup
  method.predict <-
    replaceVal(x$method.predict, "", gs("method.predict"))
  ##
  levs <- bylevs(subgroup)
  n.levs <- length(levs)
  ##
  methci <- paste(x$method.random.ci,
                  toupper(substring(x$adhoc.hakn.ci, 1, 2)),
                  sep = "-")
  methci <- gsub("-$", "", methci)
  ##
  methpi <- paste(method.predict,
                  toupper(substring(x$adhoc.hakn.pi, 1, 2)),
                  sep = "-")
  methpi <- gsub("-$", "", methpi)
  ##
  if (!(length(subgroup) > 0)) {
    warning("Argument 'subgroup' is missing.")
    return(NULL)
  }
  ##
  if (!is.null(seed) & any(method.predict %in% "NNF")) {
    if (length(seed) == 1)
      seed <- rep_len(seed, n.levs)
    chknumeric(seed, length = n.levs, name = "seed.predict.subgroup")
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
  cor.prop.mean.rate <- cor | prop | mean | rate
  ##
  is.glmm <- x$method == "GLMM"
  is.mlm <- !is.null(x$three.level) && x$three.level
  ##
  sumNA <- function(x, exclude = NULL) {
    if (all(is.na(x)))
      return(NA)
    else {
      if (is.null(exclude))
        return(sum(x, na.rm = TRUE))
      else
        return(sum(x[!exclude], na.rm = TRUE))
    }
  }
  
  
  ##
  ##
  ## (1) Subgroup analysis without common tau2
  ##
  ##
  res.i <- vector(mode = "list")
  j <- 0
  ##
  for (i in levs) {
    j <- j + 1
    sel <- subgroup == i
    ##
    if (all(is.na(x$studlab[sel])))
      stop("No data available for subgroup = ", i)
    ##
    if (!is.null(seed))
      seed.i <- seed[j]
    else
      seed.i <- NULL
    ##
    if (bin) {
      if (x$method.incr == "user") {
        incr <- NULL
        incr.e <- x$incr.e[sel]
        incr.c <- x$incr.c[sel]
      }
      else {
        incr <- x$incr
        if (length(incr) > 1)
          incr <- incr[sel]
        #
        incr.e <- NULL
        incr.c <- NULL
      }
      #
      meta1 <- metabin(x$event.e[sel], x$n.e[sel],
                       x$event.c[sel], x$n.c[sel],
                       studlab = x$studlab[sel],
                       exclude = x$exclude[sel],
                       cluster =
                         if (!is.null(x$cluster)) x$cluster[sel] else NULL,
                       rho = x$rho,
                       ##
                       method = x$method,
                       sm = x$sm,
                       #
                       incr = incr,
                       incr.e = incr.e, incr.c = incr.c,
                       method.incr = x$method.incr, allstudies = x$allstudies,
                       #
                       MH.exact = x$MH.exact,
                       RR.Cochrane = x$RR.Cochrane,
                       Q.Cochrane = x$Q.Cochrane,
                       model.glmm = x$model.glmm,
                       ##
                       level.ma = x$level.ma,
                       method.common.ci = x$method.common.ci,
                       method.random.ci = x$method.random.ci,
                       adhoc.hakn.ci = x$adhoc.hakn.ci,
                       ##
                       level.predict = x$level.predict,
                       method.predict = method.predict,
                       adhoc.hakn.pi = x$adhoc.hakn.pi,
                       seed.predict = seed.i,
                       ##
                       method.tau = x$method.tau,
                       method.tau.ci = x$method.tau.ci,
                       level.hetstat = x$level.hetstat,
                       tau.preset = tau.preset,
                       TE.tau = x$TE.tau,
                       ##
                       keepdata = FALSE,
                       warn = x$warn,
                       control = x$control)
    }
    #
    else if (cont)
      meta1 <- metacont(x$n.e[sel], x$mean.e[sel],
                        x$sd.e[sel],
                        x$n.c[sel], x$mean.c[sel],
                        x$sd.c[sel],
                        studlab = x$studlab[sel],
                        exclude = x$exclude[sel],
                        cluster =
                          if (!is.null(x$cluster)) x$cluster[sel] else NULL,
                        rho = x$rho,
                        #
                        median.e = subsetVar(x$median.e, sel),
                        q1.e = subsetVar(x$q1.e, sel),
                        q3.e = subsetVar(x$q3.e, sel),
                        min.e = subsetVar(x$min.e, sel),
                        max.e = subsetVar(x$max.e, sel),
                        #
                        median.c = subsetVar(x$median.c, sel),
                        q1.c = subsetVar(x$q1.c, sel),
                        q3.c = subsetVar(x$q3.c, sel),
                        min.c = subsetVar(x$min.c, sel),
                        max.c = subsetVar(x$max.c, sel),
                        #
                        method.mean = subsetVar(x$method.mean, sel),
                        method.sd = subsetVar(x$method.sd, sel),
                        #
                        approx.mean.e = subsetVar(x$approx.mean.e, sel),
                        approx.sd.e = subsetVar(x$approx.sd.e, sel),
                        #
                        approx.mean.c = subsetVar(x$approx.mean.c, sel),
                        approx.sd.c = subsetVar(x$approx.sd.c, sel),
                        #
                        sm = x$sm,
                        ##
                        pooledvar = x$pooledvar,
                        method.smd = x$method.smd,
                        sd.glass = x$sd.glass,
                        exact.smd = x$exact.smd,
                        ##
                        level.ma = x$level.ma,
                        method.common.ci = x$method.common.ci,
                        method.random.ci = x$method.random.ci,
                        adhoc.hakn.ci = x$adhoc.hakn.ci,
                        ##
                        level.predict = x$level.predict,
                        method.predict = method.predict,
                        adhoc.hakn.pi = x$adhoc.hakn.pi,
                        seed.predict = seed.i,
                        ##
                        method.tau = x$method.tau,
                        method.tau.ci = x$method.tau.ci,
                        level.hetstat = x$level.hetstat,
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
                       rho = x$rho,
                       ##
                       sm = x$sm,
                       ##
                       level.ma = x$level.ma,
                       method.common.ci = x$method.common.ci,
                       method.random.ci = x$method.random.ci,
                       adhoc.hakn.ci = x$adhoc.hakn.ci,
                       ##
                       level.predict = x$level.predict,
                       method.predict = method.predict,
                       adhoc.hakn.pi = x$adhoc.hakn.pi,
                       seed.predict = seed.i,
                       ##
                       method.tau = x$method.tau,
                       method.tau.ci = x$method.tau.ci,
                       level.hetstat = x$level.hetstat,
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
                       rho = x$rho,
                       #
                       weights.common =
                         if (!is.null(x$weights.common)) x$weights.common[sel] else NULL,
                       #
                       weights.random =
                         if (!is.null(x$weights.random)) x$weights.random[sel] else NULL,
                       #
                       sm = x$sm,
                       ##
                       level.ma = x$level.ma,
                       method.common.ci = x$method.common.ci,
                       method.random.ci = x$method.random.ci,
                       adhoc.hakn.ci = x$adhoc.hakn.ci,
                       ##
                       level.predict = x$level.predict,
                       method.predict = method.predict,
                       adhoc.hakn.pi = x$adhoc.hakn.pi,
                       seed.predict = seed.i,
                       ##
                       method.tau = x$method.tau,
                       method.tau.ci = x$method.tau.ci,
                       level.hetstat = x$level.hetstat,
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
    else if (inc) {
      if (x$method.incr == "user") {
        incr <- NULL
        incr.e <- x$incr.e[sel]
        incr.c <- x$incr.c[sel]
      }
      else {
        incr <- x$incr
        if (length(incr) > 1)
          incr <- incr[sel]
        #
        incr.e <- NULL
        incr.c <- NULL
      }
      #
      meta1 <- metainc(x$event.e[sel], x$time.e[sel],
                       x$event.c[sel], x$time.c[sel],
                       studlab = x$studlab[sel],
                       exclude = x$exclude[sel],
                       cluster =
                         if (!is.null(x$cluster)) x$cluster[sel] else NULL,
                       rho = x$rho,
                       ##
                       method = x$method,
                       sm = x$sm,
                       #
                       incr = incr,
                       incr.e = incr.e, incr.c = incr.c,
                       method.incr = x$method.incr,
                       #
                       model.glmm = x$model.glmm,
                       ##
                       level.ma = x$level.ma,
                       method.common.ci = x$method.common.ci,
                       method.random.ci = x$method.random.ci,
                       adhoc.hakn.ci = x$adhoc.hakn.ci,
                       ##
                       level.predict = x$level.predict,
                       method.predict = method.predict,
                       adhoc.hakn.pi = x$adhoc.hakn.pi,
                       seed.predict = seed.i,
                       ##
                       method.tau = x$method.tau,
                       method.tau.ci = x$method.tau.ci,
                       level.hetstat = x$level.hetstat,
                       tau.preset = tau.preset,
                       TE.tau = x$TE.tau,
                       ##
                       n.e = x$n.e[sel], n.c = x$n.c[sel],
                       ##
                       keepdata = FALSE,
                       warn = x$warn,
                       control = x$control)
    }
    #
    else if (mean)
      meta1 <- metamean(x$n[sel], x$mean[sel], x$sd[sel],
                        studlab = x$studlab[sel],
                        exclude = x$exclude[sel],
                        cluster = subsetVar(x$cluster, sel),
                        rho = x$rho,
                        #
                        median = subsetVar(x$median, sel),
                        q1 = subsetVar(x$q1, sel),
                        q3 = subsetVar(x$q3, sel),
                        min = subsetVar(x$min, sel),
                        max = subsetVar(x$max, sel),
                        method.mean = subsetVar(x$method.mean, sel),
                        method.sd = subsetVar(x$method.sd, sel),
                        approx.mean = subsetVar(x$approx.mean, sel),
                        approx.sd = subsetVar(x$approx.sd, sel),
                        #
                        sm = x$sm,
                        ##
                        level.ma = x$level.ma,
                        method.common.ci = x$method.common.ci,
                        method.random.ci = x$method.random.ci,
                        adhoc.hakn.ci = x$adhoc.hakn.ci,
                        ##
                        level.predict = x$level.predict,
                        method.predict = method.predict,
                        adhoc.hakn.pi = x$adhoc.hakn.pi,
                        seed.predict = seed.i,
                        ##
                        method.tau = x$method.tau,
                        method.tau.ci = x$method.tau.ci,
                        level.hetstat = x$level.hetstat,
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
                        rho = x$rho,
                        ##
                        method = x$method,
                        sm = x$sm,
                        incr = if (length(x$incr) == 1) x$incr else x$incr[sel],
                        method.incr = x$method.incr,
                        ##
                        level.ma = x$level.ma,
                        method.common.ci = x$method.common.ci,
                        method.random.ci = x$method.random.ci,
                        adhoc.hakn.ci = x$adhoc.hakn.ci,
                        ##
                        level.predict = x$level.predict,
                        method.predict = method.predict,
                        adhoc.hakn.pi = x$adhoc.hakn.pi,
                        seed.predict = seed.i,
                        ##
                        method.tau = x$method.tau,
                        method.tau.ci = x$method.tau.ci,
                        level.hetstat = x$level.hetstat,
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
                        rho = x$rho,
                        ##
                        method = x$method,
                        sm = x$sm,
                        incr = if (length(x$incr) == 1) x$incr else x$incr[sel],
                        method.incr = x$method.incr,
                        ##
                        level.ma = x$level.ma,
                        method.common.ci = x$method.common.ci,
                        method.random.ci = x$method.random.ci,
                        adhoc.hakn.ci = x$adhoc.hakn.ci,
                        ##
                        level.predict = x$level.predict,
                        method.predict = method.predict,
                        adhoc.hakn.pi = x$adhoc.hakn.pi,
                        seed.predict = seed.i,
                        ##
                        method.tau = x$method.tau,
                        method.tau.ci = x$method.tau.ci,
                        level.hetstat = x$level.hetstat,
                        tau.preset = tau.preset,
                        TE.tau = x$TE.tau,
                        ##
                        n = if (!is.null(x$n)) x$n[sel] else NULL,
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
                       lower.tau2 =
                         if (meta1$three.level) NA else meta1$lower.tau2,
                       upper.tau2 =
                         if (meta1$three.level) NA else meta1$upper.tau2,
                       #
                       tau = sum(meta1$tau),
                       lower.tau =
                         if (meta1$three.level) NA else meta1$lower.tau,
                       upper.tau =
                         if (meta1$three.level) NA else meta1$upper.tau,
                       #
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
                       event =
                         if (prop) sumNA(meta1$event, meta1$exclude) else NA,
                       n =
                         if (cor.prop.mean.rate) sumNA(meta1$n, meta1$exclude) else NA,
                       ##
                       event.e =
                         if (bin.inc) sumNA(meta1$event.e, meta1$exclude) else NA,
                       n.e =
                         if (bin.cont.gen) sumNA(meta1$n.e, meta1$exclude) else NA,
                       event.c =
                         if (bin.inc) sumNA(meta1$event.c, meta1$exclude) else NA,
                       n.c =
                         if (bin.cont.gen) sumNA(meta1$n.c, meta1$exclude) else NA,
                       ##
                       time.e =
                         if (inc) sumNA(meta1$time.e, meta1$exclude) else NA,
                       time.c =
                         if (inc) sumNA(meta1$time.c, meta1$exclude) else NA,
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
  Q.w <- extrVec(res.i, "Q", levs, TRUE)
  ##
  tau2.w <- extrVec(res.i, "tau2", levs)
  lower.tau2.w <- extrVec(res.i, "lower.tau2", levs)
  upper.tau2.w <- extrVec(res.i, "upper.tau2", levs)
  #
  tau.w <- extrVec(res.i, "tau", levs)
  lower.tau.w <- extrVec(res.i, "lower.tau", levs)
  upper.tau.w <- extrVec(res.i, "upper.tau", levs)
  #
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
        metagen(TE.random.w,
                lower = lower.random.w[, i],
                upper = upper.random.w[, i],
                method.tau = "DL")$Q
  }
  else {
    n.random <- 1
    Q.b.random <-
      metagen(TE.random.w, seTE.random.w, method.tau = "DL")$Q
  }
  ##
  df.Q.b <- NA
  df.Q.b.common <- NA
  df.Q.b.random <- list()
  ##
  pval.Q.b.common <- NA
  pval.Q.b.random <- NA
  ##
  if (missing(subgroup.rma)) {
    df.Q.b <- if (x$k == 0) 0 else x$k - 1 - sum((k.w - 1)[!is.na(Q.w)])
    df.Q.b.common <- df.Q.b
    pval.Q.b.common <- pvalQ(Q.b.common, df.Q.b.common)
    ##
    if (n.random > 1) {
      pval.Q.b.random <- vector("numeric", n.random)
      df.Q.b.random <- vector("list", n.random)
      for (i in seq_len(n.random)) {
        df.Q.b.random[[i]] <- df.Q.b
        pval.Q.b.random[i] <- pvalQ(Q.b.random[i], df.Q.b.random[[i]])
      }
      ##
      names(df.Q.b.random) <- names(pval.Q.b.random) <-
        colnames(seTE.random.w)
    }
    else {
      df.Q.b.random <- df.Q.b
      pval.Q.b.random <- pvalQ(Q.b.random, df.Q.b.random)
    }
  }
  
  
  ##
  ##
  ## (2) Subgroup analysis with common tau-squared (three-level model)
  ##
  ##
  
  if (!missing(subgroup.rma) && is.mlm) {
    n.methci <- length(x$method.random.ci)
    ##
    mod <- as.call(~ subgroup.rma - 1)
    mod.Q <- as.call(~ subgroup.rma)
    ##
    if (is.null(x$exclude))
      x$exclude <- rep(FALSE, length(x$TE))
    ##
    list.mlm <- list(yi = x$TE[!x$exclude],
                     V = vcalc(vi = x$seTE[!x$exclude]^2,
                               cluster = x$cluster[!x$exclude],
                               obs = seq_along(x$cluster)[!x$exclude],
                               rho = replaceNULL(x$rho, 0)))
    ##
    mlm <-
      runMLM(c(list.mlm,
               list(data = 
                      data.frame(subgroup.rma = subgroup.rma,
                                 cluster = x$cluster[!x$exclude],
                                 idx = seq_along(x$cluster)[!x$exclude]))),
             method.tau = x$method.tau,
             method.random.ci = x$method.random.ci,
             level = x$level.ma,
             mods = mod,
             control = list(x$control),
             warn = FALSE)
    ##
    ## Random effects model
    ##
    sel.r <- !is.na(TE.random.w)
    TE.random.w[sel.r] <- as.numeric(mlm[[1]]$b)
    seTE.random.w[sel.r] <- as.numeric(mlm[[1]]$se)
    ##
    if (n.methci == 1) {
      lower.random.w[sel.r] <- as.numeric(mlm[[1]]$ci.lb)
      upper.random.w[sel.r] <- as.numeric(mlm[[1]]$ci.ub)
      statistic.random.w[sel.r] <- as.numeric(mlm[[1]]$zval)
      pval.random.w[sel.r] <- as.numeric(mlm[[1]]$pval)
    }
    else {
      for (i in seq_len(n.methci)) {
        lower.random.w[sel.r, i] <- as.numeric(mlm[[i]]$ci.lb)
        upper.random.w[sel.r, i] <- as.numeric(mlm[[i]]$ci.ub)
        statistic.random.w[sel.r, i] <- as.numeric(mlm[[i]]$zval)
        pval.random.w[sel.r, i] <- as.numeric(mlm[[i]]$pval)
      }
    }
    ##
    ## Heterogeneity measures
    ##
    tau2.w <- rep_len(sum(mlm[[1]]$sigma2), n.levs)
    lower.tau2.w <- tau2.w
    upper.tau2.w <- tau2.w
    lower.tau2.w[!is.na(lower.tau2.w)] <- NA
    upper.tau2.w[!is.na(upper.tau2.w)] <- NA
    #
    tau.w <- sqrt(tau2.w)
    lower.tau.w <- tau.w
    upper.tau.w <- tau.w
    lower.tau.w[!is.na(lower.tau.w)] <- NA
    upper.tau.w[!is.na(upper.tau.w)] <- NA
    #
    upper.Rb.w <- lower.Rb.w <- Rb.w <- rep_len(NA, n.levs)
    ##
    ## Prediction interval
    ##
    tau2.calc <- if (is.na(sum(mlm[[1]]$sigma2))) 0 else tau2.w
    seTE.predict.w <- sqrt(seTE.random.w^2 + tau2.w)
    ##
    n.methpi <- length(method.predict)
    ##
    for (i in seq_len(n.methpi)) {
      df.i <- set_df_predict(method.predict[i], k.w)
      #
      ci.p <- ci(TE.random.w, seTE.predict.w, x$level.predict, df = df.i)
      ##
      if (n.methpi == 1) {
        lower.predict.w <- ci.p$lower
        upper.predict.w <- ci.p$upper
      }
      else {
        lower.predict.w[, i] <- ci.p$lower
        upper.predict.w[, i] <- ci.p$upper
      }
      ##
      sel.p <-
        !(method.predict[i] == "V" & k.w >= 2 |
            method.predict[i] == "HTS" & k.w >= 3 |
            method.predict[i] == "S")
      #
      if (n.methpi == 1) {
        lower.predict.w[sel.p] <- NA
        upper.predict.w[sel.p] <- NA
      }
      else {
        lower.predict.w[sel.p, i] <- NA
        upper.predict.w[sel.p, i] <- NA
      }
    }
    ##
    seTE.predict <- seTE.predict.w
    lower.predict <- lower.predict.w
    upper.predict <- upper.predict.w
    ##
    if (length(method.predict) > 1)
      names(lower.predict) <- names(upper.predict) <-
        methpi
    ##
    ## Test for subgroup differences
    ##
    mlm.Q <-
      runMLM(c(list.mlm,
               list(data = 
                      data.frame(subgroup.rma = subgroup.rma,
                                 cluster = x$cluster[!x$exclude],
                                 idx = seq_along(x$cluster)[!x$exclude]))),
             method.tau = x$method.tau,
             method.random.ci = x$method.random.ci,
             level = x$level.ma,
             mods = mod.Q,
             control = list(x$control),
             warn = FALSE)
    ##
    ## Tests for subgroup differences
    ##
    for (i in seq_len(n.methci)) {
      Q.b.random[i] <- mlm.Q[[i]]$QM
      df.i <- mlm.Q[[i]]$QMdf
      df.Q.b.random[[i]] <- df.i[!is.na(df.i)]
      pval.Q.b.random[i] <- mlm.Q[[i]]$QMp
      ##
      if (n.methci == 1)
        df.Q.b.random <- df.Q.b.random[[1]]
    }
    ##
    Q.w.common <- NA
    Q.w.random <- NA
    df.Q.w <- NA
    pval.Q.w.common <- NA
    pval.Q.w.random <- NA
  }
  
  
  ##
  ##
  ## (3) Subgroup analysis with common tau-squared (GLMM)
  ##
  ##
  
  if (!missing(subgroup.rma) && is.glmm) {
    n.methci <- length(x$method.random.ci)
    ##
    mod <- as.call(~ subgroup.rma - 1)
    mod.Q <- as.call(~ subgroup.rma)
    ##
    if (is.null(x$exclude))
      x$exclude <- rep(FALSE, length(x$TE))
    ##
    if (prop) {
      list.glmm <-
        list(xi = x$event[!x$exclude], ni = x$n[!x$exclude],
             measure = "PLO")
      ##
      use.random <-
        sum(!x$exclude) > 1 &
        sum(x$event[!x$exclude], na.rm = TRUE) > 0 &
        any(x$event[!x$exclude] != x$n[!x$exclude])
    }
    ##
    else if (rate) {
      list.glmm <-
        list(xi = x$event[!x$exclude], ti = x$time[!x$exclude],
             measure = "IRLN")
      ##
      use.random <-
        sum(!x$exclude) > 1 &
        sum(x$event[!x$exclude], na.rm = TRUE) > 0 &
        any(x$event[!x$exclude] != x$time[!x$exclude])
    }
    ##
    else if (inc) {
      list.glmm <-
        list(x1i = x$event.e[!x$exclude], t1i = x$time.e[!x$exclude],
             x2i = x$event.c[!x$exclude], t2i = x$time.c[!x$exclude],
             measure = "IRR", model = x$model.glmm)
      ##
      use.random <-
        sum(!x$exclude) > 1 &
        !((sum(x$event.e[!x$exclude], na.rm = TRUE) == 0 &
           sum(x$event.c[!x$exclude], na.rm = TRUE) == 0) |
          (!any(x$event.e[!x$exclude] != x$time.e[!x$exclude]) |
           !any(x$event.c[!x$exclude] != x$time.c[!x$exclude])))
    }
    ##
    else if (bin) {
      list.glmm <-
        list(ai = x$event.e[!x$exclude], n1i = x$n.e[!x$exclude],
             ci = x$event.c[!x$exclude], n2i = x$n.c[!x$exclude],
             measure = "OR", model = x$model.glmm)
      ##
      use.random <-
        sum(!x$exclude) > 1 &
        !((sum(x$event.e[!x$exclude], na.rm = TRUE) == 0 &
           sum(x$event.c[!x$exclude], na.rm = TRUE) == 0) |
          (!any(x$event.e[!x$exclude] != x$n.e[!x$exclude]) |
           !any(x$event.c[!x$exclude] != x$n.c[!x$exclude])))
    }
    ##
    glmms <-
      runGLMM(list.glmm,
              method.tau = x$method.tau,
              method.random.ci = x$method.random.ci,
              level = x$level.ma,
              data = list(data = data.frame(subgroup.rma)),
              mods = mod,
              control = list(x$control),
              use.random = use.random,
              warn = FALSE)
    ##
    glmm.c <- glmms$glmm.common
    glmm.r <- glmms$glmm.random
    ##
    ## Common effect model
    ##
    sel.c <- !is.na(TE.common.w)
    TE.common.w[sel.c] <- as.numeric(glmm.c$b)
    seTE.common.w[sel.c] <- as.numeric(glmm.c$se)
    lower.common.w[sel.c] <- as.numeric(glmm.c$ci.lb)
    upper.common.w[sel.c] <- as.numeric(glmm.c$ci.ub)
    statistic.common.w[sel.c] <- as.numeric(glmm.c$zval)
    pval.common.w[sel.c] <- as.numeric(glmm.c$pval)
    ##
    ## Random effects model
    ##
    sel.r <- !is.na(TE.random.w)
    TE.random.w[sel.r] <- as.numeric(glmm.r[[1]]$b)
    seTE.random.w[sel.r] <- as.numeric(glmm.r[[1]]$se)
    ##
    if (n.methci == 1) {
      lower.random.w[sel.r] <- as.numeric(glmm.r[[1]]$ci.lb)
      upper.random.w[sel.r] <- as.numeric(glmm.r[[1]]$ci.ub)
      statistic.random.w[sel.r] <- as.numeric(glmm.r[[1]]$zval)
      pval.random.w[sel.r] <- as.numeric(glmm.r[[1]]$pval)
    }
    else {
      for (i in seq_len(n.methci)) {
        lower.random.w[sel.r, i] <- as.numeric(glmm.r[[i]]$ci.lb)
        upper.random.w[sel.r, i] <- as.numeric(glmm.r[[i]]$ci.ub)
        statistic.random.w[sel.r, i] <- as.numeric(glmm.r[[i]]$zval)
        pval.random.w[sel.r, i] <- as.numeric(glmm.r[[i]]$pval)
      }
    }
    ##
    ## Heterogeneity measures
    ##
    tau2.w <- rep_len(glmm.r[[1]]$tau2, n.levs)
    lower.tau2.w <- tau2.w
    upper.tau2.w <- tau2.w
    lower.tau2.w[!is.na(lower.tau2.w)] <- NA
    upper.tau2.w[!is.na(upper.tau2.w)] <- NA
    #
    tau.w <- sqrt(tau2.w)
    lower.tau.w <- tau.w
    upper.tau.w <- tau.w
    lower.tau.w[!is.na(lower.tau.w)] <- NA
    upper.tau.w[!is.na(upper.tau.w)] <- NA
    #
    upper.Rb.w <- lower.Rb.w <- Rb.w <- rep_len(NA, n.levs)
    ##
    ## Prediction interval
    ##
    tau2.calc <- if (is.na(glmm.r[[1]]$tau2)) 0 else tau2.w
    seTE.predict.w <- sqrt(seTE.random.w^2 + tau2.w)
    ##
    n.methpi <- length(method.predict)
    ##
    for (i in seq_len(n.methpi)) {
      df.i <- set_df_predict(method.predict[i], k.w)
      #
      ci.p <- ci(TE.random.w, seTE.predict.w, x$level.predict, df = df.i)
      ##
      if (n.methpi == 1) {
        lower.predict.w <- ci.p$lower
        upper.predict.w <- ci.p$upper
      }
      else {
        lower.predict.w[, i] <- ci.p$lower
        upper.predict.w[, i] <- ci.p$upper
      }
      ##
      sel.p <-
        !(method.predict[i] == "V" & k.w >= 2 |
            method.predict[i] == "HTS" & k.w >= 3 |
            method.predict[i] == "S")
      #
      if (n.methpi == 1) {
        lower.predict.w[sel.p] <- NA
        upper.predict.w[sel.p] <- NA
      }
      else {
        lower.predict.w[sel.p, i] <- NA
        upper.predict.w[sel.p, i] <- NA
      }
    }
    ##
    seTE.predict <- seTE.predict.w
    lower.predict <- lower.predict.w
    upper.predict <- upper.predict.w
    ##
    if (length(method.predict) > 1)
      names(lower.predict) <- names(upper.predict) <-
        methpi
    ##
    ## Test for subgroup differences
    ##
    glmms.Q <-
      runGLMM(list.glmm,
              method.tau = x$method.tau,
              method.random.ci = x$method.random.ci,
              level = x$level.ma,
              data = list(data = data.frame(subgroup.rma)),
              mods = mod.Q,
              control = list(x$control),
              use.random = use.random,
              warn = FALSE)
    ##
    glmm.c.Q <- glmms.Q$glmm.common
    glmm.r.Q <- glmms.Q$glmm.random
    ##
    ## Tests for subgroup differences
    ##
    Q.b.common <- glmm.c.Q$QM
    df.Q.b.common <- glmm.c.Q$QMdf[!is.na(glmm.c.Q$QMdf)]
    pval.Q.b.common <- glmm.c.Q$QMp
    ##
    for (i in seq_len(n.methci)) {
      Q.b.random[i] <- glmm.r.Q[[i]]$QM
      df.i <- glmm.r.Q[[i]]$QMdf
      df.Q.b.random[[i]] <- df.i[!is.na(df.i)]
      pval.Q.b.random[i] <- glmm.r.Q[[i]]$QMp
      ##
      if (n.methci == 1)
        df.Q.b.random <- df.Q.b.random[[1]]
    }
    ##
    Q.w.common <- NA
    Q.w.random <- NA
    df.Q.w <- NA
    pval.Q.w.common <- NA
    pval.Q.w.random <- NA
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
              lower.tau2.w = lower.tau2.w,
              upper.tau2.w = upper.tau2.w,
              #
              tau.w = tau.w,
              lower.tau.w = lower.tau.w,
              upper.tau.w = upper.tau.w,
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
              pval.Q.b.random = pval.Q.b.random,
              ##
              seed.predict.subgroup = seed
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
  ## Set common effect results to NA for three-level model
  ##
  if (is.mlm) {
    res$TE.common.w[!is.na(res$TE.common.w)] <- NA
    res$seTE.common.w[!is.na(res$seTE.common.w)] <- NA
    res$statistic.common.w[!is.na(res$statistic.common.w)] <- NA
    res$pval.common.w[!is.na(res$pval.common.w)] <- NA
    res$lower.common.w[!is.na(res$lower.common.w)] <- NA
    res$upper.common.w[!is.na(res$upper.common.w)] <- NA
    res$w.common.w[!is.na(res$w.common.w)] <- NA
    res$Q.w.common[!is.na(res$Q.w.common)] <- NA
    res$pval.Q.w.common[!is.na(res$pval.Q.w.common)] <- NA
    res$Q.b.common[!is.na(res$Q.b.common)] <- NA
    res$df.Q.b.common[!is.na(res$df.Q.b.common)] <- NA
    res$pval.Q.b.common[!is.na(res$pval.Q.b.common)] <- NA
  }
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
