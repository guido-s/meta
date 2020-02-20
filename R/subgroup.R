subgroup <- function(x, tau.preset = NULL, byvar.glmm, ...) {
  
  
  byvar <- x$byvar
  ##
  bylevs <- bylevs(byvar)
  n.bylevs <- length(bylevs)
  
  
  if (!(length(byvar) > 0)) {
    warning("Argument 'byvar' is missing.")
    return(NULL)
  }
  
  
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
  
  
  sumNA <- function(x)
    if (all(is.na(x)))
      NA
    else
      sum(x, na.rm = TRUE)
  
  
  res.w <- matrix(NA, ncol = 34, nrow = n.bylevs)
  add.w <- matrix(NA, ncol =  2, nrow = n.bylevs)
  j <- 0
  ##
  for (i in bylevs) {
    j <- j+1
    sel <- byvar == i
    ##
    if (all(is.na(x$studlab[sel])))
      stop("No data available for byvar = ", i)
    ##
    ##
    if (bin)
      meta1 <- metabin(x$event.e[sel], x$n.e[sel],
                       x$event.c[sel], x$n.c[sel],
                       studlab = x$studlab[sel],
                       exclude = x$exclude[sel],
                       method = x$method,
                       sm = x$sm,
                       incr = if (length(x$incr) == 1) x$incr else x$incr[sel],
                       allincr = x$allincr,
                       addincr = x$addincr,
                       allstudies = x$allstudies,
                       MH.exact = x$MH.exact,
                       RR.Cochrane = x$RR.Cochrane,
                       Q.Cochrane = x$Q.Cochrane,
                       level = x$level, level.comb = x$level.comb,
                       hakn = x$hakn,
                       method.tau = x$method.tau,
                       method.tau.ci = x$method.tau.ci,
                       tau.preset = tau.preset,
                       TE.tau = x$TE.tau,
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
                        sm = x$sm, pooledvar = x$pooledvar,
                        level = x$level, level.comb = x$level.comb,
                        hakn = x$hakn,
                        method.tau = x$method.tau,
                        method.tau.ci = x$method.tau.ci,
                        tau.preset = tau.preset, TE.tau = x$TE.tau,
                        warn = x$warn)
    ##
    else if (cor)
      meta1 <- metacor(x$cor[sel], x$n[sel],
                       sm = x$sm,
                       studlab = x$studlab[sel],
                       exclude = x$exclude[sel],
                       level = x$level, level.comb = x$level.comb,
                       hakn = x$hakn,
                       method.tau = x$method.tau,
                       method.tau.ci = x$method.tau.ci,
                       tau.preset = tau.preset, TE.tau = x$TE.tau)
    ##
    else if (gen)
      meta1 <- metagen(x$TE[sel], x$seTE[sel],
                       sm = x$sm,
                       studlab = x$studlab[sel],
                       exclude = x$exclude[sel],
                       level = x$level, level.comb = x$level.comb,
                       hakn = x$hakn,
                       method.tau = x$method.tau,
                       method.tau.ci = x$method.tau.ci,
                       tau.preset = tau.preset, TE.tau = x$TE.tau,
                       n.e = x$n.e[sel], n.c = x$n.c[sel],
                       warn = x$warn,
                       control = x$control)
    ##
    else if (inc)
      meta1 <- metainc(x$event.e[sel], x$time.e[sel],
                       x$event.c[sel], x$time.c[sel],
                       studlab = x$studlab[sel],
                       exclude = x$exclude[sel],
                       method = x$method,
                       sm = x$sm,
                       incr = if (length(x$incr) == 1) x$incr else x$incr[sel],
                       allincr = x$allincr,
                       addincr = x$addincr,
                       level = x$level, level.comb = x$level.comb,
                       hakn = x$hakn,
                       method.tau = x$method.tau,
                       method.tau.ci = x$method.tau.ci,
                       tau.preset = tau.preset,
                       TE.tau = x$TE.tau,
                       warn = x$warn,
                       control = x$control)
    ##
    else if (mean)
      meta1 <- metamean(x$n[sel], x$mean[sel], x$sd[sel],
                        sm = x$sm,
                        studlab = x$studlab[sel],
                        exclude = x$exclude[sel],
                        level = x$level, level.comb = x$level.comb,
                        hakn = x$hakn,
                        method.tau = x$method.tau,
                        method.tau.ci = x$method.tau.ci,
                        tau.preset = tau.preset, TE.tau = x$TE.tau,
                        warn = x$warn)
    ##
    else if (prop)
      meta1 <- metaprop(x$event[sel], x$n[sel],
                        sm = x$sm,
                        studlab = x$studlab[sel],
                        exclude = x$exclude[sel],
                        level = x$level, level.comb = x$level.comb,
                        incr = if (length(x$incr) == 1) x$incr else x$incr[sel],
                        allincr = x$allincr,
                        addincr = x$addincr,
                        hakn = x$hakn,
                        method = x$method,
                        method.tau = x$method.tau,
                        method.tau.ci = x$method.tau.ci,
                        tau.preset = tau.preset, TE.tau = x$TE.tau,
                        warn = x$warn)
    ##
    else if (rate)
      meta1 <- metarate(x$event[sel], x$time[sel],
                        sm = x$sm,
                        studlab = x$studlab[sel],
                        exclude = x$exclude[sel],
                        level = x$level, level.comb = x$level.comb,
                        incr = if (length(x$incr) == 1) x$incr else x$incr[sel],
                        allincr = x$allincr,
                        addincr = x$addincr,
                        hakn = x$hakn,
                        method.tau = x$method.tau,
                        method.tau.ci = x$method.tau.ci,
                        tau.preset = tau.preset, TE.tau = x$TE.tau,
                        warn = x$warn,
                        control = x$control)
    ##
    else
      stop("No meta-analysis object used for subgroup analysis.")
    ##
    ##
    res.w[j,] <- c(meta1$TE.fixed,                            #  1
                   meta1$seTE.fixed,                          #  2
                   meta1$Q,                                   #  3
                   meta1$k,                                   #  4
                   length(meta1$TE),                          #  5
                   meta1$TE.random,                           #  6
                   meta1$seTE.random,                         #  7
                   meta1$H,                                   #  8
                   meta1$lower.H,                             #  9
                   meta1$upper.H,                             # 10
                   meta1$I2,                                  # 11
                   meta1$lower.I2,                            # 12
                   meta1$upper.I2,                            # 13
                   meta1$tau2,                                # 14
                   meta1$lower.tau2,                          # 15
                   meta1$upper.tau2,                          # 16
                   meta1$tau,                                 # 17
                   meta1$lower.tau,                           # 18
                   meta1$upper.tau,                           # 19
                   1 / mean(1 / x$n[sel]),                    # 20
                   sum(x$w.fixed[sel]),                       # 21
                   sum(x$w.random[sel]),                      # 22
                   if (bin.inc) sumNA(meta1$event.e) else NA, # 23
                   if (bin.cont.gen) sumNA(meta1$n.e) else NA,# 24
                   if (bin.inc) sumNA(meta1$event.c) else NA, # 25
                   if (bin.cont.gen) sumNA(meta1$n.c) else NA,# 26
                   if (prop) sumNA(meta1$event) else NA,      # 27
                   if (cor.prop.mean) sumNA(meta1$n) else NA, # 28
                   if (inc) sumNA(meta1$time.e) else NA,      # 29
                   if (inc) sumNA(meta1$time.c) else NA,      # 30
                   1 / mean(1 / x$time[sel]),                 # 31
                   meta1$Rb,                                  # 32
                   meta1$lower.Rb,                            # 33
                   meta1$upper.Rb                             # 34
                   )
    ##
    add.w[j, ] <- c(meta1$sign.lower.tau, # 1
                    meta1$sign.upper.tau  # 2
                    )
  }
  ##
  TE.fixed.w    <- res.w[, 1]
  seTE.fixed.w  <- res.w[, 2]
  Q.w           <- res.w[, 3]
  k.w           <- res.w[, 4]
  k.all.w       <- res.w[, 5]
  TE.random.w   <- res.w[, 6]
  seTE.random.w <- res.w[, 7]
  ##
  H.w     <- res.w[, 8]
  H.w.low <- res.w[, 9]
  H.w.upp <- res.w[, 10]
  ##
  I2.w     <- res.w[, 11]
  I2.w.low <- res.w[, 12]
  I2.w.upp <- res.w[, 13]
  ##
  tau2.w <- res.w[, 14]
  lower.tau2.w <- res.w[, 15]
  upper.tau2.w <- res.w[, 16]
  tau.w <- res.w[, 17]
  lower.tau.w <- res.w[, 18]
  upper.tau.w <- res.w[, 19]
  ##
  sign.lower.tau.w <- add.w[, 1]
  sign.upper.tau.w <- add.w[, 2]
  ##
  n.harmonic.mean.w <- res.w[, 20]
  ##
  w.fixed.w  <- res.w[, 21]
  w.random.w <- res.w[, 22]
  ##
  event.e.w <- res.w[, 23]
  n.e.w     <- res.w[, 24]
  event.c.w <- res.w[, 25]
  n.c.w     <- res.w[, 26]
  event.w   <- res.w[, 27]
  n.w       <- res.w[, 28]
  ##
  time.e.w <- res.w[, 29]
  time.c.w <- res.w[, 30]
  t.harmonic.mean.w <- res.w[, 31]
  ##
  Rb.w     <- res.w[, 32]
  Rb.w.low <- res.w[, 33]
  Rb.w.upp <- res.w[, 34]
  ##
  ## GLMM with common tau-squared
  ##
  if (x$method == "GLMM" & !missing(byvar.glmm)) {
    mod <- as.call(~ byvar.glmm - 1)
    mod.Q <- as.call(~ byvar.glmm)
    ##
    if (prop) {
      glmm.fixed <- rma.glmm(xi = x$event, ni = x$n,
                             mods = mod,
                             method = "FE", test = ifelse(x$hakn, "t", "z"),
                             level = 100 * x$level.comb,
                             measure = "PLO", ...)
      ##
      glmm.random <- rma.glmm(xi = x$event, ni = x$n,
                              mods = mod,
                              method = x$method.tau,
                              test = ifelse(x$hakn, "t", "z"),
                              level = 100 * x$level.comb,
                              measure = "PLO", ...)
      ##
      ## Test for subgroup differences
      ##
      oldopts <- options(warn = -1)
      glmm.fixed.Q <- rma.glmm(xi = x$event, ni = x$n,
                               mods = mod.Q,
                               method = "FE", test = ifelse(x$hakn, "t", "z"),
                               level = 100 * x$level.comb,
                               measure = "PLO", ...)
      ##
      glmm.random.Q <- rma.glmm(xi = x$event, ni = x$n,
                                mods = mod.Q,
                                method = x$method.tau,
                                test = ifelse(x$hakn, "t", "z"),
                                level = 100 * x$level.comb,
                                measure = "PLO", ...)
      options(oldopts)
    }
    ##
    else if (rate) {
      glmm.fixed <- rma.glmm(xi = x$event, ti = x$time,
                             mods = mod,
                             method = "FE", test = ifelse(x$hakn, "t", "z"),
                             level = 100 * x$level.comb,
                             measure = "IRLN", control = x$control, ...)
      ##
      glmm.random <- rma.glmm(xi = x$event, ti = x$time,
                              mods = mod,
                              method = x$method.tau,
                              test = ifelse(x$hakn, "t", "z"),
                              level = 100 * x$level.comb,
                              measure = "IRLN", control = x$control, ...)
      ##
      ## Test for subgroup differences
      ##
      oldopts <- options(warn = -1)
      glmm.fixed.Q <- rma.glmm(xi = x$event, ti = x$time,
                               mods = mod.Q,
                               method = "FE", test = ifelse(x$hakn, "t", "z"),
                               level = 100 * x$level.comb,
                               measure = "IRLN", control = x$control, ...)
      ##
      glmm.random.Q <- rma.glmm(xi = x$event, ti = x$time,
                                mods = mod.Q,
                                method = x$method.tau,
                                test = ifelse(x$hakn, "t", "z"),
                                level = 100 * x$level.comb,
                                measure = "IRLN", control = x$control, ...)
      options(oldopts)
    }
    ##
    else if (inc) {
      glmm.fixed <- rma.glmm(x1i = x$event.e, t1i = x$time.e,
                             x2i = x$event.c, t2i = x$time.c,
                             mods = mod,
                             method = "FE", test = ifelse(x$hakn, "t", "z"),
                             model = x$model.glmm,
                             level = 100 * x$level.comb,
                             measure = "IRR", control = x$control, ...)
      ##
      glmm.random <- rma.glmm(x1i = x$event.e, t1i = x$time.e,
                              x2i = x$event.c, t2i = x$time.c,
                              mods = mod,
                              method = x$method.tau,
                              model = x$model.glmm,
                              test = ifelse(x$hakn, "t", "z"),
                              level = 100 * x$level.comb,
                              measure = "IRR", control = x$control, ...)
      ##
      ## Test for subgroup differences
      ##
      oldopts <- options(warn = -1)
      glmm.fixed.Q <- rma.glmm(x1i = x$event.e, t1i = x$time.e,
                               x2i = x$event.c, t2i = x$time.c,
                               mods = mod.Q,
                               method = "FE", test = ifelse(x$hakn, "t", "z"),
                               model = x$model.glmm,
                               level = 100 * x$level.comb,
                               measure = "IRR", control = x$control, ...)
      ##
      glmm.random.Q <- rma.glmm(x1i = x$event.e, t1i = x$time.e,
                                x2i = x$event.c, t2i = x$time.c,
                                mods = mod.Q,
                                method = x$method.tau,
                                model = x$model.glmm,
                                test = ifelse(x$hakn, "t", "z"),
                                level = 100 * x$level.comb,
                                measure = "IRR", control = x$control, ...)
      options(oldopts)
    }
    ##
    else if (bin) {
      glmm.fixed <- rma.glmm(ai = x$event.e, n1i = x$n.e,
                             ci = x$event.c, n2i = x$n.c,
                             mods = mod,
                             method = "FE", test = ifelse(x$hakn, "t", "z"),
                             model = x$model.glmm,
                             level = 100 * x$level.comb,
                             measure = "OR", control = x$control, ...)
      ##
      glmm.random <- rma.glmm(ai = x$event.e, n1i = x$n.e,
                              ci = x$event.c, n2i = x$n.c,
                              mods = mod,
                              method = x$method.tau,
                              model = x$model.glmm,
                              test = ifelse(x$hakn, "t", "z"),
                              level = 100 * x$level.comb,
                              measure = "OR", control = x$control, ...)
      ##
      ## Test for subgroup differences
      ##
      oldopts <- options(warn = -1)
      glmm.fixed.Q <- rma.glmm(ai = x$event.e, n1i = x$n.e,
                               ci = x$event.c, n2i = x$n.c,
                               mods = mod.Q,
                               method = "FE", test = ifelse(x$hakn, "t", "z"),
                               model = x$model.glmm,
                               level = 100 * x$level.comb,
                               measure = "OR", control = x$control, ...)
      ##
      glmm.random.Q <- rma.glmm(ai = x$event.e, n1i = x$n.e,
                                ci = x$event.c, n2i = x$n.c,
                                mods = mod.Q,
                                method = x$method.tau,
                                model = x$model.glmm,
                                test = ifelse(x$hakn, "t", "z"),
                                level = 100 * x$level.comb,
                                measure = "OR", control = x$control, ...)
      options(oldopts)
    }
    ##
    TE.fixed.w   <- as.numeric(glmm.fixed$b)
    seTE.fixed.w <- as.numeric(glmm.fixed$se)
    TE.random.w   <- as.numeric(glmm.random$b)
    seTE.random.w <- as.numeric(glmm.random$se)
    ##
    tau2.w <- rep_len(glmm.random$tau2, n.bylevs)
    lower.tau2.w <- upper.tau2.w <- rep_len(NA, n.bylevs)
    tau.w <- rep_len(sqrt(glmm.random$tau2), n.bylevs)
    lower.tau.w <- upper.tau.w <- rep_len(NA, n.bylevs)
    sign.lower.tau.w <- sign.upper.tau.w <- rep_len("", n.bylevs)
    ##
    Rb.w     <- rep_len(NA, n.bylevs)
    Rb.w.low <- rep_len(NA, n.bylevs)
    Rb.w.upp <- rep_len(NA, n.bylevs)
  }
  ##
  ci.fixed.w  <- ci(TE.fixed.w, seTE.fixed.w, x$level.comb)
  ##
  if (!is.null(x$hakn) && x$hakn)
    ci.random.w <- ci(TE.random.w, seTE.random.w, x$level.comb, df = k.w - 1)
  else
    ci.random.w <- ci(TE.random.w, seTE.random.w, x$level.comb)
  ##
  ## Tests for subgroup differences
  ##
  if (x$method == "GLMM" & !missing(byvar.glmm)) {
    Q.w.fixed <- glmm.fixed.Q$QE.Wld
    df.Q.w <- glmm.fixed.Q$k.eff - glmm.fixed.Q$p.eff
    pval.Q.w.fixed  <- glmm.fixed.Q$QEp.Wld
    ##
    Q.b.fixed  <- glmm.fixed.Q$QM
    Q.b.random <- glmm.random.Q$QM
    ##
    df.Q.b <- glmm.fixed.Q$p.eff - 1
    ##
    pval.Q.b.fixed  <- glmm.fixed.Q$QMp
    pval.Q.b.random <- glmm.random.Q$QMp
  }
  else {
    Q.w.fixed <- sum(Q.w, na.rm = TRUE)
    df.Q.w <- sum((k.w - 1)[!is.na(Q.w)])
    pval.Q.w.fixed  <- pvalQ(Q.w.fixed, df.Q.w)
    ##
    Q.b.fixed  <- metagen(TE.fixed.w, seTE.fixed.w)$Q
    Q.b.random <- metagen(TE.random.w, seTE.random.w)$Q
    ##
    df.Q.b <- ifelse(x$k == 0, 0, x$k - 1 - sum((k.w - 1)[!is.na(Q.w)]))
    ##
    pval.Q.b.fixed  <- pvalQ(Q.b.fixed, df.Q.b)
    pval.Q.b.random <- pvalQ(Q.b.random, df.Q.b)
  }
  
  
  res <- list(bylevs = bylevs,
              ##
              TE.fixed.w = ci.fixed.w$TE,
              seTE.fixed.w = ci.fixed.w$seTE,
              lower.fixed.w = ci.fixed.w$lower,
              upper.fixed.w = ci.fixed.w$upper,
              zval.fixed.w = ci.fixed.w$z,
              pval.fixed.w = ci.fixed.w$p,
              w.fixed.w = w.fixed.w,
              ##
              TE.random.w = ci.random.w$TE,
              seTE.random.w = ci.random.w$seTE,
              lower.random.w = ci.random.w$lower,
              upper.random.w = ci.random.w$upper,
              zval.random.w = ci.random.w$z,
              pval.random.w = ci.random.w$p,
              df.hakn.w = ci.random.w$df,
              w.random.w = w.random.w,
              ##
              n.harmonic.mean.w = n.harmonic.mean.w,
              t.harmonic.mean.w = t.harmonic.mean.w,
              ##
              event.e.w = event.e.w,
              time.e.w = time.e.w,
              n.e.w = n.e.w,
              event.c.w = event.c.w,
              time.c.w = time.c.w,
              n.c.w = n.c.w,
              n.w = n.w,
              event.w = event.w,
              ##
              k.w = k.w,
              k.all.w = k.all.w,
              Q.w = Q.w,
              pval.Q.w = pvalQ(Q.w, k.w - 1),
              ##
              tau2.w = tau2.w,
              lower.tau2.w = lower.tau2.w,
              upper.tau2.w = upper.tau2.w,
              tau.w = tau.w,
              lower.tau.w = lower.tau.w,
              upper.tau.w = upper.tau.w,
              sign.lower.tau.w = sign.lower.tau.w,
              sign.upper.tau.w = sign.upper.tau.w,
              ##
              H.w = H.w,
              lower.H.w = H.w.low,
              upper.H.w = H.w.upp,
              ##
              I2.w = I2.w,
              lower.I2.w = I2.w.low,
              upper.I2.w = I2.w.upp,
              ##
              Rb.w = Rb.w,
              lower.Rb.w = Rb.w.low,
              upper.Rb.w = Rb.w.upp,
              ##
              Q.w.fixed = Q.w.fixed,
              Q.w.random = NA,
              df.Q.w = df.Q.w,
              pval.Q.w.fixed = pval.Q.w.fixed,
              pval.Q.w.random = NA,
              ##
              Q.b.fixed = Q.b.fixed,
              Q.b.random = Q.b.random,
              df.Q.b = df.Q.b,
              pval.Q.b.fixed = pval.Q.b.fixed,
              pval.Q.b.random = pval.Q.b.random
              )
  
  
  res
}
