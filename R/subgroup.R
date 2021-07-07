subgroup <- function(x, tau.preset = NULL, byvar.rma, ...) {
  
  
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
  
  
  res.w <- matrix(NA, ncol = 48, nrow = n.bylevs)
  add.w <- matrix("", ncol =  2, nrow = n.bylevs)
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
                       adhoc.hakn = x$adhoc.hakn,
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
                        id = if (!is.null(x$id)) x$id[sel] else NULL,
                        sm = x$sm, pooledvar = x$pooledvar,
                        level = x$level, level.comb = x$level.comb,
                        hakn = x$hakn,
                        adhoc.hakn = x$adhoc.hakn,
                        method.tau = x$method.tau,
                        method.tau.ci = x$method.tau.ci,
                        tau.preset = tau.preset, TE.tau = x$TE.tau,
                        warn = x$warn,
                        control = x$control)
    ##
    else if (cor)
      meta1 <- metacor(x$cor[sel], x$n[sel],
                       sm = x$sm,
                       studlab = x$studlab[sel],
                       exclude = x$exclude[sel],
                       level = x$level, level.comb = x$level.comb,
                       hakn = x$hakn,
                       adhoc.hakn = x$adhoc.hakn,
                       method.tau = x$method.tau,
                       method.tau.ci = x$method.tau.ci,
                       tau.preset = tau.preset, TE.tau = x$TE.tau,
                       null.effect = x$null.effect,
                       control = x$control)
    ##
    else if (gen)
      meta1 <- metagen(x$TE[sel], x$seTE[sel],
                       sm = x$sm,
                       studlab = x$studlab[sel],
                       exclude = x$exclude[sel],
                       id = if (!is.null(x$id)) x$id[sel] else NULL,
                       level = x$level, level.comb = x$level.comb,
                       hakn = x$hakn,
                       adhoc.hakn = x$adhoc.hakn,
                       method.tau = x$method.tau,
                       method.tau.ci = x$method.tau.ci,
                       tau.preset = tau.preset, TE.tau = x$TE.tau,
                       null.effect = x$null.effect,
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
                       adhoc.hakn = x$adhoc.hakn,
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
                        adhoc.hakn = x$adhoc.hakn,
                        method.tau = x$method.tau,
                        method.tau.ci = x$method.tau.ci,
                        tau.preset = tau.preset, TE.tau = x$TE.tau,
                        null.effect = x$null.effect,
                        warn = x$warn,
                        control = x$control)
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
                        adhoc.hakn = x$adhoc.hakn,
                        method = x$method,
                        method.tau = x$method.tau,
                        method.tau.ci = x$method.tau.ci,
                        tau.preset = tau.preset, TE.tau = x$TE.tau,
                        null.effect = x$null.effect,
                        warn = x$warn,
                        control = x$control)
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
                        adhoc.hakn = x$adhoc.hakn,
                        method.tau = x$method.tau,
                        method.tau.ci = x$method.tau.ci,
                        tau.preset = tau.preset, TE.tau = x$TE.tau,
                        null.effect = x$null.effect,
                        warn = x$warn,
                        control = x$control)
    ##
    else
      stop("No meta-analysis object used for subgroup analysis.")
    ##
    n.tau <- length(meta1$tau)
    n.tau.ci <- length(meta1$lower.tau)
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
                   meta1$tau2,                                # 14-15
                   if (n.tau == 1) NA,                        #
                   meta1$lower.tau2,                          # 16-17
                   if (n.tau.ci == 1) NA,                     #
                   meta1$upper.tau2,                          # 18-19
                   if (n.tau.ci == 1) NA,                     #
                   meta1$tau,                                 # 20-21
                   if (n.tau == 1) NA,                        #
                   meta1$lower.tau,                           # 22-23
                   if (n.tau.ci == 1) NA,                     #
                   meta1$upper.tau,                           # 24-25
                   if (n.tau.ci == 1) NA,                     #
                   1 / mean(1 / x$n[sel]),                    # 26
                   sum(x$w.fixed[sel]),                       # 27
                   sum(x$w.random[sel]),                      # 28
                   if (bin.inc) sumNA(meta1$event.e) else NA, # 29
                   if (bin.cont.gen) sumNA(meta1$n.e) else NA,# 30
                   if (bin.inc) sumNA(meta1$event.c) else NA, # 31
                   if (bin.cont.gen) sumNA(meta1$n.c) else NA,# 32
                   if (prop) sumNA(meta1$event) else NA,      # 33
                   if (cor.prop.mean) sumNA(meta1$n) else NA, # 34
                   if (inc) sumNA(meta1$time.e) else NA,      # 35
                   if (inc) sumNA(meta1$time.c) else NA,      # 36
                   1 / mean(1 / x$time[sel]),                 # 37
                   meta1$Rb,                                  # 38
                   meta1$lower.Rb,                            # 39
                   meta1$upper.Rb,                            # 40
                   meta1$lower.fixed,                         # 41
                   meta1$upper.fixed,                         # 42
                   meta1$statistic.fixed,                     # 41
                   meta1$pval.fixed,                          # 42
                   meta1$lower.random,                        # 45
                   meta1$upper.random,                        # 46
                   meta1$statistic.random,                    # 47
                   meta1$pval.random                          # 48
                   )
    ##
    if (n.tau.ci == 1)
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
  tau2.1.w <- res.w[, 14]
  tau2.2.w <- res.w[, 15]
  lower.tau2.1.w <- res.w[, 16]
  lower.tau2.2.w <- res.w[, 17]
  upper.tau2.1.w <- res.w[, 18]
  upper.tau2.2.w <- res.w[, 19]
  ##
  tau.1.w <- res.w[, 20]
  tau.2.w <- res.w[, 21]
  lower.tau.1.w <- res.w[, 22]
  lower.tau.2.w <- res.w[, 23]
  upper.tau.1.w <- res.w[, 24]
  upper.tau.2.w <- res.w[, 25]
  ##
  tau2.w <- ifelse(is.na(tau2.2.w), tau2.1.w, tau2.1.w + tau2.2.w)
  lower.tau2.w <- ifelse(is.na(tau2.2.w), lower.tau2.1.w, NA)
  upper.tau2.w <- ifelse(is.na(tau2.2.w), upper.tau2.1.w, NA)
  ##
  tau.w <- ifelse(is.na(tau2.2.w), tau.1.w, sqrt(tau.1.w^2 + tau.2.w^2))
  lower.tau.w <- ifelse(is.na(tau2.2.w), lower.tau.1.w, NA)
  upper.tau.w <- ifelse(is.na(tau2.2.w), upper.tau.1.w, NA)
  ##
  sign.lower.tau.w <- add.w[, 1]
  sign.upper.tau.w <- add.w[, 2]
  ##
  n.harmonic.mean.w <- res.w[, 26]
  ##
  w.fixed.w  <- res.w[, 27]
  w.random.w <- res.w[, 28]
  ##
  event.e.w <- res.w[, 29]
  n.e.w     <- res.w[, 30]
  event.c.w <- res.w[, 31]
  n.c.w     <- res.w[, 32]
  event.w   <- res.w[, 33]
  n.w       <- res.w[, 34]
  ##
  time.e.w <- res.w[, 35]
  time.c.w <- res.w[, 36]
  t.harmonic.mean.w <- res.w[, 37]
  ##
  Rb.w     <- res.w[, 38]
  Rb.w.low <- res.w[, 39]
  Rb.w.upp <- res.w[, 40]
  ##
  lower.fixed.w <- res.w[, 41]
  upper.fixed.w <- res.w[, 42]
  statistic.fixed.w <- res.w[, 43]
  pval.fixed.w <- res.w[, 44]
  ##
  lower.random.w <- res.w[, 45]
  upper.random.w <- res.w[, 46]
  statistic.random.w <- res.w[, 47]
  pval.random.w <- res.w[, 48]
  ##
  ## Three-level model with common tau-squared
  ##
  if (!is.null(x$three.level) && x$three.level && !missing(byvar.rma)) {
    mod <- as.call(~ byvar.rma - 1)
    mod.Q <- as.call(~ byvar.rma)
    ##
    no.with.id <- x$id
    runID <- seq_along(no.with.id)
    ##
    mv.random <- rma.mv(x$TE, x$seTE^2,
                        mods = mod,
                        random = ~ 1 | no.with.id / runID,
                        method = x$method.tau,
                        test = ifelse(x$hakn, "t", "z"),
                        level = 100 * x$level.comb)
    ##
    mv.random.Q <-
      suppressWarnings(rma.mv(x$TE, x$seTE^2,
                              mods = mod.Q,
                              random = ~ 1 | no.with.id / runID,
                              method = x$method.tau,
                              test = ifelse(x$hakn, "t", "z"),
                              level = 100 * x$level.comb))
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
    tau2.w <- rep_len(sum(mv.random$sigma2), n.bylevs)
    lower.tau2.w <- upper.tau2.w <- rep_len(NA, n.bylevs)
    tau.w <- rep_len(sqrt(sum(mv.random$sigma2)), n.bylevs)
    lower.tau.w <- upper.tau.w <- rep_len(NA, n.bylevs)
    sign.lower.tau.w <- sign.upper.tau.w <- rep_len("", n.bylevs)
    ##
    tau2.1.w <- tau2.w
    tau2.2.w <- rep_len(NA, n.bylevs)
    lower.tau2.1.w <- lower.tau2.w
    lower.tau2.2.w <- rep_len(NA, n.bylevs)
    upper.tau2.1.w <- upper.tau2.w
    upper.tau2.2.w <- rep_len(NA, n.bylevs)
    ##
    Rb.w     <- rep_len(NA, n.bylevs)
    Rb.w.low <- rep_len(NA, n.bylevs)
    Rb.w.upp <- rep_len(NA, n.bylevs)
    ##
    ci.fixed.w  <- ci(TE.fixed.w, seTE.fixed.w, x$level.comb)
    ##
    if (!is.null(x$hakn) && x$hakn)
      ci.random.w <- ci(TE.random.w, seTE.random.w, x$level.comb, df = k.w - 1)
    else
      ci.random.w <- ci(TE.random.w, seTE.random.w, x$level.comb)
    ##
    lower.fixed.w <- ci.fixed.w$lower
    upper.fixed.w <- ci.fixed.w$upper
    statistic.fixed.w <- ci.fixed.w$statistic
    pval.fixed.w <- ci.fixed.w$p
    ##
    lower.random.w <- ci.random.w$lower
    upper.random.w <- ci.random.w$upper
    statistic.random.w <- ci.random.w$statistic
    pval.random.w <- ci.random.w$p
  }
  ##
  ## GLMM with common tau-squared
  ##
  if (x$method == "GLMM" & !missing(byvar.rma)) {
    mod <- as.call(~ byvar.rma - 1)
    mod.Q <- as.call(~ byvar.rma)
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
      glmm.fixed.Q <-
        suppressWarnings(rma.glmm(xi = x$event, ni = x$n,
                                  mods = mod.Q,
                                  method = "FE",
                                  test = ifelse(x$hakn, "t", "z"),
                                  level = 100 * x$level.comb,
                                  measure = "PLO", ...))
      ##
      glmm.random.Q <-
        suppressWarnings(rma.glmm(xi = x$event, ni = x$n,
                                  mods = mod.Q,
                                  method = x$method.tau,
                                  test = ifelse(x$hakn, "t", "z"),
                                  level = 100 * x$level.comb,
                                  measure = "PLO", ...))
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
      glmm.fixed.Q <-
        suppressWarnings(rma.glmm(xi = x$event, ti = x$time,
                                  mods = mod.Q,
                                  method = "FE",
                                  test = ifelse(x$hakn, "t", "z"),
                                  level = 100 * x$level.comb,
                                  measure = "IRLN", control = x$control, ...))
      ##
      glmm.random.Q <-
        suppressWarnings(rma.glmm(xi = x$event, ti = x$time,
                                  mods = mod.Q,
                                  method = x$method.tau,
                                  test = ifelse(x$hakn, "t", "z"),
                                  level = 100 * x$level.comb,
                                  measure = "IRLN", control = x$control, ...))
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
      glmm.fixed.Q <-
        suppressWarnings(rma.glmm(x1i = x$event.e, t1i = x$time.e,
                                  x2i = x$event.c, t2i = x$time.c,
                                  mods = mod.Q,
                                  method = "FE",
                                  test = ifelse(x$hakn, "t", "z"),
                                  model = x$model.glmm,
                                  level = 100 * x$level.comb,
                                  measure = "IRR", control = x$control, ...))
      ##
      glmm.random.Q <-
        suppressWarnings(rma.glmm(x1i = x$event.e, t1i = x$time.e,
                                  x2i = x$event.c, t2i = x$time.c,
                                  mods = mod.Q,
                                  method = x$method.tau,
                                  model = x$model.glmm,
                                  test = ifelse(x$hakn, "t", "z"),
                                  level = 100 * x$level.comb,
                                  measure = "IRR", control = x$control, ...))
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
      glmm.fixed.Q <-
        suppressWarnings(rma.glmm(ai = x$event.e, n1i = x$n.e,
                                  ci = x$event.c, n2i = x$n.c,
                                  mods = mod.Q,
                                  method = "FE",
                                  test = ifelse(x$hakn, "t", "z"),
                                  model = x$model.glmm,
                                  level = 100 * x$level.comb,
                                  measure = "OR", control = x$control, ...))
      ##
      glmm.random.Q <-
        suppressWarnings(rma.glmm(ai = x$event.e, n1i = x$n.e,
                                  ci = x$event.c, n2i = x$n.c,
                                  mods = mod.Q,
                                  method = x$method.tau,
                                  model = x$model.glmm,
                                  test = ifelse(x$hakn, "t", "z"),
                                  level = 100 * x$level.comb,
                                  measure = "OR", control = x$control, ...))
    }
    ##
    if (length(TE.fixed.w) != length(as.numeric(glmm.fixed$b))) {
      TE.fixed.w[!is.na(TE.fixed.w)] <- as.numeric(glmm.fixed$b)
      seTE.fixed.w[!is.na(seTE.fixed.w)] <- as.numeric(glmm.fixed$se)
      TE.random.w[!is.na(TE.random.w)] <- as.numeric(glmm.random$b)
      seTE.random.w[!is.na(seTE.random.w)] <- as.numeric(glmm.random$se)
    }
    else {
      TE.fixed.w   <- as.numeric(glmm.fixed$b)
      seTE.fixed.w <- as.numeric(glmm.fixed$se)
      TE.random.w   <- as.numeric(glmm.random$b)
      seTE.random.w <- as.numeric(glmm.random$se)
    }
    ##
    tau2.w <- rep_len(glmm.random$tau2, n.bylevs)
    lower.tau2.w <- upper.tau2.w <- rep_len(NA, n.bylevs)
    tau.w <- rep_len(sqrt(glmm.random$tau2), n.bylevs)
    lower.tau.w <- upper.tau.w <- rep_len(NA, n.bylevs)
    sign.lower.tau.w <- sign.upper.tau.w <- rep_len("", n.bylevs)
    ##
    tau2.1.w <- tau2.w
    tau2.2.w <- rep_len(NA, n.bylevs)
    lower.tau2.1.w <- lower.tau2.w
    lower.tau2.2.w <- rep_len(NA, n.bylevs)
    upper.tau2.1.w <- upper.tau2.w
    upper.tau2.2.w <- rep_len(NA, n.bylevs)
    ##
    Rb.w     <- rep_len(NA, n.bylevs)
    Rb.w.low <- rep_len(NA, n.bylevs)
    Rb.w.upp <- rep_len(NA, n.bylevs)
    ##
    ci.fixed.w  <- ci(TE.fixed.w, seTE.fixed.w, x$level.comb)
    ##
    if (!is.null(x$hakn) && x$hakn)
      ci.random.w <- ci(TE.random.w, seTE.random.w, x$level.comb, df = k.w - 1)
    else
      ci.random.w <- ci(TE.random.w, seTE.random.w, x$level.comb)
    ##
    lower.fixed.w <- ci.fixed.w$lower
    upper.fixed.w <- ci.fixed.w$upper
    statistic.fixed.w <- ci.fixed.w$statistic
    pval.fixed.w <- ci.fixed.w$p
    ##
    lower.random.w <- ci.random.w$lower
    upper.random.w <- ci.random.w$upper
    statistic.random.w <- ci.random.w$statistic
    pval.random.w <- ci.random.w$p
  }
  ##
  ## Tests for subgroup differences
  ##
  if (x$method == "GLMM" & !missing(byvar.rma)) {
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
  else if (!is.null(x$three.level) && x$three.level && !missing(byvar.rma)) {
    Q.w.fixed <- NA
    df.Q.w <- mv.random.Q$k.eff - mv.random.Q$p.eff
    pval.Q.w.fixed <- NA
    ##
    Q.b.fixed  <- NA
    Q.b.random <- mv.random.Q$QM
    ##
    df.Q.b <- mv.random.Q$p.eff - 1
    ##
    pval.Q.b.fixed  <- NA
    pval.Q.b.random <- mv.random.Q$QMp
  }
  else {
    Q.w.fixed <- sum(Q.w, na.rm = TRUE)
    df.Q.w <- sum((k.w - 1)[!is.na(Q.w)])
    pval.Q.w.fixed  <- pvalQ(Q.w.fixed, df.Q.w)
    ##
    Q.b.fixed  <- metagen(TE.fixed.w, seTE.fixed.w, method.tau = "DL")$Q
    Q.b.random <- metagen(TE.random.w, seTE.random.w, method.tau = "DL")$Q
    ##
    df.Q.b <- ifelse(x$k == 0, 0, x$k - 1 - sum((k.w - 1)[!is.na(Q.w)]))
    ##
    pval.Q.b.fixed  <- pvalQ(Q.b.fixed, df.Q.b)
    pval.Q.b.random <- pvalQ(Q.b.random, df.Q.b)
  }
  ##
  ## Prediction interval
  ##
  seTE.predict.w <- sqrt(seTE.random.w^2 + tau2.w)
  ci.p.w <- ci(TE.random.w, seTE.predict.w, x$level.predict, k.w - 2)
  ##
  p.lower.w <- ci.p.w$lower
  p.upper.w <- ci.p.w$upper
  ##
  p.lower.w[k.w < 3] <- NA
  p.upper.w[k.w < 3] <- NA
  ##
  ## Degrees of freedom of Hartung-Knapp method
  ##
  df.hakn.w <- k.w - 1
  if (!x$hakn)
    df.hakn.w[!is.na(df.hakn.w)] <- NA
  
  
  res <- list(bylevs = bylevs,
              ##
              TE.fixed.w = TE.fixed.w,
              seTE.fixed.w = seTE.fixed.w,
              lower.fixed.w = lower.fixed.w,
              upper.fixed.w = upper.fixed.w,
              statistic.fixed.w = statistic.fixed.w,
              zval.fixed.w = statistic.fixed.w,
              pval.fixed.w = pval.fixed.w,
              w.fixed.w = w.fixed.w,
              ##
              TE.random.w = TE.random.w,
              seTE.random.w = seTE.random.w,
              lower.random.w = lower.random.w,
              upper.random.w = upper.random.w,
              statistic.random.w = statistic.random.w,
              zval.random.w = statistic.random.w,
              pval.random.w = pval.random.w,
              w.random.w = w.random.w,
              ##
              seTE.predict.w = seTE.predict.w,
              lower.predict.w = p.lower.w, upper.predict.w = p.upper.w,
              ##
              df.hakn.w = df.hakn.w,
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
              tau2.1.w = tau2.1.w,
              lower.tau2.1.w = lower.tau2.1.w,
              upper.tau2.1.w = upper.tau2.1.w,
              tau.1.w = tau.1.w,
              lower.tau.1.w = lower.tau.1.w,
              upper.tau.1.w = upper.tau.1.w,
              ##
              tau2.2.w = tau2.2.w,
              lower.tau2.2.w = lower.tau2.2.w,
              upper.tau2.2.w = upper.tau2.2.w,
              tau.2.w = tau.2.w,
              lower.tau.2.w = lower.tau.2.w,
              upper.tau.2.w = upper.tau.2.w,
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
