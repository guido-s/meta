meta2meth <- function(x, outclab = NULL) {
  list(sm = x$sm,
       method = x$method,
       method.random = x$method.random,
       three.level = replaceNULL(x$three.level, NA),
       ##
       k = x$k,
       k.study = x$k.study,
       k.all = x$k.all,
       k.TE = x$k.TE,
       ##
       level = x$level.ma,
       level.ma = x$level.ma,
       level.predict = x$level.predict,
       ##
       common = x$common,
       random = x$random,
       overall = x$overall,
       overall.hetstat = x$overall.hetstat,
       ##
       method.common.ci = x$method.common.ci,
       method.random.ci = x$method.random.ci,
       adhoc.hakn.ci = x$adhoc.hakn.ci,
       method.tau = x$method.tau,
       tau.preset = replaceNULL(x$tau.preset),
       TE.tau = replaceNULL(x$TE.tau),
       tau.common = replaceNULL(x$tau.common, FALSE),
       #
       method.I2 = replaceNULL(x$method.I2, "Q"),
       #
       prediction = x$prediction,
       prediction.subgroup =
         replaceNULL(x$prediction.subgroup, FALSE),
       method.predict = x$method.predict,
       adhoc.hakn.pi = x$adhoc.hakn.pi,
       ##
       method.bias = "",
       null.effect = x$null.effect,
       ##
       title = x$title,
       complab = x$complab,
       outclab = replaceNULL(outclab, x$outclab),
       ##
       label.e = x$label.e,
       label.c = x$label.c,
       label.left = x$label.left,
       label.right = x$label.right,
       ##
       print.subgroup.name = FALSE,
       sep.subgroup = "",
       warn = replaceNULL(x$warn, FALSE),
       ##
       backtransf = x$backtransf,
       pscale = replaceNULL(x$pscale, 1),
       irscale = replaceNULL(x$irscale, 1),
       irunit = replaceNULL(x$ir.unit))
}


overall2meta <- function(x, common, random, prediction, name) {
  res <-
    list(
      TE.common = if (common) x$TE.common else NULL,
      seTE.common = if (common) x$seTE.common else NULL,
      lower.common = if (common) x$lower.common else NULL,
      upper.common = if (common) x$upper.common else NULL,
      statistic.common = if (common) x$statistic.common else NULL,
      pval.common = if (common) x$pval.common else NULL,
      text.common = if (common) x$text.common else NULL,
      ##
      TE.random = if (random) x$TE.random else NULL,
      seTE.random = if (random) x$seTE.random else NULL,
      lower.random = if (random) x$lower.random else NULL,
      upper.random = if (random) x$upper.random else NULL,
      statistic.random = if (random) x$statistic.random else NULL,
      pval.random = if (random) x$pval.random else NULL,
      text.random = if (random) x$text.random else NULL,
      ##
      df.random = if (random) replaceNULL(x$df.random) else NULL,
      df.hakn = if (random) replaceNULL(x$df.hakn) else NULL,
      df.kero = if (random) replaceNULL(x$df.kero) else NULL,
      ##
      n.harmonic.mean.ma = 1 / mean(1 / replaceNULL(x$n)),
      t.harmonic.mean.ma = 1 / mean(1 / replaceNULL(x$time)),
      ##
      seTE.predict = if (prediction) x$seTE.predict else NULL,
      df.predict = if (prediction) x$df.predict else NULL,
      lower.predict = if (prediction) x$lower.predict else NULL,
      upper.predict = if (prediction) x$upper.predict else NULL,
      ##
      n.e = sum(replaceNULL(x$n.e)),
      n.c = sum(replaceNULL(x$n.c)),
      ##
      Q = x$Q[1],
      df.Q = x$df.Q[1],
      pval.Q = x$pval.Q[1],
      ##
      tau2 = x$tau2,
      lower.tau2 = x$lower.tau2,
      upper.tau2 = x$upper.tau2,
      se.tau2 = replaceNULL(x$se.tau2),
      tau = x$tau,
      lower.tau = x$lower.tau,
      upper.tau = x$upper.tau,
      ##
      H = x$H,
      lower.H = x$lower.H,
      upper.H = x$upper.H,
      I2 = x$I2,
      lower.I2 = x$lower.I2,
      upper.I2 = x$upper.I2,
      Rb = x$Rb,
      lower.Rb = x$lower.Rb,
      upper.Rb = x$upper.Rb,
      ##
      Q.b.common = NA,
      Q.b.random = NA,
      df.Q.b = NA,
      pval.Q.b.common = NA,
      pval.Q.b.random = NA
      )
  ##
  res
}


subgr2meta <- function(x, common, random, prediction, name) {
  n.subgr <- length(x$subgroup.levels)
  ##
  res <- list(
    studlab =
      c(if (common)
          paste0(x$subgroup.levels,
                if (random | prediction) " (common effect)"),
        if (random)
          paste0(x$subgroup.levels,
                if (common | prediction) " (random effects)"),
        if (prediction)
          paste0(x$subgroup.levels,
                if (common | random) " (prediction)")
        ),
    ##
    subgroup = NA,
    ##
    subgroup.levels = name,
    ##
    TE =
      c(if (common) x$TE.common.w,
        if (random) x$TE.random.w,
        if (prediction) rep(NA, n.subgr)),
    ##
    seTE =
      c(if (common) x$seTE.common.w,
        if (random) x$seTE.random.w,
        if (prediction) rep(NA, n.subgr)),
    lower =
      c(if (common) x$lower.common.w,
        if (random) x$lower.random.w,
        if (prediction) x$lower.predict.w),
    upper =
      c(if (common) x$upper.common.w,
        if (random) x$upper.random.w,
        if (prediction) x$upper.predict.w),
    statistic =
      c(if (common) x$statistic.common.w,
        if (random) x$statistic.random.w,
        if (prediction) rep(NA, n.subgr)),
    pval =
      c(if (common) x$pval.common.w,
        if (random) x$pval.random.w,
        if (prediction) rep(NA, n.subgr)),
    ##
    TE.common = x$TE.common,
    seTE.common = x$seTE.common,
    lower.common = x$lower.common,
    upper.common = x$upper.common,
    statistic.common = x$statistic.common,
    pval.common = x$pval.common,
    text.common = if (common) x$text.common else NULL,
    ##
    TE.random = x$TE.random,
    seTE.random = x$seTE.random,
    lower.random = x$lower.random,
    upper.random = x$upper.random,
    statistic.random = x$statistic.random,
    pval.random = x$pval.random,
    text.random = if (random) x$text.random else NULL,
    ##
    lower.predict = x$lower.predict,
    upper.predict = x$upper.predict,
    text.predict = if (prediction) x$text.predict else NULL,
    ##
    df.random = replaceNULL(x$df.random),
    df.hakn = replaceNULL(x$df.hakn),
    df.kero = replaceNULL(x$df.kero),
    ##
    n.harmonic.mean = replaceNULL(x$n.harmonic.mean.w),
    t.harmonic.mean = replaceNULL(x$t.harmonic.mean.w),
    ##
    n.e = replaceNULL(x$n.e),
    n.c = replaceNULL(x$n.c),
    ##
    Q = x$Q[1],
    df.Q = x$df.Q[1],
    pval.Q = x$pval.Q[1],
    ##
    tau2 = x$tau2,
    tau = x$tau,
    H = x$H,
    lower.H = x$lower.H,
    upper.H = x$upper.H,
    I2 = x$I2,
    lower.I2 = x$lower.I2,
    upper.I2 = x$upper.I2,
    Rb = x$Rb,
    lower.Rb = x$lower.Rb,
    upper.Rb = x$upper.Rb,
    ##
    Q.b.common = NA,
    Q.b.random = NA,
    df.Q.b = NA,
    pval.Q.b.common = NA,
    pval.Q.b.random = NA)
  ##
  res$subgroup <- rep_len(name, length(res$TE))
  ##
  res
}


overall2subgr <- function(x) {
  res <- list()
  ##
  vars.c <- c("TE.common", "seTE.common", "lower.common", "upper.common",
              "statistic.common", "pval.common")
  ##
  for (i in vars.c)
    res[[paste0(i, ".w")]] <- replaceNULL(x[[i]])
  ##
  res[["w.common.w"]] <- rep(0, length(res[["lower.common.w"]]))
  ##
  vars.r <- c("TE.random", "seTE.random", "lower.random", "upper.random",
              "statistic.random", "pval.random",
              "df.random", "df.hakn", "df.kero")
  ##
  for (i in vars.r)
    res[[paste0(i, ".w")]] <- replaceNULL(x[[i]])
  ##
  res[["w.random.w"]] <- rep(0, length(res[["lower.random.w"]]))
  ##
  res[["n.harmonic.mean.ma"]] <- 1 / mean(1 / replaceNULL(x$n))
  res[["t.harmonic.mean.ma"]] <- 1 / mean(1 / replaceNULL(x$time))
  #
  res[["n.harmonic.mean.w"]] <- 1 / mean(1 / replaceNULL(x$n))
  res[["t.harmonic.mean.w"]] <- 1 / mean(1 / replaceNULL(x$time))
  ##
  res[["n.e.w"]] <- sum(replaceNULL(x$n.e))
  res[["n.c.w"]] <- sum(replaceNULL(x$n.c))
  ##
  vars <- c("k", "k.study", "k.all", "k.TE",
            "Q", "df.Q", "pval.Q",
            "tau2", "tau",
            "H", "lower.H", "upper.H",
            "I2", "lower.I2", "upper.I2",
            "Rb", "lower.Rb", "upper.Rb")
  ##
  for (i in vars)
    res[[paste0(i, ".w")]] <- replaceNULL(x[[i]])
  #
  res$Q.w <- res$Q.w[1]
  res$df.Q.w <- res$df.Q.w[1]
  res$pval.Q.w <- res$pval.Q.w[1]
  #
  res
}


subgr2data <- function(x, common, random, prediction, name, debug = FALSE) {
  n.subgr <- length(x$subgroup.levels)
  ##
  NAs <- rep_len(NA, n.subgr - 1)
  fillNAs <- rep_len(NA, n.subgr * (common + random + prediction - 1))
  ##
  res <- list()
  ##
  res$name <-
    c(if (common)
        paste(rep(name, n.subgr),
              if (random | prediction)
                "(common effect)"),
      if (random)
        paste(rep(name, n.subgr),
              if (common | prediction)
                "(random effects)", n.subgr),
      if (prediction)
        paste(rep(name, n.subgr),
              if (common | random)
                "(prediction)"))
  ##
  res$studlab <-
      c(if (common)
          paste0(x$subgroup.levels,
                if (random | prediction) " (common effect)"),
        if (random)
          paste0(x$subgroup.levels,
                if (common | prediction) " (random effects)"),
        if (prediction)
          paste0(x$subgroup.levels,
                if (common | random) " (prediction)"))
  ##
  res$subgroup <-
    c(if (common) rep(name, n.subgr),
      if (random) rep(name, n.subgr),
      if (prediction) rep(name, n.subgr))
  ##
  res$type <-
    c(if (common) rep("square", n.subgr),
      if (random) rep("square", n.subgr),
      if (prediction) rep("predict", n.subgr))
  ##
  res$model <-
    c(if (common) rep("common", n.subgr),
      if (random) rep("random", n.subgr),
      if (prediction) rep("predict", n.subgr))
  ##
  res$method <-
    c(if (common) rep(x$method, n.subgr),
      if (random) rep(x$method.random, n.subgr),
      if (prediction) rep(x$method.predict, n.subgr))
  ##
  res$method.tau <-
    c(if (common) rep("", n.subgr),
      if (random) rep(x$method.tau, n.subgr),
      if (prediction) rep("", n.subgr))
  ##
  res$three.level <-
    c(if (common) rep(FALSE, n.subgr),
      if (random) rep(x$three.level, n.subgr),
      if (prediction) rep(FALSE, n.subgr))
  ##
  res$method.tau.ci <-
    c(if (common) rep("", n.subgr),
      if (random) rep(x$method.tau.ci, n.subgr),
      if (prediction) rep("", n.subgr))
  ##
  res$tau.preset <-
    c(if (common) rep(NA, n.subgr),
      if (random) rep(replaceNULL(x$tau.preset, NA), n.subgr),
      if (prediction) rep(NA, n.subgr))
  #
  res$method.I2 <-
    c(if (common) rep(x$method.I2, n.subgr),
      if (random) rep(x$method.I2, n.subgr),
      if (prediction) rep("", n.subgr))
  #
  res$method.common.ci <-
    c(if (common) rep(x$method.common.ci, n.subgr),
      if (random) rep("", n.subgr),
      if (prediction) rep("", n.subgr))
  #
  res$method.random.ci <-
    c(if (common) rep("", n.subgr),
      if (random) rep(x$method.random.ci, n.subgr),
      if (prediction) rep("", n.subgr))
  ##
  res$df.random <-
    c(if (common) rep(NA, n.subgr),
      if (random) rep(x$df.random, n.subgr),
      if (prediction) rep(NA, n.subgr))
  ##
  res$adhoc.hakn.ci <-
    c(if (common) rep("", n.subgr),
      if (random) rep(x$adhoc.hakn.ci, n.subgr),
      if (prediction) rep("", n.subgr))
  ##
  res$method.predict <-
    c(if (common) rep("", n.subgr),
      if (random) rep("", n.subgr),
      if (prediction) rep(x$method.predict, n.subgr))
  ##
  res$df.predict <-
    c(if (common) rep(NA, n.subgr),
      if (random) rep(NA, n.subgr),
      if (prediction) rep(x$df.predict, n.subgr))
  ##
  res$adhoc.hakn.pi <-
    c(if (common) rep("", n.subgr),
      if (random) rep("", n.subgr),
      if (prediction) rep(x$adhoc.hakn.pi, n.subgr))
  ##
  res$rho <-
    c(if (common) rep(NA, n.subgr),
      if (random) rep(replaceNULL(x$rho), n.subgr),
      if (prediction) rep(NA, n.subgr))
  ##
  res$n.e <-
    c(if (common) replaceNULL(x$n.e.w, rep(NA, n.subgr)),
      if (random) replaceNULL(x$n.e.w, rep(NA, n.subgr)),
      if (prediction) replaceNULL(x$n.e.w, rep(NA, n.subgr)))
  ##
  res$n.c <-
    c(if (common) replaceNULL(x$n.c.w, rep(NA, n.subgr)),
      if (random) replaceNULL(x$n.c.w, rep(NA, n.subgr)),
      if (prediction) replaceNULL(x$n.c.w, rep(NA, n.subgr)))
  ##
  res$TE <-
    c(if (common) x$TE.common.w,
      if (random) x$TE.random.w,
      if (prediction) rep(NA, n.subgr))
  ##
  res$seTE <-
    c(if (common) x$seTE.common.w,
      if (random) x$seTE.random.w,
      if (prediction) rep(NA, n.subgr))
  ##
  res$lower <-
    c(if (common) x$lower.common.w,
      if (random) x$lower.random.w,
      if (prediction) x$lower.predict.w)
  ##
  res$upper <-
    c(if (common) x$upper.common.w,
      if (random) x$upper.random.w,
      if (prediction) x$upper.predict.w)
  ##
  res$statistic <-
    c(if (common) x$statistic.common.w,
      if (random) x$statistic.random.w,
      if (prediction) rep(NA, n.subgr))
  ##
  res$pval <-
      c(if (common) x$pval.common.w,
        if (random) x$pval.random.w,
        if (prediction) rep(NA, n.subgr))
  ##
  res$w.random <- res$w.common <-
    c(if (common) x$w.common.w,
      if (random) x$w.random.w,
      if (prediction) rep(NA, n.subgr))
  ##
  if (!is.null(names(res$TE)))
    names(res$seTE) <- names(res$lower) <- names(res$upper) <-
      names(res$statistic) <- names(res$pval) <- names(res$TE)
  ##
  res$k <- crp(x$k.w, common, random, prediction)
  res$k.study <- crp(x$k.study.w, common, random, prediction)
  res$k.all <- crp(x$k.all.w, common, random, prediction)
  res$k.TE <- crp(x$k.TE.w, common, random, prediction)
  ##
  res$Q <- c(x$Q.w, fillNAs)
  res$df.Q <- c(x$k.w - 1, fillNAs)
  res$pval.Q <- c(x$pval.Q.w, fillNAs)
  ##
  res$tau2 <- c(x$tau2.w, fillNAs)
  res$tau <- c(x$tau.w, fillNAs)
  ##
  res$H <- c(x$H.w, fillNAs)
  res$lower.H <- c(x$lower.H.w, fillNAs)
  res$upper.H <- c(x$upper.H.w, fillNAs)
  ##
  res$I2 <- c(x$I2.w, fillNAs)
  res$lower.I2 <- c(x$lower.I2.w, fillNAs)
  res$upper.I2 <- c(x$upper.I2.w, fillNAs)
  ##
  res$Rb <- c(x$Rb.w, fillNAs)
  res$lower.Rb <- c(x$lower.Rb.w, fillNAs)
  res$upper.Rb <- c(x$upper.Rb.w, fillNAs)
  ##
  res$Q.b <-
    c(if (common) addNAs2var(x$Q.b.common, n.subgr),
      if (random) addNAs2var(x$Q.b.random, n.subgr),
      if (prediction) rep(NA, n.subgr))
  ##
  res$df.Q.b <-
    c(if (common)
        addNAs2var(unlist(lapply(x$df.Q.b.common, compr)), n.subgr),
      if (random)
        addNAs2var(unlist(lapply(x$df.Q.b.random, compr)), n.subgr),
      if (prediction)
        rep(NA, n.subgr))
  ##
  res$pval.Q.b <-
    c(if (common) addNAs2var(x$pval.Q.b.common, n.subgr),
      if (random) addNAs2var(x$pval.Q.b.random, n.subgr),
      if (prediction) rep(NA, n.subgr))
  #
  if (debug)
    print(res)
  #
  res <- as.data.frame(res, row.names = seq_len(length(res$TE)))
  #
  res
}


overall2data <- function(x, common, random, prediction, name, subgroup,
                         debug = FALSE) {
  res <- list()
  ##
  res$studlab <-
    c(if (common)
        paste0(name,
              if (random | prediction) " (common effect)"),
      if (random)
        paste0(name,
              if (common | prediction) " (random effects)"),
      if (prediction)
        paste0(name,
              if (common | random) " (prediction)")
      )
  ##
  res$type <-
    c(if (common) "square",
      if (random) "square",
      if (prediction) "predict")
  ##
  res$model <-
    c(if (common) "common",
      if (random) "random",
      if (prediction) "predict")
  ##
  res$method <-
    c(if (common) x$method,
      if (random) x$method.random,
      if (prediction) x$method.predict)
  ##
  res$method.tau <-
    c(if (common) "",
      if (random) x$method.tau,
      if (prediction) "")
  ##
  res$three.level <-
    c(if (common) FALSE,
      if (random) x$three.level,
      if (prediction) FALSE)
  ##
  res$method.tau.ci <-
    c(if (common) "",
      if (random) x$method.tau.ci,
      if (prediction) "")
  ##
  res$tau.preset <-
    c(if (common) NA,
      if (random) replaceNULL(x$tau.preset, NA),
      if (prediction) NA)
  #
  res$method.I2 <-
    c(if (common) x$method.I2,
      if (random) x$method.I2,
      if (prediction) "")
  #
  res$method.common.ci <-
    c(if (common) x$method.common.ci,
      if (random) "",
      if (prediction) "")
  #
  res$method.random.ci <-
    c(if (common) "",
      if (random) x$method.random.ci,
      if (prediction) "")
  ##
  res$df.random <-
    c(if (common) NA,
      if (random) x$df.random,
      if (prediction) NA)
  ##
  res$adhoc.hakn.ci <-
    c(if (common) "",
      if (random) x$adhoc.hakn.ci,
      if (prediction) "")
  ##
  res$method.predict <-
    c(if (common) "",
      if (random) "",
      if (prediction) x$method.predict)
  ##
  res$df.predict <-
    c(if (common) NA,
      if (random) NA,
      if (prediction) x$df.predict)
  ##
  res$adhoc.hakn.pi <-
    c(if (common) "",
      if (random) "",
      if (prediction) x$adhoc.hakn.pi)
  ##
  res$rho <-
    c(if (common) NA,
      if (random) replaceNULL(x$rho, NA),
      if (prediction) NA)
  ##
  res$n.e <- sum(x$n.e, na.rm = TRUE)
  ##
  res$n.c <- sum(x$n.c, na.rm = TRUE)
  ##
  res$TE <-
    c(if (common) x$TE.common,
      if (random) x$TE.random,
      if (prediction) NA)
  ##
  res$seTE <-
    c(if (common) x$seTE.common,
      if (random) x$seTE.random,
      if (prediction) NA)
  ##
  res$lower <-
    c(if (common) x$lower.common,
      if (random) x$lower.random,
      if (prediction) x$lower.predict)
  ##
  res$upper <-
    c(if (common) x$upper.common,
      if (random) x$upper.random,
      if (prediction) x$upper.predict)
  ##
  res$statistic <-
    c(if (common) x$statistic.common,
      if (random) x$statistic.random,
      if (prediction) NA)
  ##
  res$pval <-
    c(if (common) x$pval.common,
      if (random) x$pval.random,
      if (prediction) NA)
  ##
  res$k <- crp(x$k, common, random, prediction)
  res$k.study <- crp(x$k.study, common, random, prediction)
  res$k.all <- crp(x$k.all, common, random, prediction)
  res$k.TE <- crp(x$k.TE, common, random, prediction)
  res$k.MH <- crp(x$k.MH, common, random, prediction)
  res$k0 <- crp(x$k0, common, random, prediction)
  ##
  sumcrp <- common + random + prediction
  ##
  res$Q <- c(x$Q[1], rep(NA, sumcrp - 1))
  res$df.Q <- c(x$df.Q[1], rep(NA, sumcrp - 1))
  res$pval.Q <- c(x$pval.Q[1], rep(NA, sumcrp - 1))
  ##
  res$tau2 <- c(x$tau2, rep(NA, sumcrp - 1))
  res$tau <- c(x$tau, rep(NA, sumcrp - 1))
  ##
  res$H <- c(x$H, rep(NA, sumcrp - 1))
  res$lower.H <- c(x$lower.H, rep(NA, sumcrp - 1))
  res$upper.H <- c(x$upper.H, rep(NA, sumcrp - 1))
  ##
  res$I2 <- c(x$I2, rep(NA, sumcrp - 1))
  res$lower.I2 <- c(x$lower.I2, rep(NA, sumcrp - 1))
  res$upper.I2 <- c(x$upper.I2, rep(NA, sumcrp - 1))
  ##
  res$Rb <- c(x$Rb, rep(NA, sumcrp - 1))
  res$lower.Rb <- c(x$lower.Rb, rep(NA, sumcrp - 1))
  res$upper.Rb <- c(x$upper.Rb, rep(NA, sumcrp - 1))
  #
  res$subgroup <- subgroup
  #
  if (debug)
    print(res)
  #
  res <- as.data.frame(res)
  #
  res
}


makeunique <- function(x, val = NA) {
  res <- unique(x)
  ##
  if (length(res) != 1)
    res <- val
  ##
  res
}


addNAs2var <- function(x, n) {
  if (n == 1)
    res <- x
  else {
    res <- vector("numeric", 0)
    ##
    for (i in seq_along(x))
      res <- c(res, c(x[i], rep(NA, n - 1)))
  }
  ##
  res
}


crp <- function(x, common, random, prediction, replace = NA) {
  x <- replaceNULL(x, replace)
  c(if (common) x, if (random) x, if (prediction) x)
}


compr <- function(x) {
  res <- paste(x, collapse = ", ")
  ##
  if (all(!grepl(",", res)))
    res <- as.numeric(res)
  ##
  res
}
