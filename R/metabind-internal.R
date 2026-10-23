extrMA <- function(id, x, name, type = "overall", outclab = NULL) {
  x <- x[[id]]
  #
  type <-
    setchar(type, c("overall", "prediction", "common", "random"), length = 1)
  #
  common <- type == "common"
  random <- type == "random"
  #
  metamiss <- inherits(x, "metamiss")
  metabin <- inherits(x, "metabin")
  metacont <- inherits(x, "metacont")
  metainc <- inherits(x, "metainc")
  metaprop <- inherits(x, "metaprop")
  metarate <- inherits(x, "metarate")
  #
  if (type == "overall") {
    res <- data.frame(
      id,
      class = setClass(x),
      name = name[id],
      subgroup = NA,
      #
      sm = x$sm,
      #
      common = x$common,
      random = x$random,
      prediction = x$prediction,
      #
      outclab = replaceNULL(outclab, x$outclab),
      #
      level.ma = x$level.ma,
      level.predict = x$level.predict,
      #
      three.level = replaceNULL(x$three.level, NA),
      rho = replaceNULL(x$rho, NA),
      #
      null.effect = x$null.effect,
      #
      title = x$title,
      complab = x$complab,
      #
      label.e = x$label.e,
      label.c = x$label.c,
      label.left = x$label.left,
      label.right = x$label.right,
      #
      method.I2 = x$method.I2,
      method.tau = x$method.tau,
      method.tau.ci = x$method.tau.ci,
      #
      tau.preset = replaceNULL(x$tau.preset),
      TE.tau = replaceNULL(x$TE.tau),
      tau.common = replaceNULL(x$tau.common)
    )
    #
    if (metabin) {
      res$MH.exact <- x$MH.exact
      #
      res$RR.Cochrane <- x$RR.Cochrane
      res$Q.Cochrane <- x$Q.Cochrane
      res$model.glmm <- x$model.glmm
      #
      res$incr <- format_seq(x$incr, exclude = 0, digits = 4)
      res$method.incr <- x$method.incr
      res$allstudies <- x$allstudies
      res$sparse <- x$sparse
      res$doublezeros <- x$doublezeros
    }
    #
    else if (metacont) {
      res$pooledvar <- x$pooledvar
      res$method.smd <- x$method.smd
      res$exact.smd <- x$exact.smd
      res$sd.glass <- x$sd.glass
      res$method.mean <- x$method.mean
      res$method.sd <- x$method.sd
    }
    #
    else if (metainc) {
      res$model.glmm <- x$model.glmm
      #
      res$incr <- format_seq(x$incr, exclude = 0, digits = 4)
      res$method.incr <- x$method.incr
      res$sparse <- x$sparse
    }
    #
    else if (metaprop || is_prop(x$sm)) {
      res$pscale <- x$pscale
      #
      res$incr <- format_seq(x$incr, exclude = 0, digits = 4)
      res$method.incr <- x$method.incr
    }
    #
    else if (metarate || is_rate(x$sm)) {
      res$irscale <- x$irscale
      res$irunit <- x$irunit
      #
      res$incr <- format_seq(x$incr, exclude = 0, digits = 4)
      res$method.incr <- x$method.incr
    }
    #
    if (metamiss) {
      res$class <- "metamiss"
      res$method.miss <- x$method.miss
      #
      res$IMOR.e <- format_seq(x$IMOR.e, digits = 4)
      res$IMOR.c <- format_seq(x$IMOR.c, digits = 4)
    }
    #
    res$k <- x$k
    res$k.study <- x$k.study
    res$k.all <- x$k.all
    res$k.TE <- x$k.TE
    res$k.MH <- replaceNULL(x$k.MH)
    res$k0 <- replaceNULL(x$k0)
    #
    res$n.harmonic.mean <- replaceNULL(x$n.harmonic.mean)
    res$t.harmonic.mean <- replaceNULL(x$t.harmonic.mean)
    #
    res$n.e <- sum(replaceNULL(x$n.e), na.rm = TRUE)
    res$n.c <- sum(replaceNULL(x$n.c), na.rm = TRUE)
    #
    res$Q <- x$Q[1]
    res$df.Q <- x$df.Q[1]
    res$pval.Q <- x$pval.Q[1]
    #
    res$tau2 <- x$tau2
    res$lower.tau2 <- x$lower.tau2
    res$upper.tau2 <- x$upper.tau2
    res$se.tau2 <- replaceNULL(x$se.tau2)
    res$tau <- x$tau
    res$lower.tau <- x$lower.tau
    res$upper.tau <- x$upper.tau
    #
    res$H <- x$H
    res$lower.H <- x$lower.H
    res$upper.H <- x$upper.H
    res$I2 <- x$I2
    res$lower.I2 <- x$lower.I2
    res$upper.I2 <- x$upper.I2
    res$Rb <- x$Rb
    res$lower.Rb <- x$lower.Rb
    res$upper.Rb <- x$upper.Rb
  }
  else if (common || random) {
    res <- data.frame(
      id,
      #
      model = if (common) "common" else "random",
      method = if (common) x$method else x$method.random,
      #
      method.interval = if (common) x$method.common.ci else x$method.random.ci,
      adhoc.hakn = if (common) NA else x$adhoc.hakn.ci
      )
    #
    if (common) {
      res$adhoc.hakn <- NULL
      #
      res$TE <- x$TE.common
      res$seTE <- x$seTE.common
      res$lower <- x$lower.common
      res$upper <- x$upper.common
      res$statistic <- x$statistic.common
      res$pval <- x$pval.common
      res$text <- x$text.common
    }
    else {
      res$TE <- x$TE.random
      res$seTE <- x$seTE.random
      res$lower <- x$lower.random
      res$upper <- x$upper.random
      res$statistic <- x$statistic.random
      res$pval <- x$pval.random
      res$text <- x$text.random
      #
      res$df <- replaceNULL(x$df.random)
      res$df.hakn <- replaceNULL(x$df.hakn.ci)
      res$df.kero <- replaceNULL(x$df.kero)
    }
  }
  else if (type == "prediction") {
    res <- data.frame(
      id,
      model = "predict",
      #
      method.interval = x$method.predict,
      adhoc.hakn = x$adhoc.hakn.pi,
      #
      seTE = x$seTE.predict,
      lower = x$lower.predict,
      upper = x$upper.predict,
      #
      df = x$df.predict,
      df.hakn = replaceNULL(x$df.hakn.pi),
      df.kero = replaceNULL(x$df.kero)
    )
  }
  #
  res
}

subgr2meta <- function(x, common, random, prediction, name) {
  n.subgr <- length(x$subgroup.levels)
  #
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
    #
    subgroup = NA,
    #
    subgroup.levels = name,
    #
    TE =
      c(if (common) x$TE.common.w,
        if (random) x$TE.random.w,
        if (prediction) rep(NA, n.subgr)),
    #
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
    #
    TE.common = x$TE.common,
    seTE.common = x$seTE.common,
    lower.common = x$lower.common,
    upper.common = x$upper.common,
    statistic.common = x$statistic.common,
    pval.common = x$pval.common,
    text.common = if (common) x$text.common else NULL,
    #
    TE.random = x$TE.random,
    seTE.random = x$seTE.random,
    lower.random = x$lower.random,
    upper.random = x$upper.random,
    statistic.random = x$statistic.random,
    pval.random = x$pval.random,
    text.random = if (random) x$text.random else NULL,
    #
    lower.predict = x$lower.predict,
    upper.predict = x$upper.predict,
    text.predict = if (prediction) x$text.predict else NULL,
    #
    df.random = replaceNULL(x$df.random),
    df.hakn = replaceNULL(x$df.hakn),
    df.kero = replaceNULL(x$df.kero),
    #
    n.harmonic.mean = replaceNULL(x$n.harmonic.mean.w),
    t.harmonic.mean = replaceNULL(x$t.harmonic.mean.w),
    #
    n.e = replaceNULL(x$n.e),
    n.c = replaceNULL(x$n.c),
    #
    Q = x$Q[1],
    df.Q = x$df.Q[1],
    pval.Q = x$pval.Q[1],
    #
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
    #
    Q.b.common = NA,
    Q.b.random = NA,
    df.Q.b = NA,
    pval.Q.b.common = NA,
    pval.Q.b.random = NA)
  #
  res$subgroup <- rep_len(name, length(res$TE))
  #
  res
}

makeunique <- function(x, val = NA) {
  res <- unique(x)
  #
  if (length(res) != 1)
    res <- val
  #
  res
}

addNAs2var <- function(x, n) {
  if (n == 1)
    res <- x
  else {
    res <- vector("numeric", 0)
    #
    for (i in seq_along(x))
      res <- c(res, c(x[i], rep(NA, n - 1)))
  }
  #
  res
}

setClass <- function(x) {
  if (inherits(x, "trimfill"))
    class <- "trimfill"
  else
    class <- class(x)
  #
  class <- class[class != "metamiss"]
  class <- class[class != "meta"]
  #
  if (length(class) > 1) {
    warning("Meta-analysis has ", length(class), " classes: ",
            paste(class, collapse = ", "), "\n Keeping the first class: ",
            class[1], call. = FALSE)
    #
    class <- class[1]
  }
  #
  class
}

drop_allNAs <- function(x) {
  all_NA <- colSums(is.na(x)) == nrow(x)
  #
  x[, !all_NA, drop = FALSE]
}
