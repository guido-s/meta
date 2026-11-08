extrMA <- function(id, x, name, type = "overall", outclab = NULL) {
  x <- x[[id]]
  #
  type <-
    setchar(type, c("overall", "common", "random", "prediction"), length = 1)
  #
  common <- type == "common"
  random <- type == "random"
  #
  metagen <- inherits(x, "metagen")
  metabin <- inherits(x, "metabin")
  metacont <- inherits(x, "metacont")
  metainc <- inherits(x, "metainc")
  metaprop <- inherits(x, "metaprop")
  metarate <- inherits(x, "metarate")
  trimfill <- inherits(x, "trimfill")
  metamiss <- inherits(x, "metamiss")
  #
  if (type == "overall") {
    res <- data.frame(
      id,
      name = name[id],
      class = setClass(x),
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
      level.hetstat = x$level.hetstat,
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
    if (metagen) {
      res$func.transf <- replaceNULL(x$func.transf)
      res$args.transf <- replaceNULL(x$args.transf)
      res$func.backtransf <- replaceNULL(x$func.backtransf)
      res$args.backtransf <- replaceNULL(x$args.backtransf)
    }
    else if (metabin) {
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
      #
      res$phi <- replaceNULL(x$phi)
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
      res$sparse <- x$sparse
    }
    #
    else if (metarate || is_rate(x$sm)) {
      res$irscale <- x$irscale
      res$irunit <- x$irunit
      #
      res$incr <- format_seq(x$incr, exclude = 0, digits = 4)
      res$method.incr <- x$method.incr
      res$sparse <- x$sparse
    }
    #
    if (metamiss) {
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
    #
    if (trimfill) {
      res$k0 <- replaceNULL(x$k0)
      res$type.tf <- replaceNULL(x$type, "")
    }
    #
    res$n.harmonic.mean <- replaceNULL(x$n.harmonic.mean)
    res$t.harmonic.mean <- replaceNULL(x$t.harmonic.mean)
    #
    res$n.e <- sum_or_na(x$n.e)
    res$n.c <- sum_or_na(x$n.c)
    res$n <- sum_or_na(x$n)
    #
    res$event.e <- sum_or_na(x$event.e)
    res$event.c <- sum_or_na(x$event.c)
    res$event <- sum_or_na(x$event)
    #
    res$time.e <- sum_or_na(x$time.e)
    res$time.c <- sum_or_na(x$time.c)
    res$time <- sum_or_na(x$time)
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
  else if (common) {
    res <- data.frame(
      id,
      type = "Common effect model",
      method = x$method,
      #
      method.interval = x$method.common.ci,
      adhoc.hakn = NA,
      #
      df = NA,
      df.hakn = NA,
      df.kero = NA,
      #
      TE = x$TE.common,
      seTE = x$seTE.common,
      lower = x$lower.common,
      upper = x$upper.common,
      statistic = x$statistic.common,
      pval = x$pval.common,
      #
      text = x$text.common
    )
  }
  else if (random) {
    res <- data.frame(
      id,
      type = "Random effects model",
      method = x$method.random,
      #
      method.interval = x$method.random.ci,
      adhoc.hakn = x$adhoc.hakn.ci,
      #
      df = replaceNULL(x$df.random),
      df.hakn = replaceNULL(x$df.hakn.ci),
      df.kero = replaceNULL(x$df.kero),
      #
      TE = x$TE.random,
      seTE = x$seTE.random,
      lower = x$lower.random,
      upper = x$upper.random,
      statistic = x$statistic.random,
      pval = x$pval.random,
      #
      text = x$text.random
    )
  }
  else if (type == "prediction") {
    res <- data.frame(
      id,
      type = "Prediction interval",
      method = x$method.random,
      #
      method.interval = x$method.predict,
      adhoc.hakn = x$adhoc.hakn.pi,
      #
      df = x$df.predict,
      df.hakn = replaceNULL(x$df.hakn.pi),
      df.kero = replaceNULL(x$df.kero),
      #
      seTE = x$seTE.predict,
      lower = x$lower.predict,
      upper = x$upper.predict,
      #
      text = x$text.predict
    )
  }
  #
  res
}

extrSub <- function(id, x, type = "common") {
  x <- x[[id]]
  #
  n.by <- length(x$subgroup.levels)
  #
  type <- setchar(type, c("common", "random", "prediction"), length = 1)
  #
  if (is.numeric(x$subgroup.levels))
    x$subgroup.levels <- as.character(x$subgroup.levels)
  #
  if (type == "common") {
    res <- data.frame(
      id = id,
      #
      subgroup.name = x$subgroup.name,
      subgroup.levels = x$subgroup.levels,
      k.w = x$k.w,
      #
      TE = matrix(list2mat(x$TE.common.w), ncol = 1),
      seTE = matrix(list2mat(x$seTE.common.w), ncol = 1),
      lower = matrix(list2mat(x$lower.common.w), ncol = 1),
      upper = matrix(list2mat(x$upper.common.w), ncol = 1),
      statistic = matrix(list2mat(x$statistic.common.w), ncol = 1),
      pval = matrix(list2mat(x$pval.common.w), ncol = 1),
      #
      Q.b = x$Q.b.common,
      df.Q.b = x$df.Q.b.common,
      pval.Q.b = x$pval.Q.b.common
    )
  }
  else if (type == "random") {
    method <- rep(colnames(x$lower.random.w),
                  rep(nrow(x$lower.random.w),
                      ncol(x$lower.random.w)))
    method <- replaceNULL(method, x$method.random.ci)
    #
    res <- data.frame(
      id = id,
      #
      method.interval = method,
      adhoc.hakn = x$adhoc.hakn.ci,
      #
      subgroup.name = x$subgroup.name,
      subgroup.levels = x$subgroup.levels,
      k.w = x$k.w,
      #
      TE = matrix(list2mat(x$TE.random.w), ncol = 1),
      seTE = matrix(list2mat(x$seTE.random.w), ncol = 1),
      lower = matrix(list2mat(x$lower.random.w), ncol = 1),
      upper = matrix(list2mat(x$upper.random.w), ncol = 1),
      df = matrix(list2mat(x$df.random.w), ncol = 1),
      statistic = matrix(list2mat(x$statistic.random.w), ncol = 1),
      pval = matrix(list2mat(x$pval.random.w), ncol = 1),
      #
      Q.b = x$Q.b.random,
      df.Q.b = x$df.Q.b.random,
      pval.Q.b = x$pval.Q.b.random
    )
  }
  else if (type == "prediction") {
    method <- rep(colnames(x$lower.predict.w),
                  rep(nrow(x$lower.predict.w),
                      ncol(x$lower.predict.w)))
    method <- replaceNULL(method, x$method.predict)
    #
    res <- data.frame(
      id = id,
      #
      method.interval = method,
      adhoc.hakn = x$adhoc.hakn.pi,
      #
      subgroup.name = x$subgroup.name,
      subgroup.levels = x$subgroup.levels,
      k.w = x$k.w,
      #
      lower = matrix(list2mat(x$lower.predict.w), ncol = 1),
      upper = matrix(list2mat(x$upper.predict.w), ncol = 1),
      df = matrix(list2mat(x$df.predict.w), ncol = 1)
    )
  }
  #
  res$n.harmonic.mean <- replaceNULL(x$n.harmonic.mean.w)
  res$t.harmonic.mean <- replaceNULL(x$t.harmonic.mean.w)
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
    n = replaceNULL(x$n),
    #
    event.e = replaceNULL(x$event.e),
    event.c = replaceNULL(x$event.c),
    event = replaceNULL(x$event),
    #
    time.e = replaceNULL(x$time.e),
    time.c = replaceNULL(x$time.c),
    time = replaceNULL(x$time),
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

makeunique <- function(x, val = NA, orig = FALSE) {
  res <- unique(x)
  #
  if (length(res) != 1) {
    if (orig)
      return(res)
    else
      res <- val
  }
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
  else if (inherits(x, "metamiss"))
    class <- "metamiss"
  else
    class <- class(x)
  #
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
