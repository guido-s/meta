setlistname <- function(x, name) {
  if (!is.list(x))
    x <- list(x)
  #
  if (!is.null(name) & is.null(names(x)))
    names(x) <- name
  #
  x
}

setname <- function(x, name) {
  if (!is.null(x) &&
      is.null(names(x)) &&
      !(is.null(name) || all(name == ""))) {
    nx <- length(x)
    #
    if (length(name) == nx)
      names(x) <- name
    else if (length(name) == 1)
      names(x) <- rep_len(name, nx)
  }
  #
  x
}

samedata <- function(x, y, stop = TRUE) {
  mismatch <- FALSE
  #
  if (!is.null(x$data) & !is.null(y$data)) {
    #
    if (nrow(x$data) != nrow(y$data)) {
      if (stop)
        stop("Meta-analyses based on different data sets.",
             call. = FALSE)
      else
        return(TRUE)
    }
    #
    if (inherits(x, "metabin")) {
      mismatch1 <-
        mismatch(x, y, ".event.e")
      mismatch2 <-
        mismatch(x, y, ".n.e")
      mismatch3 <-
        mismatch(x, y, ".event.c")
      mismatch4 <-
        mismatch(x, y, ".n.c")
      #
      mismatch <- mismatch1 | mismatch2 | mismatch3 | mismatch4
    }
    else if (inherits(x, "metacont")) {
      mismatch1 <-
        mismatch(x, y, ".n.e")
      mismatch2 <-
        mismatch(x, y, ".mean.e")
      mismatch3 <-
        mismatch(x, y, ".sd.e")
      mismatch4 <-
        mismatch(x, y, ".n.c")
      mismatch5 <-
        mismatch(x, y, ".mean.c")
      mismatch6 <-
        mismatch(x, y, ".sd.c")
      #
      mismatch <-
        mismatch1 | mismatch2 | mismatch3 | mismatch4 |
        mismatch5 | mismatch6
    }
    else if (inherits(x, "metacor")) {
      mismatch1 <-
        mismatch(x, y, ".cor")
      mismatch2 <-
        mismatch(x, y, ".n")
      #
      mismatch <- mismatch1 | mismatch2
    }
    else if (inherits(x, "metagen")) {
      #
      mismatch1 <-
        mismatch(x, y, ".n.e")
      mismatch2 <-
        mismatch(x, y, ".n.c")
      #
      mismatch <- mismatch1 | mismatch2
    }
    else if (inherits(x, "metainc")) {
      mismatch1 <-
        mismatch(x, y, ".time.e")
      mismatch2 <-
        mismatch(x, y, ".n.e")
      mismatch3 <-
        mismatch(x, y, ".time.c")
      mismatch4 <-
        mismatch(x, y, ".n.c")
      #
      mismatch <- mismatch1 | mismatch2 | mismatch3 | mismatch4
    }
    else if (inherits(x, "metamean")) {
      mismatch1 <-
        mismatch(x, y, ".n")
      mismatch2 <-
        mismatch(x, y, ".mean")
      mismatch3 <-
        mismatch(x, y, ".sd")
      #
      mismatch <-
        mismatch1 | mismatch2 | mismatch3
    }
    else if (inherits(x, "metaprop")) {
      mismatch1 <-
        mismatch(x, y, ".event")
      mismatch2 <-
        mismatch(x, y, ".n")
      #
      mismatch <- mismatch1 | mismatch2
    }
    else if (inherits(x, "metarate")) {
      mismatch1 <-
        mismatch(x, y, ".time")
      mismatch2 <-
        mismatch(x, y, ".n")
      #
      mismatch <- mismatch1 | mismatch2
    }
    #
    if (mismatch & stop)
      stop("Meta-analyses based on different data sets.",
           call. = FALSE)
  }
  #
  invisible(mismatch)
}

samesm <- function(x, y) {
  if (!is.null(x$sm) & !is.null(y$sm)) {
    if (is_relative_effect(x$sm) != is_relative_effect(y$sm))
      stop("Summary measures used in meta-analyses do not fit.",
           call. = FALSE)
  }
  #
  if (inherits(x, "metabin")) {
    if ((x$sm != y$sm) &
        any(c(x$sm, y$sm) %in% c("RD", "ASD")))
      stop("Summary measures used in meta-analyses do not fit.",
           call. = FALSE)
  }
  #
  invisible(NULL)
}

samesubgroups <- function(x, y) {
  if (!is.null(x$subgroup.levels) &
      !is.null(y$subgroup.levels)) {
    #
    if (length(x$subgroup.levels) !=
        length(y$subgroup.levels))
      stop("Meta-analyses have different number of subgroups.",
           call. = FALSE)
    #
    if (any(x$subgroup.levels != y$subgroup.levels))
      stop("Meta-analyses based on different subgroup-analyses.",
           call. = FALSE)
  }
  #
  invisible(NULL)
}

updateobj <- function(x,
                      label.common, label.random, label.predict,
                      hetlabel, taulabel, label.subgroup,
                      text.common, text.random, text.predict,
                      text.w.common, text.w.random) {
  
  is.merge <- inherits(x, "metamerge")
  is.copas <- !is.merge & inherits(x, "copas")
  is.limit <- !is.merge & inherits(x, "limitmeta")
  is.robu  <- !is.merge & inherits(x, "robu")
  is.tf    <- !is.merge & inherits(x, "trimfill")
  #
  res <- x
  #
  if (is.null(x$hetlabel) || all(x$hetlabel == ""))
    res$hetlabel <- replaceNULL(hetlabel, "")
  #
  if (is.null(x$detail.tau) || all(x$detail.tau == ""))
    res$detail.tau <- replaceNULL(taulabel, "")
  #
  if (!is.merge) {    
    res$TE.common <- setname(res$TE.common, label.common)
    res$seTE.common <- setname(res$seTE.common, label.common)
    res$statistic.common <- setname(res$statistic.common, label.common)
    res$pval.common <- setname(res$pval.common, label.common)
    res$lower.common <- setname(res$lower.common, label.common)
    res$upper.common <- setname(res$upper.common, label.common)
    res$zval.common <- setname(res$zval.common, label.common)
    #
    res$TE.random <- setname(res$TE.random, label.random)
    res$seTE.random <- setname(res$seTE.random, label.random)
    res$statistic.random <- setname(res$statistic.random, label.random)
    res$pval.random <- setname(res$pval.random, label.random)
    res$method.random.ci <- setname(res$method.random.ci, label.random)
    res$df.random <- setname(res$df.random, label.random)
    res$lower.random <- setname(res$lower.random, label.random)
    res$upper.random <- setname(res$upper.random, label.random)
    res$zval.random <- setname(res$zval.random, label.random)
    res$seTE.classic <- setname(res$seTE.classic, label.random)
    res$adhoc.hakn.ci <- setname(res$adhoc.hakn.ci, label.random)
    res$df.hakn <- setname(res$df.hakn, label.random)
    res$seTE.hakn.ci <- setname(res$seTE.hakn.ci, label.random)
    res$seTE.hakn.adhoc.ci <- setname(res$seTE.hakn.adhoc.ci, label.random)
    res$df.kero <- setname(res$df.kero, label.random)
    res$seTE.kero <- setname(res$seTE.kero, label.random)
    #
    res$seTE.predict <- setname(res$seTE.predict, label.predict)
    res$df.predict <- setname(res$df.predict, label.predict)
    res$lower.predict <- setname(res$lower.predict, label.predict)
    res$upper.predict <- setname(res$upper.predict, label.predict)
    res$adhoc.hakn.pi <- setname(res$adhoc.hakn.pi, label.predict)
    res$seTE.hakn.pi <- setname(res$seTE.hakn.pi, label.predict)
    res$seTE.hakn.adhoc.pi <- setname(res$seTE.hakn.adhoc.pi, label.predict)
    #
    res$Q <- setname(res$Q, hetlabel)
    res$df.Q <- setlistname(res$df.Q, hetlabel)
    res$pval.Q <- setname(res$pval.Q, hetlabel)
    #
    res$I2 <- setname(res$I2, hetlabel)
    res$lower.I2 <- setname(res$lower.I2, hetlabel)
    res$upper.I2 <- setname(res$upper.I2, hetlabel)
    #
    res$H <- setname(res$H, hetlabel)
    res$lower.H <- setname(res$lower.H, hetlabel)
    res$upper.H <- setname(res$upper.H, hetlabel)
    #
    res$Rb <- setname(res$Rb, hetlabel)
    res$lower.Rb <- setname(res$lower.Rb, hetlabel)
    res$upper.Rb <- setname(res$upper.Rb, hetlabel)
    #
    res$TE.common.w <- setlistname(res$TE.common.w, label.subgroup)
    res$seTE.common.w <- setlistname(res$seTE.common.w, label.subgroup)
    res$statistic.common.w <-
      setlistname(res$statistic.common.w, label.subgroup)
    res$pval.common.w <- setlistname(res$pval.common.w, label.subgroup)
    res$lower.common.w <- setlistname(res$lower.common.w, label.subgroup)
    res$upper.common.w <- setlistname(res$upper.common.w, label.subgroup)
    #
    res$w.common.w <- setlistname(res$w.common.w, label.subgroup)
    #
    res$TE.random.w <- setlistname(res$TE.random.w, label.subgroup)
    res$seTE.random.w <- setlistname(res$seTE.random.w, label.subgroup)
    res$statistic.random.w <-
      setlistname(res$statistic.random.w, label.subgroup)
    res$pval.random.w <- setlistname(res$pval.random.w, label.subgroup)
    res$df.random.w <- setlistname(res$df.random.w, label.subgroup)
    res$lower.random.w <- setlistname(res$lower.random.w, label.subgroup)
    res$upper.random.w <- setlistname(res$upper.random.w, label.subgroup)
    res$df.hakn.w <- setlistname(res$df.hakn.w, label.subgroup)
    res$df.kero.w <- setlistname(res$df.kero.w, label.subgroup)
    #
    res$w.random.w <- setlistname(res$w.random.w, label.subgroup)
    #
    res$seTE.predict.w <- setlistname(res$seTE.predict.w, label.subgroup)
    res$df.predict.w <- setlistname(res$df.predict.w, label.subgroup)
    res$lower.predict.w <- setlistname(res$lower.predict.w, label.subgroup)
    res$upper.predict.w <- setlistname(res$upper.predict.w, label.subgroup)
    #
    res$Q.w <- setname(res$Q.w, label.subgroup)
    res$tau.w <- setlistname(res$tau.w, label.subgroup)
    #
    res$Q.b.common <- setname(res$Q.b.common, label.subgroup)
    res$df.Q.b.common <- setlistname(res$df.Q.b.common, label.subgroup)
    res$pval.Q.b.common <- setname(res$pval.Q.b.common, label.subgroup)
    #
    res$Q.b.random <- setname(res$Q.b.random, label.subgroup)
    res$df.Q.b.random <- setlistname(res$df.Q.b.random, label.subgroup)
    res$pval.Q.b.random <- setname(res$pval.Q.b.random, label.subgroup)
  }
  #
  # Act upon ordinary meta-analysis object
  #
  if (!(is.copas | is.limit | is.robu | is.tf)) {
    if (!(is.null(label.common) || label.common == "") &
        is.null(text.common))
      res$text.common <- paste0(res$text.common, " (", label.common, ")")
    else if (!is.null(text.common))
      res$text.common <- text.common
    #
    if (!is.null(label.random) & is.null(text.random))
      res$text.random <- paste0(res$text.random, " (", label.random, ")")
    else if (!is.null(text.random))
      res$text.random <- text.random
    #
    if (!is.null(label.predict) & is.null(text.predict))
      res$text.predict <- paste0(res$text.predict, " (", label.predict, ")")
    else if (!is.null(text.predict))
      res$text.predict <- text.predict
    #
    if (!is.null(label.common) & is.null(text.w.common))
      res$text.w.common <- label.common
    else if (!is.null(text.w.common))
      res$text.w.common <- text.w.common
    #
    if (!is.null(label.random) & is.null(text.w.random))
      res$text.w.random <- label.random
    else if (!is.null(text.w.random))
      res$text.w.random <- text.w.random
    #
    return(res)
  }
  #
  # Other objects
  #
  res$text.predict <- ""
  #
  if (is.copas | is.limit | is.robu) {
    res$method.random <- "Inverse"
    #
    res$method.random.ci <- "classic"
    res$adhoc.hakn.ci <- ""
    res$df.random <- Inf
    #
    res$method.predict <- NULL
    res$adhoc.hakn.pi <- NULL
    res$seTE.predict <- NULL
    res$df.predict <- NULL
    res$level.predict <- NULL
    res$lower.predict <- NULL
    res$upper.predict <- NULL
    res$seTE.hakn.pi <- NULL
    res$seTE.hakn.adhoc.pi <- NULL
  }
  #
  if (is.copas) {
    if (is.null(label.random))
      label.random <- "copas"
    #
    if (is.null(hetlabel))
      hetlabel <- "copas"
    #
    if (is.null(taulabel))
      taulabel <- "copas"
    #
    res$TE.random <- setname(res$TE.adjust, label.random)
    res$seTE.random <- setname(res$seTE.adjust, label.random)
    res$lower.random <- setname(res$lower.adjust, label.random)
    res$upper.random <- setname(res$upper.adjust, label.random)
    res$statistic.random <- setname(res$statistic.adjust, label.random)
    res$pval.random <- setname(res$pval.adjust, label.random)
    #
    res$w.random <- rep(0, length(res$w.random))
    #
    res$tau <- setname(res$tau.adjust, taulabel)
    res$lower.tau <- NA
    res$upper.tau <- NA
    res$tau2 <- setname(res$tau.adjust^2, taulabel)
    res$lower.tau2 <- NA
    res$upper.tau2 <- NA
    res$se.tau <- NA
    #
    res$method.tau <- "ML"
    #
    res$hetlabel <- hetlabel
    res$detail.tau <- taulabel
    #
    if (is.null(text.random))
      res$text.random <- "Copas selection model"
    else
      res$text.random <- text.random
    #
    if (is.null(text.w.random))
      res$text.w.random <- "Copas"
    else
      res$text.w.random <- text.w.random
    #
    res$k.study <- sum(is.finite(res$TE))
  }
  else if (is.limit) {
    if (is.null(label.random))
      label.random <- "limit"
    #
    if (is.null(hetlabel))
      hetlabel <- "limit"
    #
    if (is.null(taulabel))
      taulabel <- "limit"
    #
    res$TE.random <- setname(res$TE.adjust, label.random)
    res$seTE.random <- setname(res$seTE.adjust, label.random)
    res$lower.random <- setname(res$lower.adjust, label.random)
    res$upper.random <- setname(res$upper.adjust, label.random)
    res$statistic.random <- setname(res$statistic.adjust, label.random)
    res$pval.random <- setname(res$pval.adjust, label.random)
    #
    res$w.random <- rep(0, length(res$w.random))
    #
    res$tau <- setname(res$tau, taulabel)
    res$lower.tau <- NA
    res$upper.tau <- NA
    res$tau2 <- setname(res$tau2, taulabel)
    res$lower.tau2 <- NA
    res$upper.tau2 <- NA
    res$se.tau <- NA
    #
    res$method.tau <- res$x$method.tau
    #
    res$hetlabel <- hetlabel
    res$detail.tau <- taulabel
    #
    if (is.null(text.random))
      res$text.random <- "Limit meta-analysis"
    else
      res$text.random <- text.random
    #
    if (is.null(text.w.random))
      res$text.w.random <- "limit"
    else
      res$text.w.random <- text.w.random
    #
    res$k.study <- res$k
  }
  else if (is.robu) {
    if (is.null(label.random))
      label.random <- "RVE"
    #
    if (is.null(hetlabel))
      hetlabel <- "RVE"
    #
    if (is.null(taulabel))
      taulabel <- "RVE"
    #
    res$TE.random <- setname(res$reg_table$b.r[1], label.random)
    res$seTE.random <- setname(res$reg_table$SE[1], label.random)
    res$lower.random <- setname(res$reg_table$CI.L[1], label.random)
    res$upper.random <- setname(res$reg_table$CI.U[1], label.random)
    res$statistic.random <- setname(res$reg_table$t[1], label.random)
    res$pval.random <- setname(res$reg_table$prob[1], label.random)
    #
    res$w.random <- res$data.full$r.weights
    #
    res$level.ma <- 0.95
    #
    res$tau <- setname(sqrt(res$mod_info$tau.sq), taulabel)
    res$lower.tau <- NA
    res$upper.tau <- NA
    res$tau2 <- setname(res$mod_info$tau.sq, taulabel)
    res$lower.tau2 <- NA
    res$upper.tau2 <- NA
    res$se.tau <- NA
    #
    res$method.tau <- "DL"
    #
    res$hetlabel <- hetlabel
    res$detail.tau <- taulabel
    #
    if (is.null(text.random))
      res$text.random <- "RVE model"
    else
      res$text.random <- text.random
    #
    if (is.null(text.w.random))
      res$text.w.random <- "RVE"
    else
      res$text.w.random <- text.w.random
  }
  else if (is.tf) {
    if (is.null(hetlabel))
      hetlabel <- "TF"
    #
    if (is.null(taulabel))
      taulabel <- "TF"
    #
    if (!(is.null(label.common) || label.common == "") &
        is.null(text.common))
      res$text.common <- paste0(res$text.common, " (", label.common, ")")
    else if (!is.null(text.common))
      res$text.common <- text.common
    else
      res$text.common <- "Trim-and-fill method (CE)"
    #
    if (!is.null(label.random) & is.null(text.random))
      res$text.random <- paste0(res$text.random, " (", label.random, ")")
    else if (!is.null(text.random))
      res$text.random <- text.random
    else
      res$text.random <- "Trim-and-fill method (RE)"
    #
    res$hetlabel <- hetlabel
    res$detail.tau <- taulabel
    #
    if ((!is.null(label.common) & is.null(text.w.common)) ||
        !is.null(text.w.common))
      res$text.w.common <- text.w.common
    else
      res$text.w.common <- "CE-TF"
    #
    if ((!is.null(label.random) & is.null(text.w.random)) ||
        !is.null(text.w.random))
      res$text.w.random <- text.w.random
    else
      res$text.w.random <- "RE-TF"
  }
  
  res
}

dropcommon <- function(x) {
  res <- x
  #
  res$method <- NULL
  #
  res$w.common <- NULL
  res$TE.common <- NULL
  res$seTE.common <- NULL
  res$statistic.common <- NULL
  res$pval.common <- NULL
  res$lower.common <- NULL
  res$upper.common <- NULL
  res$zval.common <- NULL
  res$text.common <- NULL
  #
  res$k.MH <- NULL
  #
  res$TE.common.w <- NULL
  res$seTE.common.w <- NULL
  res$statistic.common.w <- NULL
  res$pval.common.w <- NULL
  res$lower.common.w <- NULL
  res$upper.common.w <- NULL
  res$w.common.w <- NULL
  #
  res$Q.w.common <- NULL
  res$pval.Q.w.common <- NULL
  #
  res$Q.b.common <- NULL
  res$df.Q.b.common <- NULL
  res$pval.Q.b.common <- NULL
  #
  res$zval.common.w <- NULL
  res$TE.fixed.w <- NULL
  res$seTE.fixed.w <- NULL
  res$lower.fixed.w <- NULL
  res$upper.fixed.w <- NULL
  res$statistic.fixed.w <- NULL
  res$pval.fixed.w <- NULL
  res$zval.fixed.w <- NULL
  res$w.fixed.w <- NULL
  res$Q.w.fixed <- NULL
  res$pval.Q.w.fixed <- NULL
  res$Q.b.fixed <- NULL
  res$pval.Q.b.fixed <- NULL
  #
  res
}

droprandom <- function(x) {
  res <- x
  #
  res$method.random <- NULL
  #
  res$w.random <- NULL
  res$TE.random <- NULL
  res$seTE.random <- NULL
  res$statistic.random <- NULL
  res$pval.random <- NULL
  res$method.random.ci <- NULL
  res$df.random <- NULL
  res$lower.random <- NULL
  res$upper.random <- NULL
  res$zval.random <- NULL
  #
  res$seTE.classic <- NULL
  #
  res$adhoc.hakn.ci <- NULL
  res$df.hakn <- NULL
  res$seTE.hakn.ci <- NULL
  res$seTE.hakn.adhoc.ci <- NULL
  #
  res$df.kero <- NULL
  res$seTE.kero <- NULL
  #
  res$text.random <- NULL
  #
  res$cluster <- FALSE
  res$three.level <- FALSE
  #
  res$k.study <- res$k
  #
  res$TE.random.w <- NULL
  res$seTE.random.w <- NULL
  res$statistic.random.w <- NULL
  res$pval.random.w <- NULL
  res$df.random.w <- NULL
  res$lower.random.w <- NULL
  res$upper.random.w <- NULL
  res$w.random.w <- NULL
  #
  res$seTE.classic.w <- NULL
  #
  res$df.hakn.ci.w <- NULL
  res$seTE.hakn.ci.w <- NULL
  res$seTE.hakn.adhoc.ci.w <- NULL
  #
  res$df.kero.w <- NULL
  res$seTE.kero.w <- NULL
  #
  res$Q.w.random <- NULL
  res$pval.Q.w.random <- NULL
  #
  res$Q.b.random <- NULL
  res$df.Q.b.random <- NULL
  res$pval.Q.b.random <- NULL
  #
  res$zval.random.w <- NULL
  #
  res
}

droppredict <- function(x) {
  res <- x
  #
  res$method.predict <- NULL
  res$adhoc.hakn.pi <- NULL
  res$df.hakn.pi <- NULL
  res$seTE.predict <- NULL
  res$df.predict <- NULL
  res$lower.predict <- NULL
  res$upper.predict <- NULL
  res$seTE.hakn.pi <- NULL
  res$seTE.hakn.adhoc.pi <- NULL
  res$text.predict <- NULL
  #      
  res$seTE.predict.w <- NULL
  res$df.predict.w <- NULL
  res$lower.predict.w <- NULL
  res$upper.predict.w <- NULL
  res$seTE.hakn.pi.w <- NULL
  res$seTE.hakn.adhoc.pi.w <- NULL
  #      
  res$prediction.subgroup <- FALSE
  #
  res
}

expandmerge <- function(x, y, nc1 = 0, nr1 = 0, nc2 = 0, nr2 = 0) {
  n1 <- nc1 + nr1
  n2 <- nc2 + nr2
  #
  if ((is.null(x) & is.null(y)))
    return(NULL)
  else if (length(x) == n1 & length(y) == n2)
    return(c(x[seq_len(nc1)], y[seq_len(nc2)],
             x[nc1 + seq_len(nr1)], y[nc2 + seq_len(nr2)]))
  #
  if (length(x) == 1 & n1 > 1)
    x <- rep(x, n1)
  if (length(y) == 1 & n2 > 1)
    y <- rep(y, n2)
  #
  if (length(x) == 0) {
    if (is.character(y))
      x <- rep_len("", n1)
    else
      x <- rep_len(NA, n1)
  }
  if (length(y) == 0) {
    if (is.character(x))
      y <- rep_len("", n2)
    else
      y <- rep_len(NA, n2)
  }
  #
  res <- c(x[seq_len(nc1)], y[seq_len(nc2)],
           x[nc1 + seq_len(nr1)], y[nc2 + seq_len(nr2)])
  res
}

normalize_overall_results <- function(x) {
  for (model in c("common", "random")) {
    lower <- paste0("lower.", model)
    upper <- paste0("upper.", model)
    #
    n.result <- length(x[[lower]])
    #
    if (n.result == 0L)
      next
    #
    if (length(x[[upper]]) != n.result)
      stop("Overall lower and upper confidence limits have different lengths.",
           call. = FALSE)
    #
    vars <-
      paste0(c("TE", "seTE", "statistic", "pval", "text"), ".", model)
    #
    character.vars <- paste0("text.", model)
    #
    if (model == "random") {
      vars <- c(vars, "df.random", "method.random.ci", "adhoc.hakn.ci")
      character.vars <- c(character.vars,
                            "method.random.ci", "adhoc.hakn.ci")
    }
    #
    for (i in vars) {
      var.i <- x[[i]]
      #
      if (is.null(var.i))
        x[[i]] <-
          rep(if (i %in% character.vars) NA_character_ else NA_real_, n.result)
      else if (length(var.i) == 1L && n.result > 1L)
        x[[i]] <- rep(var.i, n.result)
      else if (length(var.i) != n.result)
        stop("Length of '", i,
             "' does not match the number of overall confidence intervals.",
             call. = FALSE)
    }
  }
  #
  x
}

result_characteristics <- function(x) {
  vars <- c("k", "n", "n.e", "n.c", "n.harmonic.mean", "t.harmonic.mean")
  #
  empty <-
    as.data.frame(setNames(rep(list(numeric(0)), length(vars)), vars))
  #
  if (inherits(x, "metamerge")) {
    res <- x[["result.characteristics"]]
    #
    if (is.null(res))
      stop("Result characteristics are missing from the merged object.",
           call. = FALSE)
    #
    for (model in c("common", "random", "predict")) {
      nr <- length(x[[paste0("lower.", model)]])
      #
      if (nr == 0L)
        res[[model]] <- list(overall = empty, subgroup = NULL)
      else if (nrow(res[[model]]$overall) != nr)
        stop("Result characteristics do not match the stored results.",
             call. = FALSE)
    }
    #
    return(res)
  }
  #
  total <- function(x) {
    if (is.null(x)) NULL else sum(x)
  }
  #
  expand <- function(x, nr) {
    if (nr == 0L)
      numeric(0)
    else if (is.null(x))
      rep(NA_real_, nr)
    else if (length(x) == 1L)
      rep(unname(x), nr)
    else if (length(x) == nr)
      unname(x)
    else
      stop("Characteristic does not match the number of results.",
           call. = FALSE)
  }
  #
  overall <- list(
    k = x[["k"]],
    n = if (!is.null(x[["n"]]))
      total(x[["n"]])
    else if (!is.null(x[["n.e"]]) && !is.null(x[["n.c"]]))
      total(x[["n.e"]] + x[["n.c"]])
    else
      NULL,
    n.e = total(x[["n.e"]]),
    n.c = total(x[["n.c"]]),
    #
    n.harmonic.mean = x[["n.harmonic.mean"]],
    t.harmonic.mean = x[["t.harmonic.mean"]])
  #
  subgroup <- list(
    k = x[["k.w"]],
    n = if (!is.null(x[["n.w"]]))
      x[["n.w"]]
    else if (!is.null(x[["n.e.w"]]) && !is.null(x[["n.c.w"]]))
      x[["n.e.w"]] + x[["n.c.w"]]
    else
      NULL,
    n.e = x[["n.e.w"]],
    n.c = x[["n.c.w"]],
    n.harmonic.mean = x[["n.harmonic.mean.w"]],
    t.harmonic.mean = x[["t.harmonic.mean.w"]])
  #
  res <- list()
  #
  for (model in c("common", "random", "predict")) {
    nr <- length(x[[paste0("lower.", model)]])
    #
    if (nr == 0L) {
      res[[model]] <- list(overall = empty, subgroup = NULL)
      next
    }
    #
    overall.model <- as.data.frame(lapply(overall, expand, nr = nr))
    reference <- list2mat(x[[paste0("lower.", model, ".w")]])
    #
    subgroup.model <- NULL
    #
    if (!is.null(reference)) {
      subgroup.model <- lapply(subgroup, function(x) {
        if (length(x) == 1L && is.na(x))
          x <- NULL
        #
        subgroup_result_matrix(x, reference)
      })
    }
    #
    res[[model]] <- list(
      overall = overall.model,
      subgroup = subgroup.model)
  }
  #
  res
}

merge_result_characteristics <- function(x, y) {
  first <- result_characteristics(x)
  second <- result_characteristics(y)
  #
  res <- list()
  #
  for (model in c("common", "random", "predict")) {
    overall <- rbind(first[[model]]$overall, second[[model]]$overall)
    rownames(overall) <- NULL
    #
    a <- first[[model]]$subgroup
    b <- second[[model]]$subgroup
    subgroup <- NULL
    #
    if (!is.null(a) || !is.null(b)) {
      vars <- names(overall)
      subgroup <- setNames(lapply(vars, function(field) {
        left <- a[[field]]
        right <- b[[field]]
        #
        if (!is.null(left) && !is.null(right) &&
            !identical(rownames(left), rownames(right)))
          stop("Subgroup characteristics have different row ordering.",
               call. = FALSE)
        #
        cbind(left, right)
      }), vars)
    }
    #
    res[[model]] <- list(overall = overall, subgroup = subgroup)
  }
  #
  res
}

characteristic_rows <- function(x) {
  if (is.null(x))
    return(NULL)
  #
  as.data.frame(lapply(x, function(y) as.vector(t(y))))
}
