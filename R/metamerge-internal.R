samedata <- function(x, y) {
  if (!is.null(x$data) & !is.null(y$data)) {
    ##
    if (nrow(x$data) != nrow(y$data))
      stop("Meta-analyses based on different data sets.",
           call. = FALSE)
    ##
    mismatch <- FALSE
    ##
    if (inherits(x, "metabin")) {
      mismatch1 <-
        mismatch(x, y, ".event.e")
      mismatch2 <-
        mismatch(x, y, ".n.e")
      mismatch3 <-
        mismatch(x, y, ".event.c")
      mismatch4 <-
        mismatch(x, y, ".n.c")
      ##
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
      ##
      mismatch <-
        mismatch1 | mismatch2 | mismatch3 | mismatch4 |
        mismatch5 | mismatch6
    }
    else if (inherits(x, "metacor")) {
      mismatch1 <-
        mismatch(x, y, ".cor")
      mismatch2 <-
        mismatch(x, y, ".n")
      ##
      mismatch <- mismatch1 | mismatch2
    }
    else if (inherits(x, "metagen")) {
      ##
      mismatch1 <-
        mismatch(x, y, ".n.e")
      mismatch2 <-
        mismatch(x, y, ".n.c")
      ##
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
      ##
      mismatch <- mismatch1 | mismatch2 | mismatch3 | mismatch4
    }
    else if (inherits(x, "metamean")) {
      mismatch1 <-
        mismatch(x, y, ".n")
      mismatch2 <-
        mismatch(x, y, ".mean")
      mismatch3 <-
        mismatch(x, y, ".sd")
      ##
      mismatch <-
        mismatch1 | mismatch2 | mismatch3
    }
    else if (inherits(x, "metaprop")) {
      mismatch1 <-
        mismatch(x, y, ".event")
      mismatch2 <-
        mismatch(x, y, ".n")
      ##
      mismatch <- mismatch1 | mismatch2
    }
    else if (inherits(x, "metarate")) {
      mismatch1 <-
        mismatch(x, y, ".time")
      mismatch2 <-
        mismatch(x, y, ".n")
      ##
      mismatch <- mismatch1 | mismatch2
    }
    ##
    if (mismatch)
      stop("Meta-analyses based on different data sets.",
           call. = FALSE)
  }
  ##
  invisible(NULL)
}


samesm <- function(x, y) {
  if (!is.null(x$sm) & !is.null(y$sm)) {
    if (is_relative_effect(x$sm) != is_relative_effect(y$sm))
      stop("Summary measures used in meta-analyses do not fit.",
           call. = FALSE)
  }
  ##
  if (inherits(x, "metabin")) {
    if ((x$sm != y$sm) &
        any(c(x$sm, y$sm) %in% c("RD", "ASD")))
      stop("Summary measures used in meta-analyses do not fit.",
           call. = FALSE)
  }
  ##
  invisible(NULL)
}


samesubgroups <- function(x, y) {
  if (!is.null(x$subgroup.levels) &
      !is.null(y$subgroup.levels)) {
    ##
    if (length(x$subgroup.levels) !=
        length(y$subgroup.levels))
      stop("Meta-analyses have different number of subgroups.",
           call. = FALSE)
    ##
    if (any(x$subgroup.levels != y$subgroup.levels))
      stop("Meta-analyses based on different subgroup-analyses.",
           call. = FALSE)
  }
  ##
  invisible(NULL)
}


updateobj <- function(x,
                      text.common, text.random, text.predict,
                      text.w.common, text.w.random,
                      hetlabel, taulabel) {
  
  is.copas <- inherits(x, "copas")
  is.limit <- inherits(x, "limitmeta")
  is.robu <- inherits(x, "robu")
  is.trimfill <- inherits(x, "trimfill")
  ##
  res <- x
  ##
  if (is.null(x$hetlabel) || all(x$hetlabel == ""))
    res$hetlabel <- hetlabel
  ##
  if (all(x$detail.tau == ""))
    res$detail.tau <- taulabel
  ##
  ## Act upon ordinary meta-analysis object
  ##
  if (!(is.copas | is.limit | is.robu | is.trimfill)) {
    if (text.common != "")
      res$text.common <- text.common
    if (text.random != "")
      res$text.random <- text.random
    if (text.predict != "")
      res$text.predict <- text.predict
    ##
    if (text.w.common != "")
      res$text.w.common <- text.w.common
    if (text.w.random != "")
      res$text.w.random <- text.w.random
    ##
    return(res)
  }
  ##
  ## Other objects
  ##
  res$text.predict <- ""
  ##
  if (is.copas | is.limit | is.robu) {
    res$method.random <- "Inverse"
    ##
    res$method.random.ci <- "classic"
    res$adhoc.hakn.ci <- ""
    res$df.random <- Inf
    ##
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
  ##
  if (is.copas) {
    res$TE.random <- res$TE.adjust
    res$seTE.random <- res$seTE.adjust
    res$lower.random <- res$lower.adjust
    res$upper.random <- res$upper.adjust
    res$statistic.random <- res$statistic.adjust
    res$pval.random <- res$pval.adjust
    ##
    res$w.random <- rep(0, length(res$w.random))
    ##
    res$tau <- res$tau.adjust
    res$lower.tau <- NA
    res$upper.tau <- NA
    res$tau2 <- res$tau.adjust^2
    res$lower.tau2 <- NA
    res$upper.tau2 <- NA
    res$se.tau <- NA
    ##
    res$method.tau <- "ML"
    ##
    if (hetlabel == "")
      res$hetlabel <- "copas"
    if (taulabel == "")
      res$detail.tau <- "copas"
    ##
    if (text.random != "")
      res$text.random <- text.random
    else
      res$text.random <- "Copas selection model"
    ##
    if (text.w.random != "")
      res$text.w.random <- text.w.random
    else
      res$text.w.random <- "Copas"
  }
  else if (is.limit) {
    res$TE.random <- res$TE.adjust
    res$seTE.random <- res$seTE.adjust
    res$lower.random <- res$lower.adjust
    res$upper.random <- res$upper.adjust
    res$statistic.random <- res$statistic.adjust
    res$pval.random <- res$pval.adjust
    ##
    res$w.random <- rep(0, length(res$w.random))
    ##
    if (hetlabel == "")
      res$hetlabel <- "limit"
    if (taulabel == "")
      res$detail.tau <- "limit"
    ##
    if (text.random != "")
      res$text.random <- text.random
    else
      res$text.random <- "Limit meta-analysis"
    ##
    if (text.w.random != "")
      res$text.w.random <- text.w.random
    else
      res$text.w.random <- "limit"
  }
  else if (is.robu) {
    res$TE.random <- res$reg_table$b.r[1]
    res$seTE.random <- res$reg_table$SE[1]
    res$lower.random <- res$reg_table$CI.L[1]
    res$upper.random <- res$reg_table$CI.U[1]
    res$statistic.random <- res$reg_table$t[1]
    res$pval.random <- res$reg_table$prob[1]
    ##
    res$w.random <- res$data.full$r.weights
    ##
    res$level.ma <- 0.95
    res$tau <- sqrt(res$mod_info$tau.sq)
    res$lower.tau <- NA
    res$upper.tau <- NA
    res$tau2 <- res$mod_info$tau.sq
    res$lower.tau2 <- NA
    res$upper.tau2 <- NA
    res$se.tau <- NA
    ##
    res$method.tau <- "DL"
    ##
    if (hetlabel == "")
      res$hetlabel <- "RVE"
    if (taulabel == "")
      res$detail.tau <- "RVE"
    ##
    if (text.random != "")
      res$text.random <- text.random
    else
      res$text.random <- "RVE model"
    ##
    if (text.w.random != "")
      res$text.w.random <- text.w.random
    else
      res$text.w.random <- "RVE"
  }
  else if (is.trimfill) {
    if (hetlabel == "")
      res$hetlabel <- "TF"
    if (taulabel == "")
      res$detail.tau <- "TF"
    ##
    if (text.common != "")
      res$text.common <- text.common
    else
      res$text.common <- "Trim-and-fill method (CE)"
    ##
    if (text.random != "")
      res$text.random <- text.random
    else
      res$text.random <- "Trim-and-fill method (RE)"
    ##
    if (text.w.common != "")
      res$text.w.common <- text.w.common
    else
      res$text.w.common <- "trim-fill"
    ##
    if (text.w.random != "")
      res$text.w.random <- text.w.random
    else
      res$text.w.random <- "trim-fill"
  }
  
  res
}

dropcommon <- function(x, dropsubgroup) {
  res <- x
  ##
  res$method <- NULL
  ##
  res$w.common <- NULL
  res$TE.common <- NULL
  res$seTE.common <- NULL
  res$statistic.common <- NULL
  res$pval.common <- NULL
  res$lower.common <- NULL
  res$upper.common <- NULL
  res$zval.common <- NULL
  res$text.common <- NULL
  ##
  res$k.MH <- NULL
  ##
  if (missing(dropsubgroup))
    dropsubgroup <- !is.null(res$subgroup)
  if (dropsubgroup) {
    res$TE.common.w <- NULL
    res$seTE.common.w <- NULL
    res$statistic.common.w <- NULL
    res$pval.common.w <- NULL
    res$lower.common.w <- NULL
    res$upper.common.w <- NULL
    res$w.common.w <- NULL
    ##
    res$Q.w.common <- NULL
    res$pval.Q.w.common <- NULL
    ##
    res$Q.b.common <- NULL
    res$df.Q.b.common <- NULL
    res$pval.Q.b.common <- NULL
    ##
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
  }
  ##
  res
}

droprandom <- function(x, dropsubgroup) {
  res <- x
  ##
  res$method.random <- NULL
  ##
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
  ##
  res$seTE.classic <- NULL
  ##
  res$adhoc.hakn.ci <- NULL
  res$df.hakn <- NULL
  res$seTE.hakn.ci <- NULL
  res$seTE.hakn.adhoc.ci <- NULL
  ##
  res$df.kero <- NULL
  res$seTE.kero <- NULL
  ##
  res$text.random <- NULL
  ##
  res$cluster <- FALSE
  res$three.level <- FALSE
  ##
  res$k.study <- NULL
  ##
  if (missing(dropsubgroup))
    dropsubgroup <- !is.null(res$subgroup)
  if (dropsubgroup) {
    res$TE.random.w <- NULL
    res$seTE.random.w <- NULL
    res$statistic.random.w <- NULL
    res$pval.random.w <- NULL
    res$df.random.w <- NULL
    res$lower.random.w <- NULL
    res$upper.random.w <- NULL
    res$w.random.w <- NULL
    ##
    res$seTE.classic.w <- NULL
    ##
    res$df.hakn.ci.w <- NULL
    res$seTE.hakn.ci.w <- NULL
    res$seTE.hakn.adhoc.ci.w <- NULL
    ##
    res$df.kero.w <- NULL
    res$seTE.kero.w <- NULL
    ##
    res$Q.w.random <- NULL
    res$pval.Q.w.random <- NULL
    ##
    res$Q.b.random <- NULL
    res$df.Q.b.random <- NULL
    res$pval.Q.b.random <- NULL
    ##
    res$zval.random.w <- NULL
  }
  ##
  res
}

droppredict <- function(x) {
  res <- x
  ##
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
  ##      
  if (!is.null(res$subgroup)) {
    res$seTE.predict.w <- NULL
    res$df.predict.w <- NULL
    res$lower.predict.w <- NULL
    res$upper.predict.w <- NULL
    res$seTE.hakn.pi.w <- NULL
    res$seTE.hakn.adhoc.pi.w <- NULL
  }
  ##
  res
}


mergevars <- function(x, y, name.x = NULL, name.y = NULL,
                      replace = NULL) {
  n1 <- length(x)
  n2 <- length(y)
  ##
  if (!is.null(name.x) & !is.null(name.y)) {
    if (is.null(names(x)))
      names(x) <- rep(name.x, n1)
    else if (any(name.x != "") & any(name.x != names(x)))
      names(x) <- paste(name.x, names(x), sep = "-")
    ##
    if (is.null(names(y)))
      names(y) <- rep(name.y, n2)
    else if (any(name.y != "") & any(name.y != names(y)))
      names(y) <- paste(name.y, names(y), sep = "-")
  }
  if (!is.null(replace)) {
    x <- replaceNULL(x, rep(replace, n1))
    y <- replaceNULL(y, rep(replace, n2))
  }
  ##
  if (n1 > 1 | n2 > 1)
    res <- list(x, y)
  else
    res <- c(x, y)
  ##
  res
}


expandmerge <- function(x, y, nc1 = 0, nr1 = 0, nc2 = 0, nr2 = 0) {
  n1 <- nc1 + nr1
  n2 <- nc2 + nr2
  ##
  if ((is.null(x) & is.null(y)))
    return(NULL)
  else if (length(x) == n1 & length(y) == n2)
    return(c(x[seq_len(nc1)], y[seq_len(nc2)],
             x[nc1 + seq_len(nr1)], y[nc2 + seq_len(nr2)]))
  ##
  if (length(x) == 1 & n1 > 1)
    x <- rep(x, n1)
  if (length(y) == 1 & n2 > 1)
    y <- rep(y, n2)
  ##
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
  ##
  res <- c(x[seq_len(nc1)], y[seq_len(nc2)],
           x[nc1 + seq_len(nr1)], y[nc2 + seq_len(nr2)])
  res
}
