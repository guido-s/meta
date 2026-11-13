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
      func = setClass(x),
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
      res$args.transf <- list(x$args.transf)
      res$func.backtransf <- replaceNULL(x$func.backtransf)
      res$args.backtransf <- list(x$args.backtransf)
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
    res$n.harmonic.mean <- replaceNULL(x[["n.harmonic.mean"]])
    res$t.harmonic.mean <- replaceNULL(x[["t.harmonic.mean"]])
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
    res$tau2 <- list(x$tau2)
    res$lower.tau2 <- list(x$lower.tau2)
    res$upper.tau2 <- list(x$upper.tau2)
    res$se.tau2 <- list(replaceNULL(x$se.tau2))
    res$tau <- list(x$tau)
    res$lower.tau <- list(x$lower.tau)
    res$upper.tau <- list(x$upper.tau)
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
  if (random || type == "prediction")
    res$interval.slot <- seq_len(nrow(res))
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
  if (type != "prediction") {
    n.method <- if (type == "common") 1L else length(x$method.random.ci)
    Q <- replaceNULL(x[[paste0("Q.b.", type)]], numeric())
    n.test <- length(Q)
    if (n.test > n.method)
      stop("More subgroup tests than interval methods.", call. = FALSE)
    #
    Q.by.method <- pval.by.method <- rep(NA_real_, n.method)
    df.by.method <- rep(list(NA_real_), n.method)
    #
    if (n.test > 0L) {
      dfs <- subgroup_test_df(
        replaceNULL(x[[paste0("df.Q.b.", type)]], x[["df.Q.b"]]), n.test)
      pvals <- replaceNULL(x[[paste0("pval.Q.b.", type)]], pval_subgroup(Q, dfs))
      slots <- seq_len(n.test)
      Q.by.method[slots] <- Q
      df.by.method[slots] <- dfs
      pval.by.method[slots] <- pvals[slots]
    }
  }
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
      Q.b = Q.by.method,
      df.Q.b = NA_real_,
      pval.Q.b = pval.by.method
    )
  }
  else if (type == "random") {
    res <- data.frame(
      id = id,
      #
      interval.slot = rep(seq_len(n.method), each = n.by),
      method.interval = rep(unname(x$method.random.ci), each = n.by),
      adhoc.hakn = rep(unname(x$adhoc.hakn.ci), each = n.by),
      #
      subgroup.name = x$subgroup.name,
      subgroup.levels = rep(x$subgroup.levels, times = n.method),
      k.w = rep(unname(x$k.w), times = n.method),
      #
      TE = matrix(list2mat(x$TE.random.w), ncol = 1),
      seTE = matrix(list2mat(x$seTE.random.w), ncol = 1),
      lower = matrix(list2mat(x$lower.random.w), ncol = 1),
      upper = matrix(list2mat(x$upper.random.w), ncol = 1),
      df = matrix(list2mat(x$df.random.w), ncol = 1),
      statistic = matrix(list2mat(x$statistic.random.w), ncol = 1),
      pval = matrix(list2mat(x$pval.random.w), ncol = 1),
      #
      Q.b = rep(Q.by.method, each = n.by),
      df.Q.b = NA_real_,
      pval.Q.b = rep(pval.by.method, each = n.by)
    )
  }
  else if (type == "prediction") {
    n.method <- length(x$method.predict)
    #
    res <- data.frame(
      id = id,
      #
      interval.slot = rep(seq_len(n.method), each = n.by),
      method.interval = rep(unname(x$method.predict), each = n.by),
      adhoc.hakn = rep(unname(x$adhoc.hakn.pi), each = n.by),
      #
      subgroup.name = x$subgroup.name,
      subgroup.levels = rep(x$subgroup.levels, times = n.method),
      k.w = rep(unname(x$k.w), times = n.method),
      #
      lower = matrix(list2mat(x$lower.predict.w), ncol = 1),
      upper = matrix(list2mat(x$upper.predict.w), ncol = 1),
      df = matrix(list2mat(x$df.predict.w), ncol = 1)
    )
  }
  #
  if (type != "prediction") {
    df.rows <- rep(df.by.method, each = n.by)
    #
    if (length(df.rows) != nrow(res))
      stop("Subgroup-test DF entries do not match the result rows.",
           call. = FALSE)
    #
    res$df.Q.b <- df.rows
    #
    res$test.Q.b <- vapply(df.rows, function(df) {
      if (length(df) == 0L || anyNA(df))
        NA_character_
      else if (length(df) == 1L)
        "Q"
      else if (length(df) == 2L)
        "F"
      else
        NA_character_
    }, character(1))
  }
  #
  res$n.harmonic.mean <- replaceNULL(x$n.harmonic.mean.w)
  res$t.harmonic.mean <- replaceNULL(x$t.harmonic.mean.w)
  #
  res
}

subgroup_test_df <- function(df, n.method) {
  if (is.null(df))
    df <- rep(list(NA_real_), n.method)
  else if (!is.list(df)) {
    if (n.method == 1L)
      df <- list(df)
    else if (length(df) == 1L)
      df <- rep(list(df), n.method)
    else
      stop("Multiple subgroup tests require a list of ",
           "method-specific degrees of freedom.",
           call. = FALSE)
  }
  #
  if (length(df) != n.method)
    stop("Number of subgroup-test DF entries does not match ",
         "the number of interval methods.",
         call. = FALSE)
  #
  unname(df)
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
            class[1],
            call. = FALSE)
    #
    class <- class[1]
  }
  #
  class
}

format_args <- function(x, lab.NA = ".") {
  if (is.null(x) || length(x) == 0)
    return(lab.NA)
  #
  values <- vapply(x, function(value)
    paste(deparse(value, width.cutoff = 500L), collapse = " "), character(1))
  #
  arg.names <- names(x)
  #
  if (!is.null(arg.names)) {
    named <- !is.na(arg.names) & nzchar(arg.names)
    values[named] <- paste0(arg.names[named], " = ", values[named])
  }
  paste(values, collapse = ", ")
}

metabind_interval_occurrence <- function(dat, values) {
  occurrence <- integer(nrow(dat))
  ids <- replaceNULL(dat$id, rep(1L, nrow(dat)))
  # Shared summaries acquire missing IDs when bound to source subgroup rows.
  ids <- match(ids, unique(ids))
  #
  for (id in unique(ids)) {
    rows <- which(ids == id)
    if (isCol(dat, "interval.slot"))
      rows <- rows[order(dat$interval.slot[rows])]
    #
    while (length(rows)) {
      first <- !duplicated(values[rows, , drop = FALSE])
      occurrence[rows] <- occurrence[rows] + 1L
      rows <- rows[!first]
    }
  }
  #
  occurrence
}

map_metabind_subgroups <- function(ma, x, dat, display, shared.overall,
                                   common, random, prediction, overall,
                                   shared.heterogeneity) {
  groups <- as.character(ma$subgroup.levels)
  n.by <- length(groups)
  metadata <- match(groups, as.character(x$overall$id))
  #
  args <- list(common = common, random = random,
               prediction = prediction && shared.overall,
               overall = overall && shared.overall,
               subgroup = rep(TRUE, n.by),
               subgroup.hetstat = !shared.overall, resid.hetstat = FALSE,
               test.overall = shared.overall, test.subgroup = FALSE,
               test.effect.subgroup = !shared.overall,
               text = list(subgroup = x$overall$name[metadata]))
  #
  global.rows <- group.rows <- list()
  #
  families <- c(common = "Common effect model",
                random = "Random effects model",
                predict = "Prediction interval")
  #
  selected <- c(common = common, random = random, predict = prediction)
  #
  fields <- c("TE", "seTE", "lower", "upper", "statistic", "pval", "df")
  #
  for (family in names(families)) {
    rows <- which(dat$type == families[[family]] & selected[[family]])
    keys <- intersect(c("method.interval", "adhoc.hakn",
                        if (family == "predict")
                          "level.predict"
                        else
                          "level.ma",
                        "text"), names(dat))
    #
    slots <- dat[rows, keys, drop = FALSE]
    #
    if (family == "common")
      slots <- slots[seq_len(min(nrow(slots), 1L)), , drop = FALSE]
    else {
      occurrence <- metabind_interval_occurrence(dat[rows, , drop = FALSE], slots)
      slots$.interval.occurrence <- occurrence
      slots <- unique(slots)
    }
    #
    ns <- max(nrow(slots), 1L)
    #
    slot <- integer(length(rows))
    #
    if (family == "common")
      slot[] <- 1L
    else {
      for (j in seq_len(nrow(slots))) {
        same <- occurrence == slots$.interval.occurrence[j]
        for (key in keys) {
          value <- dat[[key]][rows]
          target <- slots[[key]][j]
          same <- same & ((is.na(value) & is.na(target)) |
                            (!is.na(value) & !is.na(target) & value == target))
        }
        slot[same] <- j
      }
    }
    #
    global <- rep(NA_integer_, ns)
    grouped <- matrix(NA_integer_, n.by, ns,
                      dimnames = list(groups, NULL))
    #
    if (shared.overall) {
      global[slot] <- rows
    }
    else if (length(rows)) {
      group <- match(as.character(dat$id[rows]), groups)
      keep <- !is.na(group)
      grouped[cbind(group[keep], slot[keep])] <- rows[keep]
    }
    #
    labels <- if (nrow(slots) && "text" %in% keys)
      as.character(slots$text) else rep(NA_character_, ns)
    #
    fallback <- switch(family,
                       common = "Common effect model",
                       random = "Random effects model",
                       predict = "Prediction interval")
    #
    labels[is.na(labels) | !nzchar(labels)] <- fallback
    #
    repeated <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
    #
    if (any(repeated) && nrow(slots)) {
      codes <-
        vapply(seq_len(nrow(slots)),
               function(j) {
                 values <-
                   unlist(slots[j, setdiff(keys, "text"), drop = FALSE],
                          use.names = FALSE)
                 paste(values[!is.na(values) & nzchar(as.character(values))],
                       collapse = ", ")
               }, character(1))
      labels[repeated] <- paste0(labels[repeated], " (", codes[repeated], ")")
    }
    #
    colnames(grouped) <- labels
    #
    global.rows[[family]] <- global
    group.rows[[family]] <- grouped
    #
    args$text[[family]] <- labels
    args$text[[paste0(family, ".w")]] <- labels
    #
    ma[[paste0("text.", family)]] <- labels
    #
    for (field in fields) {
      values <- dat[[field]]
      if (is.null(values))
        values <- rep(NA_real_, nrow(dat))
      #
      ma[[paste0(field, ".", family)]] <- unname(values[global])
      ma[[paste0(field, ".", family, ".w")]] <-
        matrix(values[grouped], n.by, ns, dimnames = dimnames(grouped))
    }
    #
    method <- if (nrow(slots) && "method.interval" %in% keys)
      as.character(slots$method.interval) else "classic"
    #
    adhoc <- if (nrow(slots) && "adhoc.hakn" %in% keys)
      as.character(slots$adhoc.hakn) else ""
    adhoc[is.na(adhoc)] <- ""
    #
    method.field <- switch(family,
                           common = "method.common.ci",
                           random = "method.random.ci",
                           predict = "method.predict")
    #
    ma[[method.field]] <- method
    #
    if (family == "predict")
      level.field <- "level.predict"
    else
      level.field <- "level.ma"
    #
    if (nrow(slots) && level.field %in% keys)
      ma[[level.field]] <- makeunique(slots[[level.field]], gs(level.field))
    #
    if (length(rows) && "method" %in% names(dat) &&
        (family != "predict" || !random)) {
      name <- if (family == "common") "method" else "method.random"
      #
      ma[[name]] <- dat$method[rows[1]]
      #
      ma[[paste0(name, ".w")]] <-
        matrix(dat$method[grouped], n.by, ns, dimnames = dimnames(grouped))
    }
    if (family != "common") {
      suffix <- if (family == "random") "ci" else "pi"
      ma[[paste0("adhoc.hakn.", suffix)]] <- adhoc
      #
      for (field in c("df.hakn", "df.kero")) {
        values <- dat[[field]]
        #
        if (is.null(values))
          values <- rep(NA_real_, nrow(dat))
        #
        name <- if (field == "df.hakn") paste0(field, ".", suffix) else field
        # KR degrees of freedom belong to the random CI, not a PI
        if (field == "df.kero" && family == "predict")
          next
        #
        ma[[name]] <- unname(values[global])
        #
        ma[[paste0(name, ".w")]] <-
          matrix(values[grouped], n.by, ns, dimnames = dimnames(grouped))
      }
    }
    #
    mask <- !is.na(grouped)
    #
    flag <- switch(family,
                   common = "common.subgroup",
                   random = "random.subgroup",
                   predict = "prediction.subgroup")
    #
    args[[flag]] <- if (family == "common") unname(mask[, 1]) else mask
    #
    if (family == "common" && !shared.overall && "text" %in% names(dat)) {
      text <- as.character(dat$text[grouped[, 1]])
      text[is.na(text) | !nzchar(text)] <- fallback
      args$text$common.w <- text
    }
  }
  #
  null.effect <- replaceNULL(x$overall[["null.effect"]][metadata], rep(NA_real_, n.by))
  has.null <- !is.na(null.effect)
  args$test.effect.subgroup.common <- args$common.subgroup & has.null
  args$test.effect.subgroup.random <- unname(rowSums(args$random.subgroup) > 0) & has.null
  #
  args$label <- list(
    test.effect.subgroup.common =
      paste0("Test for overall effect",
             if (common && random)
               paste0(" (", gs("text.w.common"), " effect)"),
             ": "),
    test.effect.subgroup.random =
      paste0("Test for overall effect",
             if (common && random)
               paste0(" (", gs("text.w.random"), " effects)"),
             ": "))
  #
  heterogeneity <- c("Q", "df.Q", "pval.Q", "tau2", "lower.tau2", "upper.tau2",
                     "se.tau2", "tau", "lower.tau", "upper.tau", "H", "lower.H",
                     "upper.H", "I2", "lower.I2", "upper.I2", "Rb", "lower.Rb",
                     "upper.Rb")
  #
  settings <- c("method.I2", "method.tau", "method.tau.ci", "level.hetstat")
  #
  for (field in c(heterogeneity, settings)) {
    values <- x$overall[[field]]
    if (is.null(values))
      next
    #
    if (is.list(values)) {
      value <- values[[metadata[1]]]
      # Native forest annotations use total variance, not component-wise CIs.
      #
      group.values <- vapply(values[metadata], function(value) {
        if (length(value) == 1L)
          value
        else if (field == "tau2")
          sum(value)
        else if (field == "tau")
          sqrt(sum(value^2))
        else
          NA_real_
      }, numeric(1))
    }
    else {
      value <- values[metadata[1]]
      group.values <- values[metadata]
    }
    #
    if (shared.overall && shared.heterogeneity)
      ma[[field]] <- value
    else if (field %in% heterogeneity)
      ma[[field]] <- NA_real_
    else
      ma[[field]] <- value
    #
    if (shared.overall && shared.heterogeneity)
      ma[[paste0(field, ".w")]] <-
        rep(if (is.numeric(value)) NA_real_ else NA_character_, n.by)
    else
      ma[[paste0(field, ".w")]] <- group.values
    #
    names(ma[[paste0(field, ".w")]]) <- groups
  }
  #
  for (field in c("Q.b", "df.Q.b", "pval.Q.b")) {
    ma[[field]] <- NA_real_
    #
    for (family in c("common", "random"))
      ma[[paste0(field, ".", family)]] <-
        rep(NA_real_, length(global.rows[[family]]))
  }
  #
  ma$common <- common
  ma$random <- random
  ma$prediction <- args$prediction
  ma$overall <- args$overall
  ma$three.level <- any(x$overall[["three.level"]][metadata], na.rm = TRUE)
  ma$k.w.orig <- x$overall$k[metadata]
  #
  if (shared.overall)
    ma$null.effect <- null.effect[1]
  #
  ma$common.subgroup <- args$common.subgroup
  ma$random.subgroup <- args$random.subgroup
  ma$prediction.subgroup <- args$prediction.subgroup
  #
  index <- c(unlist(global.rows, use.names = FALSE),
             unlist(lapply(group.rows, function(rows) as.vector(t(rows))),
                    use.names = FALSE))
  #
  pooled <- display[index, , drop = FALSE]
  #
  if ("seTE.plot" %in% names(ma$data) && !"seTE.plot" %in% names(pooled) &&
      "seTE" %in% names(pooled))
    pooled$seTE.plot <- pooled$seTE
  #
  for (field in setdiff(names(ma$data), names(pooled)))
    pooled[[field]] <- ma$data[[field]][rep(NA_integer_, length(index))]
  #
  rownames(pooled) <- NULL
  #
  args$data.pooled <- pooled
  #
  list(ma = ma, args = args)
}

format_metabind_variance <- function(x, digits, lab.NA = ".", big.mark = "",
                                     upper = NULL) {
  if (is.null(upper))
    vapply(x, function(value) {
      if (length(value) == 0L)
        return(lab.NA)
      #
      paste(trimws(formatPT(value, digits = digits, lab.NA = lab.NA,
                            big.mark = big.mark)),
            collapse = ", ")
    }, character(1), USE.NAMES = FALSE)
  else
    mapply(function(lower, upper) {
      if (length(lower) == 0L || length(upper) == 0L)
        return(lab.NA)
      #
      paste(trimws(formatCI(lower, upper, digits = digits, lab.NA = lab.NA,
                            na.lab = TRUE, big.mark = big.mark)),
            collapse = ", ")
    }, x, upper, USE.NAMES = FALSE)
}

format_metabind_heterogeneity <- function(dat, print, digits,
                                         ids = NULL, pooled = FALSE,
                                         lab.NA = ".", big.mark = gs("big.mark"),
                                         scientific.pval = gs("scientific.pval"),
                                         zero.pval = gs("zero.pval"),
                                         JAMA.pval = gs("JAMA.pval"),
                                         text.tau2 = gs("text.tau2"),
                                         text.tau = gs("text.tau"),
                                         text.I2 = gs("text.I2")) {
  if (!is.null(ids))
    dat <- dat[dat$id %in% ids, , drop = FALSE]
  #
  columns <- c(if (!pooled) "name", "func", "sm",
               "method.tau", "method.tau.ci", "method.I2", "level.hetstat",
               "Q", "df.Q", "pval.Q", "tau2", "lower.tau2", "upper.tau2",
               "se.tau2", "tau", "lower.tau", "upper.tau",
               "H", "lower.H", "upper.H", "I2", "lower.I2", "upper.I2")
  dat <- dat[, intersect(columns, names(dat)), drop = FALSE]
  #
  if (pooled)
    dat <- unique(dat)
  #
  if (isCol(dat, "Q"))
    dat$Q <- formatN(dat$Q, digits = digits[["Q"]], lab.NA = lab.NA,
                     big.mark = big.mark, keep.allNA = TRUE)
  #
  if (isCol(dat, "pval.Q"))
    dat$pval.Q <- formatPT(dat$pval.Q, digits = digits[["pval.Q"]],
                           lab.NA = lab.NA,
                           scientific = scientific.pval, zero = zero.pval,
                           JAMA = JAMA.pval, keep.allNA = TRUE)
  #
  for (field in c("tau2", "tau")) {
    if (isCol(dat, field))
      dat[[field]] <- format_metabind_variance(
        dat[[field]], digits = digits[[field]],
        lab.NA = lab.NA, big.mark = big.mark)
    #
    lower <- paste0("lower.", field)
    upper <- paste0("upper.", field)
    #
    if (!(allNA(unlist(dat[[lower]])) && allNA(unlist(dat[[upper]]))))
      dat[[lower]] <- format_metabind_variance(
        dat[[lower]], upper = dat[[upper]], digits = digits[[field]],
        lab.NA = lab.NA, big.mark = big.mark)
  }
  #
  if (!allNA(dat$I2)) {
    dat$I2 <-
      percent_label(dat$I2, "", digits[["I2"]], formatN, lab.NA = lab.NA)
  }
  #
  if (!(allNA(dat$lower.I2) && allNA(dat$upper.I2)))
    dat$lower.I2 <- formatCI(dat$lower.I2, dat$upper.I2,
                             digits = digits[["I2"]],
                             scale = 100, suffix = "%", lab.NA = lab.NA,
                             na.lab = TRUE)
  #
  if (isCol(dat, "H"))
    dat$H <- formatN(dat$H, digits = digits[["I2"]], lab.NA = lab.NA,
                     big.mark = big.mark, keep.allNA = TRUE)
  #
  if (!(allNA(dat$lower.H) && allNA(dat$upper.H)))
    dat$lower.H <- formatCI(dat$lower.H, dat$upper.H, digits = digits[["I2"]],
                            lab.NA = lab.NA, na.lab = TRUE,
                            big.mark = big.mark)
  #
  dat$level.hetstat <- round(dat$level.hetstat, 4)
  dat[c("se.tau2", "upper.tau2", "upper.tau", "upper.H", "upper.I2")] <- NULL
  names(dat)[names(dat) == "level.hetstat"] <- "level"
  #
  same.level <- length(unique(dat$level)) == 1 && !anyNA(dat$level)
  show.level <- replaceNULL(print[["level.hetstat"]], !same.level)
  cilab <- if (same.level) percent_label(dat$level[1]) else "CI"
  #
  if (!show.level)
    dat$level <- NULL
  #
  for (field in c("func", "sm", "method.tau", "method.tau.ci", "method.I2")) {
    flag <- if (field %in% c("func", "sm")) field else "method"
    if (isFALSE(print[[flag]]) ||
        (is.null(print[[flag]]) && length(unique(dat[[field]])) == 1))
      dat[[field]] <- NULL
  }
  #
  if (!(print[["tau2"]] || print[["tau"]]))
    dat$method.tau <- NULL
  if (!((print[["tau2"]] && print[["tau2.ci"]]) ||
        (print[["tau"]] && print[["tau.ci"]])))
    dat$method.tau.ci <- NULL
  if (!print[["I2"]])
    dat$method.I2 <- NULL
  #
  if (!print[["Q"]])
    dat$Q <- NULL
  if (print[["Q"]] && print[["pval.Q"]]) {
    names(dat)[names(dat) == "df.Q"] <- "df"
    names(dat)[names(dat) == "pval.Q"] <- "p-value"
  }
  else {
    dat$df.Q <- NULL
    dat$pval.Q <- NULL
  }
  #
  labels <- list(tau2 = text.tau2, tau = text.tau, I2 = text.I2, H = "H")
  for (field in names(labels)) {
    lower <- paste0("lower.", field)
    if (print[[field]])
      names(dat)[names(dat) == field] <- labels[[field]]
    else
      dat[[field]] <- NULL
    #
    if (print[[field]] && print[[paste0(field, ".ci")]])
      names(dat)[names(dat) == lower] <- cilab
    else
      dat[[lower]] <- NULL
  }
  #
  if (!((print[["tau2"]] && print[["tau2.ci"]]) ||
        (print[["tau"]] && print[["tau.ci"]]) ||
        (print[["I2"]] && print[["I2.ci"]]) ||
        (print[["H"]] && print[["H.ci"]])))
    dat$level <- NULL
  #
  for (field in intersect(c("df", "level"), names(dat))) {
    if (isCol(dat, field))
      dat[[field]] <- formatN(dat[[field]],
                              digits = if (field == "level") 4 else 0,
                              lab.NA = lab.NA, big.mark = big.mark,
                              keep.allNA = TRUE)
  }
  #
  drop_allNAs(dat)
}

format_df <- function(df, lab.NA = ".", big.mark = "") {
  if (length(df) == 0L || all(is.na(df)))
    return(lab.NA)
  #
  paste(vapply(df, function(value) {
    digits <- if (is.finite(value) && is_wholenumber(value)) 0 else 2
    #
    trimws(formatN(value, digits = digits,
                   lab.NA = lab.NA, big.mark = big.mark))
  }, character(1)), collapse = ", ")
}

subgroup_test_heading <- function(test) {
  test <- unique(test[test %in% c("Q", "F")])
  #
  if (length(test) == 1L)
    test
  else
    "Q/F"
}

format_metabind_results <- function(dat, model, print, digits, sm.lab, lab.TE,
                                    with.subgroups = FALSE, lab.NA = ".",
                                    big.mark = gs("big.mark"),
                                    scientific.pval = gs("scientific.pval"),
                                    zero.pval = gs("zero.pval"),
                                    JAMA.pval = gs("JAMA.pval"),
                                    keep = character(),
                                    lab.interval = NULL, show.level = NULL) {
  # Get rid of warning 'no visible binding for global variable'
  type <- Q.b <- df.Q.b <- pval.Q.b <- NULL
  #
  model <- setchar(model, c("common", "random", "prediction"), length = 1)
  #
  result.type <- setchar(model, c("Common effect model",
                                  "Random effects model",
                                  "Prediction interval"))
  #
  dat %<>% filter(type == result.type)
  #
  if (with.subgroups && model != "prediction") {
    test.cols <- c("Q.b", "df.Q.b", "pval.Q.b", "test.Q.b")
    #
    for (field in setdiff(test.cols, names(dat)))
      dat[[field]] <-
        rep(if (field == "test.Q.b") NA_character_ else NA_real_, nrow(dat))
    #
    if (print[["pval.Q.b"]]) {
      group.cols <- intersect(
        c("id", "interval.slot", "method.interval", "adhoc.hakn", "level.ma"),
        names(dat))
      #
      dat %<>%
        group_by(across(all_of(group.cols))) %>%
        mutate(
          across(all_of(test.cols),
                 ~ replace(.x, row_number() != 1L, NA))) %>%
        ungroup() %>%
        as.data.frame()
    }
  }
  #
  if (is.null(dat$func.backtransf))
    dat$func.backtransf <- rep(NA_character_, nrow(dat))
  #
  if (is.null(dat$args.backtransf))
    dat$args.backtransf <- rep(list(NULL), nrow(dat))
  #
  has.backtransf <- !is.na(dat$func.backtransf) & dat$func.backtransf != ""
  has.args.backtransf <- lengths(dat$args.backtransf) > 0
  #
  dat$func.backtransf <- as.character(dat$func.backtransf)
  dat$func.backtransf[!has.backtransf] <- lab.NA
  dat$args.backtransf <- vapply(dat$args.backtransf, format_args,
                                character(1), lab.NA = lab.NA)
  #
  if (model == "prediction")
    level <- "level.predict"
  else
    level <- "level.ma"
  #
  columns <- c("id", "name", "subgroup.name", "subgroup.levels", "k",
               "func", "sm", "method", "method.interval",
               if (model != "common") "adhoc.hakn",
               level, "func.backtransf", "args.backtransf",
               if (model != "prediction") "TE",
               "seTE", "lower", "upper",
               if (model != "prediction") c("statistic", "pval"),
               if (model != "common") "df",
                if (model != "prediction") c("Q.b", "df.Q.b", "pval.Q.b", "test.Q.b"))
  #
  dat <- dat[, intersect(union(columns, keep), names(dat)), drop = FALSE]
  #
  if (model != "prediction") {
    if (isCol(dat, "TE"))
      dat$TE <- formatN(dat$TE, digits = digits[["digits"]], lab.NA = lab.NA,
                        big.mark = big.mark)
    #
    if (isCol(dat, "statistic"))
      dat$statistic <- formatN(dat$statistic, digits = digits[["stat"]],
                               lab.NA = lab.NA, keep.allNA = TRUE)
    #
    if (isCol(dat, "pval"))
      dat$pval <- formatPT(dat$pval, digits = digits[["pval"]], lab.NA = lab.NA,
                           scientific = scientific.pval, zero = zero.pval,
                           JAMA = JAMA.pval, keep.allNA = TRUE)
  }
  #
  if (isCol(dat, "seTE"))
    dat$seTE <- formatN(dat$seTE, digits = digits[["se"]], lab.NA = lab.NA,
                        big.mark = big.mark, keep.allNA = TRUE)
  #
  if (isCol(dat, "lower"))
    dat$lower <- formatCI(dat$lower, dat$upper, digits = digits[["digits"]],
                          lab.NA = lab.NA, na.lab = TRUE, big.mark = big.mark)
  dat$upper <- NULL
  #
  if (model != "common" && isCol(dat, "df"))
      dat$df <- round(dat$df, 2)
  #
  if (isCol(dat, level))
    dat[[level]] <- round(dat[[level]], 4)
  #
  if (model != "prediction")
    names(dat)[names(dat) == "TE"] <- sm.lab
  #
  names(dat)[names(dat) == "method.interval"] <- "interval"
  names(dat)[names(dat) == level] <- "level"
  #
  same.level <- length(unique(dat$level)) == 1 && !anyNA(dat$level)
  show.level <-
    replaceNULL(show.level, replaceNULL(print[[level]], !same.level))
  #
  if (!is.null(lab.interval))
    names(dat)[names(dat) == "lower"] <- lab.interval
  else if (same.level)
    names(dat)[names(dat) == "lower"] <-
    percent_label(dat$level[1],
                  type = if (model == "prediction") "PI" else "CI")
  else {
    names(dat)[names(dat) == "lower"] <-
      if (model == "prediction") "PI" else "CI"
  }
  #
  if (!show.level)
    dat$level <- NULL
  #
  show.backtransf <- replaceNULL(print[["backtransf"]], any(has.backtransf))
  show.args.backtransf <- show.backtransf &&
    replaceNULL(print[["args.backtransf"]], any(has.args.backtransf))
  #
  if (show.backtransf)
    names(dat)[names(dat) == "func.backtransf"] <- "backtransf"
  else
    dat$func.backtransf <- NULL
  #
  if (show.args.backtransf)
    names(dat)[names(dat) == "args.backtransf"] <- "args"
  else
    dat$args.backtransf <- NULL
  #
  interval <- dat$interval
  #
  for (field in c("func", "sm", "method", "method.interval", "k")) {
    if (field == "method.interval")
      column <- "interval"
    else
      column <- field
    #
    if (model == "prediction" && field == "method")
      constant <- all(dat[[column]] %in% "Inverse")
    else
      constant <- length(unique(dat[[column]])) == 1
    #
    if (isFALSE(print[[field]]) || (is.null(print[[field]]) && constant))
      dat[[column]] <- NULL
  }
  #
  if (model != "common") {
    if (all(is.infinite(dat$df)))
      dat$df <- NULL
    #
    if (model == "random")
      methods <- "HK"
    else
      methods <- c("HK", "HK-PR")
    #
    if (is.null(dat$interval) || !any(dat$interval %in% methods) ||
        all(is.na(dat$adhoc.hakn) | dat$adhoc.hakn == ""))
      dat$adhoc.hakn <- NULL
  }
  #
  if (model != "prediction") {
    if (!with.subgroups && print[["stat"]]) {
      stat.lab <- "z"
      if (model == "random") {
        is.t <- interval %in% c("HK", "KR", "CR0", "CR1", "CR2")
        stat.lab <- if (all(is.t)) "t" else if (all(!is.t)) "z" else "z|t"
      }
      names(dat)[names(dat) == "statistic"] <- stat.lab
    }
    else
      dat$statistic <- NULL
    #
    if (!with.subgroups && print[["pval"]])
      names(dat)[names(dat) == "pval"] <- "p-value"
    else
      dat$pval <- NULL
  }
  #
  if (print[["se"]])
    names(dat)[names(dat) == "seTE"] <- paste0("SE(", lab.TE, ")")
  else
    dat$seTE <- NULL
  #
  if (with.subgroups && model != "prediction") {
    if (print[["pval.Q.b"]]) {
      dat$df <- NULL
      #
      dat %<>%
        mutate(
          Q.b = formatN(Q.b, digits = digits[["Q"]],
                        big.mark = big.mark, lab.NA = lab.NA),
          df.Q.b = vapply(df.Q.b, format_df, character(1),
                          lab.NA = lab.NA, big.mark = big.mark),
          pval.Q.b = formatPT(pval.Q.b, digits = digits[["pval.Q"]],
                              lab.NA = lab.NA, scientific = scientific.pval,
                              zero = zero.pval, JAMA = JAMA.pval)) %>%
        rename("d.f." = df.Q.b, "p-value" = pval.Q.b)
    }
    else
      dat[c("Q.b", "df.Q.b", "pval.Q.b", "test.Q.b")] <- NULL
  }
  #
  for (field in intersect(c("k", "df", "level"), names(dat))) {
    # Subgroup level headings are resolved in print.metabind().
    if (field == "level" && with.subgroups)
      next
    #
    if (isCol(dat, field))
      dat[[field]] <- formatN(dat[[field]],
                              digits = switch(field, k = 0, df = 2, level = 4),
                              lab.NA = lab.NA, big.mark = big.mark,
                              keep.allNA = TRUE)
  }
  #
  drop_allNAs(dat)
}

backtransf_metabind <- function(dat, pscale = 1, irscale = 1, suffix = "") {
  n <- paste0("n.harmonic.mean", suffix)
  time <- paste0("t.harmonic.mean", suffix)
  lower <- paste0("lower", suffix)
  upper <- paste0("upper", suffix)
  #
  for (i in seq_len(nrow(dat))) {
    if (is_prop(dat$sm[i]) || dat$sm[i] == "RD") {
      scale <- pscale
      n.i <- dat[[n]][i]
      time.i <- NULL
    }
    else if (is_rate(dat$sm[i]) || dat$sm[i] == "IRD") {
      scale <- irscale
      n.i <- NULL
      time.i <- dat[[time]][i]
    }
    else {
      scale <- 1
      n.i <- dat[[n]][i]
      time.i <- dat[[time]][i]
    }
    #
    for (field in paste0(c("TE", "lower", "upper"), suffix))
      dat[[field]][i] <- scale *
        backtransf(dat[[field]][i], dat$sm[i], n = n.i, time = time.i,
                   func = dat$func.backtransf[i],
                   args = dat$args.backtransf[[i]])
    #
    # Switch lower and upper limits
    #
    if (!is.na(dat[[lower]][i]) && !is.na(dat[[upper]][i]) &&
        dat[[lower]][i] > dat[[upper]][i]) {
      tmp.l <- dat[[lower]][i]
      dat[[lower]][i] <- dat[[upper]][i]
      dat[[upper]][i] <- tmp.l
    }
  }
  #
  dat
}

prepare_metabind_subgroups <- function(x) {
  # Get rid of warning 'no visible binding for global variable'
  common <- random <- prediction <- id <- name <- NULL
  #
  dat = as.data.frame(x, subgroup = FALSE, drop_all_NA = FALSE)
  #
  subdat <- as.data.frame(x, subgroup = TRUE, drop_all_NA = FALSE)
  #
  same.data <- x$samedata
  #
  same.results <- FALSE
  #
  # Compare selected effects separately from original heterogeneity metadata
  ignore <- c("id", "interval.slot", "name", "subgroup", "subgroup.name",
               "subgroup.levels",
               "common", "random", "prediction", "title", "complab", "outclab",
               "label.e", "label.c", "label.left", "label.right")
  #
  heterogeneity <- c("Q", "df.Q", "pval.Q", "tau2", "lower.tau2", "upper.tau2",
                     "se.tau2", "tau", "lower.tau", "upper.tau", "H", "lower.H",
                     "upper.H", "I2", "lower.I2", "upper.I2", "Rb", "lower.Rb",
                     "upper.Rb", "method.I2", "method.tau.ci", "level.hetstat",
                     "Q.Cochrane")
  tau.settings <- c("method.tau", "tau.preset", "TE.tau", "tau.common")
  ignore <- c(ignore, heterogeneity)
  #
  if (is.null(x$random) && is.null(x$predict))
    ignore <- c(ignore, tau.settings, "three.level", "rho")
  #
  if (is.null(x$predict))
    ignore <- c(ignore, "level.predict")
  #
  if (is.null(x$common) && is.null(x$random))
    ignore <- c(ignore, "level.ma")
  #
  comparison <- dat
  #
  rows <- match(dat$id, x$overall$id)
  #
  for (i in setdiff(names(x$overall), ignore))
    comparison[[i]] <- x$overall[[i]][rows]
  #
  comparison <- comparison[, setdiff(names(comparison), ignore), drop = FALSE]
  comparison$.interval.occurrence <- metabind_interval_occurrence(dat, comparison)
  #
  ids <- unique(dat$id)
  #
  shared.heterogeneity <- same.data && length(ids) > 0 &&
    nrow(unique(x$overall[match(ids, x$overall$id),
                          intersect(c(heterogeneity, tau.settings,
                                      "three.level", "rho"), names(x$overall)),
                          drop = FALSE])) == 1
  #
  if (same.data && length(ids) > 0) {
    reference <- comparison[dat$id == ids[1], , drop = FALSE]
    #
    same.results <- all(vapply(ids[-1], function(id)
      setequal(reference, comparison[dat$id == id, , drop = FALSE]),
      logical(1)))
  }
  #
  shared.overall <- same.data && same.results
  #
  if (shared.overall)
    keep <- !duplicated(comparison)
  else
    keep <- !duplicated(cbind(id = dat$id, comparison))
  #
  dat %<>% filter(keep) %>% select(-common, -random, -prediction)
  #
  if (shared.overall)
    dat %<>% select(-id, -name)
  #
  # Rename pooled results and replace parent results with subgroup values
  #
  fields <- c("TE", "seTE", "lower", "upper", "statistic", "pval", "df",
              "k", "n.harmonic.mean", "t.harmonic.mean")
  #
  for (i in intersect(fields, names(dat)))
    names(dat)[names(dat) == i] <- paste0(i, ".pooled")
  #
  for (i in fields) {
    field.w <- paste0(i, ".w")
    #
    if (field.w %in% names(subdat))
      subdat[[i]] <- as.vector(subdat[[field.w]])
    else
      subdat[[i]] <- rep(NA_real_, nrow(subdat))
    #
    subdat[[field.w]] <- NULL
  }
  #
  list(dat = dat, subdat = subdat, same.data = same.data,
       same.results = same.results, shared.overall = shared.overall,
       shared.heterogeneity = shared.heterogeneity)
}
