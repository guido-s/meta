map_metabind_subgroups <- function(ma, x, dat, display, shared.overall,
                                   common, random, prediction, overall) {
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
  families <- c(common = "Common effect model",
                random = "Random effects model",
                predict = "Prediction interval")
  selected <- c(common = common, random = random, predict = prediction)
  fields <- c("TE", "seTE", "lower", "upper", "statistic", "pval", "df")
  #
  for (family in names(families)) {
    rows <- which(dat$type == families[[family]] & selected[[family]])
    keys <- intersect(c("method.interval", "adhoc.hakn",
                        if (family == "predict") "level.predict"
                        else "level.ma",
                        "text"), names(dat))
    slots <- dat[rows, keys, drop = FALSE]
    #
    if (family == "common")
      slots <- slots[seq_len(min(nrow(slots), 1L)), , drop = FALSE]
    else
      slots <- unique(slots)
    #
    ns <- max(nrow(slots), 1L)
    #
    slot <- integer(length(rows))
    #
    if (family == "common")
      slot[] <- 1L
    else {
      for (j in seq_len(nrow(slots))) {
        same <- rep(TRUE, length(rows))
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
    global.rows[[family]] <- global
    group.rows[[family]] <- grouped
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
  args$test.effect.subgroup.common <- args$common.subgroup
  args$test.effect.subgroup.random <- unname(rowSums(args$random.subgroup) > 0)
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
    if (shared.overall)
      ma[[field]] <- values[metadata[1]]
    else if (field %in% heterogeneity)
      ma[[field]] <- NA_real_
    else
      ma[[field]] <- values[metadata[1]]
    #
    if (shared.overall)
      ma[[paste0(field, ".w")]] <-
        rep(if (is.numeric(values)) NA_real_ else NA_character_, n.by)
    else
      ma[[paste0(field, ".w")]] <- values[metadata]
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
  ma$k.w.orig <- x$overall$k[metadata]
  #
  if (shared.overall)
    ma$null.effect <- replaceNULL(x$overall$null.effect[metadata[1]])
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

prepare_metabind_subgroups <- function(
    x, dat = as.data.frame(x, subgroup = FALSE, drop_all_NA = FALSE)) {
  # Get rid of warning 'no visible binding for global variable'
  common <- random <- prediction <- id <- name <- NULL
  #
  # Repeated parent settings must not multiply the source subgroup rows
  pooled <- x
  #
  for (i in c("common", "random", "predict")) {
    if (!is.null(pooled[[i]]))
      pooled[[i]] <- unique(pooled[[i]])
  }
  subdat <- as.data.frame(pooled, subgroup = TRUE, drop_all_NA = FALSE)
  same.data <- inherits(x$data, "data.frame") &&
    !".not_a_dataset" %in% names(x$data)
  #
  same.results <- FALSE
  #
  # Ignore labels, but compare original settings and unrounded results
  ignore <- c("id", "name", "subgroup", "subgroup.name", "subgroup.levels",
              "common", "random", "prediction", "title", "complab", "outclab",
              "label.e", "label.c", "label.left", "label.right")
  #
  comparison <- dat
  #
  rows <- match(dat$id, x$overall$id)
  #
  for (i in setdiff(names(x$overall), ignore))
    comparison[[i]] <- x$overall[[i]][rows]
  #
  comparison <- comparison[, setdiff(names(comparison), ignore), drop = FALSE]
  #
  ids <- unique(dat$id)
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
  list(dat = dat, subdat = subdat, same.data = same.data,
       same.results = same.results, shared.overall = shared.overall)
}
