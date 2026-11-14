test_that("shared summaries render stored estimates rather than pooling source rows", {
  x <- metabind_render_fixture()
  expect_true(all(x$overall$tau2 > 0))
  expect_true(all(abs(x$random$TE) > 0.01))
  drawn <- metabind_render_pdf(x, backtransf = FALSE)
  ma <- drawn$args$x
  expect_true(drawn$args$overall)
  expect_true(drawn$args$random)
  expect_equal(ma$TE.random, x$random$TE[1])
  expect_equal(ma$lower.random, x$random$lower[1])
  expect_equal(ma$upper.random, x$random$upper[1])
  expect_equal(ma$TE.common, x$common$TE[1])
  expect_equal(ma$lower.predict, x$predict$lower[1])
  expect_equal(ma$upper.predict, x$predict$upper[1])
  expect_true(all(is.na(ma$TE.random.w)))
  expect_false(any(drawn$args$random.subgroup))

  source <- drawn$result$data
  rows <- source$type == "Random effects model"
  repooled <- metagen(source$TE.orig[rows], source$seTE.orig[rows],
                       method.tau = "DL")
  expect_gt(abs(ma$TE.random - repooled$TE.random), 1e-5)
  expect_type(drawn$result, "list")
})

test_that("shared common effects do not present differing heterogeneity as shared", {
  x <- metabind_render_fixture(method.tau = "REML")
  for (hetstat in c(TRUE, FALSE)) {
    drawn <- metabind_render_pdf(
      x, random = FALSE, prediction = FALSE, backtransf = FALSE,
      hetstat = hetstat, overall.hetstat = TRUE, lab.NA = "MISSING",
      leftcols = c("studlab", "effect", "ci", "I2", "tau2"))
    expect_true(drawn$args$overall)
    expect_equal(drawn$args$x$TE.common, x$common$TE[1])
    expect_false(drawn$args$overall.hetstat)
    expect_false(drawn$args$hetstat)
    expect_false(drawn$args$subgroup.hetstat)
    expect_true(is.na(drawn$args$x$Q))
    expect_equal(unname(drawn$args$x$tau2.w), unlist(x$overall$tau2))
    expect_true(all(c("I2", "tau2") %in% names(drawn$args$data.pooled)))
    for (field in c("I2", "tau2"))
      expect_true(all(is.na(drawn$args$data.pooled[[field]]) |
                        drawn$args$data.pooled[[field]] == "MISSING"))
    labels <- vapply(drawn$labels, function(label)
      paste(deparse(label), collapse = " "), character(1))
    expect_false(any(grepl("Heterogeneity", labels, fixed = TRUE)))
  }
})

test_that("different tau methods or source data render component stored summaries", {
  for (reason in c("method.tau", "source data")) {
    x <- metabind_render_fixture(
      method.tau = if (reason == "method.tau") "REML" else "DL")
    if (reason == "source data")
      x$samedata <- FALSE
    drawn <- metabind_render_pdf(x, backtransf = FALSE)
    ma <- drawn$args$x
    expect_false(drawn$args$overall, info = reason)
    expect_true(all(is.na(ma$TE.random)), info = reason)
    expect_true(all(drawn$args$random.subgroup), info = reason)
    for (id in x$overall$id) {
      g <- match(as.character(id), as.character(ma$subgroup.levels))
      r <- match(id, x$random$id)
      expect_equal(unname(ma$TE.random.w[g, ]), x$random$TE[r], info = reason)
      expect_equal(unname(ma$lower.random.w[g, ]), x$random$lower[r])
      expect_equal(unname(ma$upper.random.w[g, ]), x$random$upper[r])
      expect_equal(unname(ma$TE.common.w[g, ]),
                   x$common$TE[match(id, x$common$id)])
    }
  }
})

test_that("explicit overall FALSE and rejection of public subgroup survive rendering", {
  x <- metabind_render_fixture()
  drawn <- metabind_render_pdf(x, overall = FALSE)
  expect_false(drawn$args$overall)
  expect_false(drawn$args$x$overall)
  expect_equal(drawn$args$x$TE.random, x$random$TE[1])
  expect_error(metabind_render_pdf(x, subgroup = type),
               "Argument 'subgroup' is unavailable", fixed = TRUE)
  expect_error(metabind_render_pdf(x, subgroup = NULL),
                "Argument 'subgroup' is unavailable", fixed = TRUE)
})

test_that("forest sharing depends on displayed models and their interval levels", {
  x <- metabind_render_fixture()
  first <- x$meta[[1]]
  second <- x$meta[[2]]
  for (difference in c("random CI", "prediction level")) {
    other <- if (difference == "random CI")
      update(second, method.random.ci = "HK")
    else
      update(second, level.predict = 0.90)
    mb <- metabind(first, other)
    drawn <- metabind_render_pdf(mb, random = FALSE, prediction = FALSE,
                                  backtransf = FALSE)
    expect_true(drawn$args$overall, info = difference)
    expect_equal(drawn$args$x$TE.common, first$TE.common, info = difference)

    hidden <- metabind_render_pdf(mb, random = FALSE, prediction = FALSE,
                                   overall = FALSE, backtransf = FALSE)
    expect_false(hidden$args$overall, info = difference)
    expect_false(hidden$args$x$overall, info = difference)

    selected <- metabind_render_pdf(
      mb, random = difference == "random CI",
      prediction = difference == "prediction level", backtransf = FALSE)
    expect_false(selected$args$overall, info = difference)
  }

  mb <- metabind(first, update(second, level.ma = 0.90))
  drawn <- metabind_render_pdf(mb, common = FALSE, random = FALSE,
                                prediction = TRUE, backtransf = FALSE)
  expect_true(drawn$args$overall)
  expect_equal(drawn$args$x$lower.predict, first$lower.predict)
  expect_equal(drawn$args$x$upper.predict, first$upper.predict)
})

test_that("source numeric effects, intervals and study counts use subgroup fields", {
  x <- metabind_render_fixture(method.tau = "REML")
  raw <- as.data.frame(x, subgroup = TRUE, drop_all_NA = FALSE)
  drawn <- metabind_render_pdf(x, backtransf = FALSE, print = list(se = TRUE))
  dat <- drawn$result$data
  expect_equal(dat$id, raw$id)
  expect_equal(as.character(dat$subgroup.levels),
               as.character(raw$subgroup.levels))
  for (field in c("TE", "seTE", "lower", "upper")) {
    expect_equal(dat[[paste0(field, ".orig")]], raw[[paste0(field, ".w")]],
                 info = field)
    expect_gt(max(abs(raw[[paste0(field, ".w")]] - raw[[field]]),
                  na.rm = TRUE), 1e-5)
  }
  expect_equal(dat$k, as.character(raw$k.w))
  expect_true(all(raw$k.w == 4))
  expect_true(all(raw$k == 8))
  prediction <- dat$type == "Prediction interval"
  expect_true(any(prediction))
  expect_true(all(is.na(drawn$args$x$TE[prediction])))
  expect_true(all(is.na(drawn$args$x$seTE[prediction])))
  expect_false(any(drawn$result$seTE.format == "1.0000"))
  expect_equal(drawn$args$x$lower[prediction], raw$lower.w[prediction])
  expect_equal(drawn$args$x$upper[prediction], raw$upper.w[prediction])
})

test_that("multiple random and PI methods retain independent stored slots", {
  x <- metabind_render_fixture(method.random.ci = c("classic", "HK"))
  x <- metabind(update(x$meta[[1]], method.predict = c("V", "HK")),
                update(x$meta[[2]], method.predict = c("V", "HK")),
                keepdata = TRUE)
  for (grouped in c(FALSE, TRUE)) {
    y <- x
    if (grouped)
      y$samedata <- FALSE
    drawn <- metabind_render_pdf(y, backtransf = FALSE)
    ma <- drawn$args$x
    expect_equal(ma$method.random.ci, c("classic", "HK"))
    expect_equal(ma$method.predict, c("V", "HK"))
    for (family in c("random", "predict")) {
      stored <- y[[family]]
      for (id in if (grouped) y$overall$id else y$overall$id[1]) {
        g <- match(as.character(id), as.character(ma$subgroup.levels))
        rows <- stored$id == id
        for (field in c(if (family == "random") "TE", "lower", "upper")) {
          actual <- if (grouped)
            unname(ma[[paste0(field, ".", family, ".w")]][g, ])
          else ma[[paste0(field, ".", family)]]
          expect_equal(actual, stored[[field]][rows],
                       info = paste(grouped, family, id, field))
        }
      }
    }
  }
})

test_that("sparse component families and prediction-only selections render", {
  x <- metabind_render_fixture(method.random.ci = c("classic", "HK"))
  x$samedata <- FALSE
  x$common <- x$common[x$common$id == 1, , drop = FALSE]
  x$subgroup$common <- subset(x$subgroup$common, id == 1)
  x$random <- subset(x$random, id == 2 & method.interval == "HK")
  x$subgroup$random <- subset(x$subgroup$random,
                              id == 2 & method.interval == "HK")
  x$predict <- subset(x$predict, id == 1)
  x$subgroup$predict <- subset(x$subgroup$predict, id == 1)
  drawn <- metabind_render_pdf(x, backtransf = FALSE)
  expect_equal(drawn$args$common.subgroup, c(TRUE, FALSE))
  expect_equal(unname(drawn$args$random.subgroup[, 1]), c(FALSE, TRUE))
  expect_equal(unname(drawn$args$prediction.subgroup[, 1]), c(TRUE, FALSE))
  expect_true(is.na(drawn$args$x$TE.common.w[2, 1]))
  expect_true(is.na(drawn$args$x$TE.random.w[1, 1]))
  expect_true(is.na(drawn$args$x$lower.predict.w[2, 1]))

  for (grouped in c(FALSE, TRUE)) {
    y <- metabind_render_fixture()
    if (grouped)
      y$samedata <- FALSE
    pi <- metabind_render_pdf(y, common = FALSE, random = FALSE,
                              prediction = TRUE, backtransf = FALSE)
    expect_false(pi$args$common)
    expect_false(pi$args$random)
    expect_true(all(pi$result$data$type == "Prediction interval"))
    expect_true(all(is.na(pi$args$x$TE)))
    expect_equal(pi$args$x$lower, pi$result$data$lower.orig)
    expect_equal(pi$args$x$upper, pi$result$data$upper.orig)
    if (grouped) {
      expect_true(all(pi$args$prediction.subgroup))
      expect_equal(unname(pi$args$x$lower.predict.w[, 1]), y$predict$lower)
    }
    else {
      expect_true(pi$args$prediction)
      expect_true(pi$args$overall)
      expect_equal(pi$args$x$lower.predict, y$predict$lower[1])
    }
  }
})

test_that("sparse method unions do not couple random CI and PI slots", {
  x <- metabind_render_fixture()
  x <- metabind(
    update(x$meta[[1]], method.random.ci = "classic",
           method.predict = c("V", "HK"),
           text.predict = c("Prediction interval (V)", "Prediction interval (HK)")),
    update(x$meta[[2]], method.random.ci = "HK", method.predict = "V",
           text.predict = "Prediction interval (V)"),
    keepdata = TRUE)
  drawn <- metabind_render_pdf(x, backtransf = FALSE)
  ma <- drawn$args$x
  expect_equal(ma$method.random.ci, c("classic", "HK"))
  expect_equal(ma$method.predict, c("V", "HK"))
  expect_equal(unname(drawn$args$random.subgroup),
               matrix(c(TRUE, FALSE, FALSE, TRUE), 2))
  expect_equal(unname(drawn$args$prediction.subgroup),
               matrix(c(TRUE, TRUE, TRUE, FALSE), 2))
  expect_equal(unname(ma$TE.random.w[cbind(1:2, 1:2)]), x$random$TE)
  expect_equal(unname(ma$lower.predict.w[cbind(c(1, 1, 2), c(1, 2, 1))]),
               x$predict$lower)
  expect_equal(unname(ma$upper.predict.w[cbind(c(1, 1, 2), c(1, 2, 1))]),
               x$predict$upper)
  expect_true(is.na(ma$lower.predict.w[2, 2]))

  pi <- metabind_render_pdf(x, common = FALSE, random = FALSE,
                            prediction = TRUE, backtransf = FALSE)
  expect_equal(pi$args$x$method.predict, c("V", "HK"))
  expect_equal(pi$args$x$lower.predict.w, ma$lower.predict.w)
  expect_equal(pi$args$x$upper.predict.w, ma$upper.predict.w)
  expect_true(all(is.na(pi$args$x$TE)))
})

test_that("pooled display columns follow global then family and component row order", {
  for (grouped in c(FALSE, TRUE)) {
    x <- metabind_render_fixture(method.random.ci = c("classic", "HK"))
    if (grouped)
      x$samedata <- FALSE
    drawn <- metabind_render_pdf(
      x, backtransf = FALSE, lab.NA = ".", digits = list(tau2 = 4),
      leftcols = c("studlab", "k", "method.interval", "tau2", "effect", "ci"),
      rightcols = FALSE, print = list(tau2 = TRUE))
    ma <- drawn$args$x
    pooled <- drawn$args$data.pooled
    expected <- c(ma$TE.common, ma$TE.random, ma$TE.predict,
                  as.vector(t(ma$TE.common.w)),
                  as.vector(t(ma$TE.random.w)),
                  as.vector(t(ma$TE.predict.w)))
    expect_equal(pooled$TE.orig, unname(expected))
    expected.lower <- c(ma$lower.common, ma$lower.random, ma$lower.predict,
                        as.vector(t(ma$lower.common.w)),
                        as.vector(t(ma$lower.random.w)),
                        as.vector(t(ma$lower.predict.w)))
    expect_equal(pooled$lower.orig, unname(expected.lower))
    expect_equal(nrow(pooled), length(expected))
    present <- !is.na(pooled$lower.orig)
    expect_equal(pooled$k[present], rep("8", sum(present)))
    expect_equal(pooled$method.interval[present],
                 if (grouped) c("classic", "classic", "classic", "HK",
                                "classic", "HK", "V", "V")
                 else c("classic", "classic", "HK", "V"))
    random <- which(pooled$type == "Random effects model")
    first.random <- random[if (grouped) c(1, 3) else 1]
    expect_equal(trimws(pooled$tau2[first.random]),
                 sprintf("%.4f", unlist(x$overall$tau2[if (grouped) 1:2 else 1])))
    expect_equal(drawn$result$data.pooled, pooled)
    expect_equal(drawn$result$leftcols,
                 c("studlab", "k", "method.interval", "tau2", "effect", "ci"))
  }
})

test_that("subgroup rendering preserves options and non-exponential transformations", {
  x <- metabind_render_fixture()
  x$overall$func.backtransf <- "plogis"
  x$overall$args.backtransf <- rep(list(list(location = 0.2, scale = 2)), 2)
  drawn <- metabind_render_pdf(
    x, col = list(square = "navy"), type = list(random = "circle"),
    digits = 4, print = list(backtransf = TRUE, args.backtransf = TRUE),
    leftcols = c("studlab", "func.backtransf", "args.backtransf", "effect", "ci"),
    rightcols = FALSE)
  dat <- drawn$result$data
  expect_equal(dat$TE, plogis(dat$TE.orig, location = 0.2, scale = 2))
  expect_equal(dat$lower, plogis(dat$lower.orig, location = 0.2, scale = 2))
  expect_equal(dat$upper, plogis(dat$upper.orig, location = 0.2, scale = 2))
  expect_false(drawn$args$x$backtransf)
  expect_true(all(dat$.type.study[dat$type == "Random effects model"] == "circle"))
  expect_true(all(drawn$args$col$square[dat$.type.study == "square"] == "navy"))
  expect_true(all(c("func.backtransf", "args.backtransf") %in%
                    drawn$result$leftcols))
  expect_equal(drawn$args$x$TE.random,
               plogis(x$random$TE[1], location = 0.2, scale = 2))
})

test_that("shared heterogeneity uses native rows and selected component metadata", {
  x <- metabind_render_fixture()
  drawn <- metabind_render_pdf(x)
  expect_true(drawn$args$hetstat)
  expect_true(drawn$args$overall.hetstat)
  expect_equal(drawn$args$x$Q, x$overall$Q[1])
  expect_equal(drawn$args$x$df.Q, x$overall$df.Q[1])
  expect_false(any(c("tau2", "tau", "I2", "Q", "pval.Q") %in%
                     drawn$result$leftcols))
  expect_false(metabind_render_pdf(x, hetstat = FALSE)$args$overall.hetstat)

  x$common <- subset(x$common, id == 2)
  x$subgroup$common <- subset(x$subgroup$common, id == 2)
  x$overall$Q[1] <- x$overall$Q[2] + 10
  drawn <- metabind_render_pdf(x, random = FALSE, prediction = FALSE)
  expect_equal(drawn$args$x$Q, x$overall$Q[2])
  expect_equal(as.character(drawn$args$x$subgroup.levels), "2")
})

test_that("generic subgroup inputs do not partially match harmonic mean fields", {
  ma <- metagen(c(-0.4, 0.8, 0.5, -0.3, 0.9, 0.2),
                 c(0.1, 0.3, 0.2, 0.2, 0.4, 0.15),
                 subgroup = rep(c("A", "B"), each = 3), method.tau = "DL")
  x <- metabind(ma, update(ma, subgroup = rep(c("C", "D"), 3)))
  drawn <- metabind_render_pdf(x)
  expect_true(all(is.na(x$overall$n.harmonic.mean)))
  expect_equal(drawn$args$x$TE.random, x$random$TE[1])
})

test_that("PFT subgroup rows use their own harmonic means", {
  ma <- metaprop(c(2, 10, 8, 30, 40, 6), c(10, 50, 35, 120, 150, 60),
                  sm = "PFT", subgroup = rep(c("A", "B"), each = 3),
                  method.tau = "DL")
  x <- metabind(ma, update(ma, subgroup = rep(c("C", "D"), 3)))
  raw <- as.data.frame(x, drop_all_NA = FALSE)
  expect_true("n.harmonic.mean.w" %in% names(raw))
  drawn <- metabind_render_pdf(x)
  dat <- drawn$result$data
  expect_equal(dat$n.harmonic.mean, raw$n.harmonic.mean.w)
  expected <- vapply(seq_len(nrow(dat)), function(i)
    backtransf(dat$TE.orig[i], "PFT", n = raw$n.harmonic.mean.w[i]), numeric(1))
  expect_equal(dat$TE, expected)
  expect_equal(drawn$args$x$TE.random,
               backtransf(x$random$TE[1], "PFT", n = ma$n.harmonic.mean))
})

test_that("source subgroup tests are displayed for common and random results", {
  cols <- c("Q.b", "df.Q.b", "pval.Q.b")
  for (grouped in c(FALSE, TRUE)) {
    x <- metabind_render_fixture(method.random.ci = c("classic", "HK"))
    if (grouped)
      x$samedata <- FALSE
    raw <- as.data.frame(x, drop_all_NA = FALSE)
    drawn <- metabind_render_pdf(x, digits = list(Q = 3, pval.Q = 4))
    expect_true(all(cols %in% drawn$result$leftcols))
    dat <- drawn$result$data
    model <- dat$type != "Prediction interval" &
      as.character(dat$subgroup.levels) %in% c("A", "C")
    expect_true(all(!is.na(raw$df.Q.b[model])))
    expect_equal(dat$Q.b[model],
                 meta:::formatN(raw$Q.b[model], digits = 3, lab.NA = ""))
    expect_equal(dat$df.Q.b[model],
                 vapply(raw$df.Q.b[model], meta:::format_df,
                        character(1), lab.NA = ""))
    expect_equal(dat$pval.Q.b[model],
                 meta:::formatPT(raw$pval.Q.b[model], digits = 4,
                                 lab.NA = "", lab = FALSE, trim = "left",
                                 scientific = FALSE, zero = TRUE, JAMA = FALSE))
    for (i in cols) {
      expect_true(all(dat[[i]][!model] == ""))
      expect_true(all(is.na(drawn$args$data.pooled[[i]]) |
                        drawn$args$data.pooled[[i]] == ""))
    }
    expect_false(drawn$args$test$subgroup)
  }
})

test_that("stored overall-effect tests appear for shared and component summaries", {
  for (grouped in c(FALSE, TRUE)) {
    x <- metabind_render_fixture(method.random.ci = c("classic", "HK"))
    if (grouped)
      x$samedata <- FALSE
    drawn <- metabind_render_pdf(x, digits = list(stat = 3, pval = 4),
                                 scientific.pval = FALSE)
    labels <- Filter(function(z) is.call(z) &&
                       grepl("Test for overall effect", as.character(z[[2]]), fixed = TRUE),
                     drawn$labels)
    expect_length(labels, if (grouped) 4 else 2)
    expect_identical(drawn$args$test$overall, !grouped)
    expect_identical(drawn$args$test$effect.subgroup, grouped)
    expect_false(drawn$args$test$subgroup)
    expect_false(any(c("statistic", "pval", "df") %in% drawn$result$leftcols))
    for (i in seq_along(labels)) {
      family <- if (grepl("random", labels[[i]][[2]])) "random" else "common"
      id <- if (grouped) ceiling(i / 2) else 1
      row <- which(x[[family]]$id == id)[1]
      expect_equal(unname(labels[[i]][[5]]),
                   unname(meta:::formatN(x[[family]]$statistic[row], 3, "NA")))
      expect_equal(trimws(unname(labels[[i]][[8]])),
                   trimws(unname(meta:::formatPT(x[[family]]$pval[row], digits = 4,
                     lab = TRUE, labval = "", zero = TRUE, JAMA = FALSE,
                     scientific = FALSE, lab.NA = "NA"))))
    }
  }
})

test_that("nonshared component groups display their stored heterogeneity", {
  for (different.data in c(FALSE, TRUE)) {
    x <- metabind_render_fixture(method.tau = "REML")
    if (different.data)
      x$samedata <- FALSE
    drawn <- metabind_render_pdf(x)
    expect_true(drawn$args$hetstat)
    expect_true(drawn$args$subgroup.hetstat)
    expect_false(drawn$args$overall.hetstat)
    ma <- drawn$args$x
    expect_equal(unname(ma$Q.w), x$overall$Q)
    expect_equal(unname(ma$I2.w), x$overall$I2)
    expect_equal(unname(ma$pval.Q.w), x$overall$pval.Q)
    expect_equal(unname(ma$k.w.orig), x$overall$k)
    expect_false(identical(ma$k.w.orig, ma$k.w))
    labels <- vapply(drawn$labels, function(z) paste(deparse(z), collapse = " "), "")
    expect_equal(sum(grepl("Heterogeneity", labels, fixed = TRUE)), 2)
    disabled <- metabind_render_pdf(x, hetstat = FALSE)
    expect_false(disabled$args$subgroup.hetstat)
  }
})

test_that("component sorting keeps headings, pooled counts and estimates aligned", {
  components <- lapply(seq_len(11), function(i) {
    n <- i + 4
    metagen(i + seq_len(n) / 10, rep(0.2, n),
             subgroup = rep(c("A", "B"), length.out = n),
             method.tau = "DL")
  })
  x <- metabind(components, name = paste0("Component", seq_len(11)))
  for (sorted in c(FALSE, TRUE)) {
    drawn <- metabind_render_pdf(x, sort.subgroup = sorted,
                                 backtransf = FALSE, details = FALSE,
                                 test.overall = FALSE, hetstat = FALSE,
                                 leftcols = c("studlab", "k", "name", "effect", "ci"))
    ids <- as.character(seq_len(11))
    if (sorted)
      ids <- sort(ids)
    rows <- match(ids, as.character(x$overall$id))
    expected.names <- x$overall$name[rows]
    expect_false(drawn$args$sort.subgroup)
    expect_identical(as.character(drawn$args$x$subgroup.levels), ids)
    expect_equal(drawn$args$text$subgroup, expected.names)
    expect_equal(unname(drawn$args$x$k.w.orig), x$overall$k[rows])
    expect_equal(unname(drawn$args$x$k.w), rep(4L, 11))

    labels <- lapply(drawn$columns$col.studlab$labels, `[[`, "label")
    headings <- unlist(Filter(function(z)
      is.character(z) && length(z) == 1 && grepl("^Component", z), labels))
    expect_equal(unname(headings), expected.names)

    labels <- lapply(drawn$columns$col.k$labels, `[[`, "label")
    counts <- unlist(Filter(function(z)
      is.character(z) && length(z) == 1 && grepl("^[0-9]+$", z), labels))
    expect_equal(unname(counts[seq_len(22)]), rep(as.character(x$overall$k[rows]), 2))
    expect_equal(drawn$args$data.pooled$name[!is.na(drawn$args$data.pooled$name)],
                 rep(expected.names, 2))
    effect <- drawn$result$effect.format
    effect <- effect[!is.na(effect) & nzchar(effect)]
    expected <- c(x$common$TE[match(ids, as.character(x$common$id))],
                  x$random$TE[match(ids, as.character(x$random$id))])
    expect_equal(unname(effect[seq_len(22)]),
                 unname(meta:::formatN(expected, digits = 2, lab.NA = "")))
  }
})
