result_characteristics_fixture <- function(rate = FALSE, excluded = FALSE,
                                             multiple = FALSE) {
  if (rate) {
    dat <- data.frame(event = c(1, 4, 3, 8, 2, 6, 7, 10),
                      time = c(10, 30, 70, 150, 20, 50, 100, 250),
                      group = rep(c("A", "B"), each = 4))
    metarate(event, time, data = dat, sm = "IRFT", subgroup = group,
               exclude = if (excluded) seq_len(nrow(dat)) == 2 else FALSE,
               method.tau = if (excluded) "DL" else "REML",
               method.random.ci = if (multiple) c("classic", "HK") else "HK",
               prediction = TRUE, method.predict = if (multiple) "V" else c("V", "HK"))
  }
  else {
    dat <- data.frame(
      event.e = c(10, 20, 30, 40, 50, 30, 40, 50, 60, 70),
      event.c = c(20, 22, 25, 30, 35, 20, 25, 30, 35, 40),
      n.e = seq(100, 190, 10), n.c = seq(120, 210, 10),
      group = rep(c("A", "B"), each = 5))
    metabin(event.e, n.e, event.c, n.c, data = dat, method = "Inverse", sm = "OR",
               subgroup = group, tau.common = TRUE,
               exclude = if (excluded) seq_len(nrow(dat)) == 2 else FALSE,
               method.tau = if (excluded) "DL" else "REML",
               method.random.ci = if (multiple) c("classic", "HK") else "HK",
               prediction = TRUE, method.predict = if (multiple) "V" else c("V", "HK"))
  }
}

characteristics_extract <- function(x, model, overall = TRUE, subgroup = TRUE,
                                    n = TRUE, backtransf = FALSE,
                                    prediction.subgroup = TRUE) {
  args <- list(x = x, study.results = FALSE, common = model == "common",
               random = model == "random", prediction = model == "predict",
               overall = overall, subgroup = subgroup, n = n, se = TRUE,
               statistic = TRUE, pval = TRUE, backtransf = backtransf,
               digits = 12, digits.se = 12, digits.stat = 12, digits.pval = 12,
               writexl = FALSE)
  if (subgroup)
    args$prediction.subgroup <- if (model == "predict") prediction.subgroup else FALSE
  do.call(estimates, args)
}

expect_characteristics_source_alignment <- function(
    merged, sources, labels, model, overall = TRUE, subgroup = TRUE,
    n = TRUE, backtransf = FALSE, selected = rep(TRUE, length(sources)),
    prediction.subgroup = TRUE) {
  actual <- characteristics_extract(merged, model, overall, subgroup, n, backtransf,
                                      prediction.subgroup)
  expected <- list()
  add <- function(rows, label) {
    rows$studlab <- paste0(rows$studlab, " (", label, ")")
    if (!"subgroup" %in% names(rows))
      rows$subgroup <- ""
    if (!"df" %in% names(rows))
      rows$df <- Inf
    if (n) {
      for (field in c("n", "n.e", "n.c")) {
        if (!field %in% names(rows))
          rows[[field]] <- NA_real_
      }
    }
    rows
  }
  if (overall) {
    for (source in which(selected))
      expected[[length(expected) + 1L]] <- add(
        characteristics_extract(sources[[source]], model, TRUE, FALSE, n, backtransf),
        labels[source])
  }
  if (subgroup) {
    for (group in merged$subgroup.levels) {
      for (source in which(selected)) {
        rows <- characteristics_extract(sources[[source]], model, FALSE, TRUE,
                                          n, backtransf, prediction.subgroup)
        rows <- rows[rows$subgroup == group, , drop = FALSE]
        if (nrow(rows) > 0L)
          expected[[length(expected) + 1L]] <- add(rows, labels[source])
      }
    }
  }
  expected <- dplyr::bind_rows(expected)
  if (!subgroup)
    expected$subgroup <- NULL
  if (!any(is.finite(expected$df)))
    expected$df <- NULL
  if (n) {
    for (field in c("n", "n.e", "n.c")) {
      if (all(is.na(expected[[field]])))
        expected[[field]] <- NULL
    }
  }
  expect_equal(nrow(actual), nrow(expected))
  expect_setequal(names(actual), names(expected))
  for (field in names(actual))
    expect_equal(unname(actual[[field]]), unname(expected[[field]]),
                  tolerance = 1e-11, info = paste(model, field))
  invisible(actual)
}

test_that("native characteristics distinguish absent overall fields from subgroup fields", {
  ma <- metagen(seq(0.1, 0.6, length.out = 6), rep(0.2, 6),
                  subgroup = rep(c("A", "B"), each = 3),
                  method.random.ci = c("classic", "HK"), prediction = TRUE)
  ma$n.harmonic.mean.w <- c(A = 11, B = 22)
  ma$t.harmonic.mean.w <- c(A = 33, B = 44)
  ma$n.e.w <- c(A = 12, B = 34)
  ma$n.c.w <- c(A = 14, B = 36)
  expect_null(ma[["n.harmonic.mean"]])
  expect_null(ma[["t.harmonic.mean"]])
  chars <- meta:::result_characteristics(ma)
  for (model in c("common", "random", "predict")) {
    for (field in c("n", "n.e", "n.c", "n.harmonic.mean", "t.harmonic.mean"))
      expect_true(all(is.na(chars[[model]]$overall[[field]])))
    expect_equal(unname(chars[[model]]$subgroup$n.harmonic.mean[, 1]), c(11, 22))
    expect_equal(unname(chars[[model]]$subgroup$t.harmonic.mean[, 1]), c(33, 44))
  }
})

test_that("merged characteristics follow family-specific CI and PI slots with exclusions", {
  first <- result_characteristics_fixture(excluded = TRUE)
  second <- result_characteristics_fixture(multiple = TRUE)
  merged <- metamerge(first, second, label1 = "Excluded", label2 = "Full")
  chars <- merged$result.characteristics
  expect_named(chars, c("common", "random", "predict"))
  expect_equal(chars$common$overall$k, c(9, 10))
  expect_equal(chars$random$overall$k, c(9, 10, 10))
  expect_equal(chars$predict$overall$k, c(9, 9, 10))
  # Preserve the native extractor's overall-total convention.
  expect_equal(chars$random$overall$n.e, rep(sum(first$n.e), 3))
  expect_equal(chars$random$overall$n.c, rep(sum(first$n.c), 3))
  expect_equal(unname(chars$random$subgroup$k),
                unname(cbind(first$k.w, second$k.w, second$k.w)))
  expect_equal(unname(chars$predict$subgroup$k),
                unname(cbind(first$k.w, first$k.w, second$k.w)))
  for (model in names(chars)) {
    reference <- meta:::list2mat(merged[[paste0("lower.", model, ".w")]])
    expect_equal(nrow(chars[[model]]$overall), length(merged[[paste0("lower.", model)]]))
    for (field in names(chars[[model]]$subgroup))
      expect_identical(dimnames(chars[[model]]$subgroup[[field]]), dimnames(reference))
    expect_characteristics_source_alignment(merged, list(first, second),
                                             c("Excluded", "Full"), model)
  }
})

test_that("prediction-only and selected subgroup extraction work without displaying counts", {
  first <- result_characteristics_fixture(excluded = TRUE)
  second <- result_characteristics_fixture(multiple = TRUE)
  merged <- metamerge(first, second, label1 = "Excluded", label2 = "Full")
  for (show in c(TRUE, FALSE)) {
    expect_characteristics_source_alignment(
      merged, list(first, second), c("Excluded", "Full"), "predict",
      overall = TRUE, subgroup = FALSE, n = show)
    extracted <- expect_characteristics_source_alignment(
      merged, list(first, second), c("Excluded", "Full"), "predict",
      overall = FALSE, subgroup = TRUE, n = show,
      prediction.subgroup = c(TRUE, FALSE))
    expect_true(all(extracted$subgroup == "A"))
  }
  generic <- metagen(seq(0.1, 0.6, length.out = 6), rep(0.2, 6),
                       subgroup = rep(c("A", "B"), each = 3), prediction = TRUE,
                       method.predict = c("V", "HK"))
  for (show in c(TRUE, FALSE)) {
    extracted <- characteristics_extract(generic, "predict", overall = FALSE, n = show)
    expect_equal(nrow(extracted), 4)
    expect_false(any(c("n", "n.e", "n.c") %in% names(extracted)))
  }
})

test_that("nested characteristics preserve source-specific counts and model selections", {
  first <- result_characteristics_fixture(excluded = TRUE)
  second <- result_characteristics_fixture(multiple = TRUE)
  third <- result_characteristics_fixture()
  combined <- metamerge(first, second, label1 = "Excluded", label2 = "Full")
  nested <- metamerge(combined, third, label1 = "Combined", label2 = "Third")
  for (model in c("common", "random", "predict"))
    expect_characteristics_source_alignment(nested, list(combined, third),
                                             c("Combined", "Third"), model)
  selected <- metamerge(
    first, second, label1 = "Excluded", label2 = "Full",
    common1 = TRUE, common2 = FALSE, random1 = FALSE, random2 = TRUE,
    prediction1 = TRUE, prediction2 = FALSE)
  for (model in c("common", "random", "predict"))
    expect_characteristics_source_alignment(
      selected, list(first, second), c("Excluded", "Full"), model,
      selected = if (model == "random") c(FALSE, TRUE) else c(TRUE, FALSE))

  only.third <- metamerge(
    combined, third, common1 = FALSE, random1 = FALSE, prediction1 = FALSE,
    common2 = TRUE, random2 = TRUE, prediction2 = TRUE)
  for (model in c("common", "random", "predict")) {
    expect_equal(only.third$result.characteristics[[model]]$overall,
                  meta:::result_characteristics(third)[[model]]$overall)
  }
})

test_that("unknown counts in one source are not borrowed from another source", {
  dat <- data.frame(TE = seq(0.1, 0.6, length.out = 6), seTE = rep(0.2, 6),
                    group = rep(c("A", "B"), each = 3),
                    n.e = seq(20, 70, 10), n.c = seq(25, 75, 10))
  unknown <- metagen(TE, seTE, data = dat, subgroup = group, keepdata = FALSE,
                       method.random.ci = "HK", prediction = TRUE)
  known <- metagen(TE, seTE, data = dat, subgroup = group, n.e = n.e, n.c = n.c,
                     method.random.ci = c("classic", "HK"), prediction = TRUE)
  merged <- metamerge(unknown, known, label1 = "Unknown", label2 = "Known")
  for (field in c("n", "n.e", "n.c")) {
    expect_true(is.na(merged$result.characteristics$random$overall[[field]][1]))
    expect_true(all(is.na(merged$result.characteristics$random$subgroup[[field]][, 1])))
  }
  for (model in c("common", "random", "predict"))
    expect_characteristics_source_alignment(merged, list(unknown, known),
                                             c("Unknown", "Known"), model)
})

test_that("rate back-transformations use source-specific overall and subgroup harmonic means", {
  first <- result_characteristics_fixture(rate = TRUE, excluded = TRUE)
  second <- result_characteristics_fixture(rate = TRUE, multiple = TRUE)
  expect_gt(abs(first$t.harmonic.mean - second$t.harmonic.mean), 0.1)
  merged <- metamerge(first, second, label1 = "Excluded", label2 = "Full")
  expect_equal(merged$result.characteristics$random$overall$t.harmonic.mean,
                c(first$t.harmonic.mean, rep(second$t.harmonic.mean, 2)))
  expect_equal(unname(merged$result.characteristics$predict$subgroup$t.harmonic.mean),
                unname(cbind(first$t.harmonic.mean.w, first$t.harmonic.mean.w,
                             second$t.harmonic.mean.w)))
  for (show in c(TRUE, FALSE)) {
    for (model in c("common", "random", "predict"))
      expect_characteristics_source_alignment(
        merged, list(first, second), c("Excluded", "Full"), model,
        n = show, backtransf = TRUE)
  }
})

test_that("study characteristics follow study sorting and are excluded from data-frame coercion", {
  ma <- metaprop(c(3, 5, 8), c(20, 50, 100), method = "Inverse", sm = "PFT",
                  subgroup = c("A", "A", "B"))
  for (show in c(TRUE, FALSE)) {
    extracted <- estimates(ma, study.results = TRUE, common = FALSE, random = FALSE,
                            prediction = FALSE, overall = FALSE, subgroup = FALSE,
                            sortvar = rev(seq_along(ma$TE)), n = show,
                            backtransf = TRUE, digits = 12, writexl = FALSE)
    expect_identical(extracted$studlab, rev(ma$studlab))
    expect_equal(extracted$estimate, round(rev(ma$event / ma$n), 12), tolerance = 1e-11)
    if (show)
      expect_equal(extracted$n, rev(ma$n))
    else
      expect_false("n" %in% names(extracted))
  }
  plain <- metagen(c(0.1, 0.3, 0.5), rep(0.1, 3))
  merged <- metamerge(plain, plain)
  before <- merged$result.characteristics
  expect_length(before, 3)
  dat <- as.data.frame(merged)
  expect_equal(nrow(dat), 3)
  expect_false(any(startsWith(names(dat), "result.characteristics")))
  expect_identical(merged$result.characteristics, before)
})
