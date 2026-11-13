subgroup_inference_fixture <- function(methods = "HK") {
  dat <- data.frame(event = c(10, 20, 30, 40, 50, 25, 35, 45, 55, 65),
                    total = 100, group = rep(c("A", "B"), each = 5))
  metaprop(event, total, data = dat, method = "GLMM",
           method.random.ci = methods, subgroup = group,
           common = TRUE, tau.common = TRUE)
}

subgroup_inference_print <- function(x, ...) {
  capture.output(print(x, details.methods = FALSE, digits.pval.Q = 4,
                       scientific.pval = FALSE, zero.pval = TRUE,
                        JAMA.pval = FALSE, ...))
}

subgroup_inference_forest <- function(x, ...) {
  grDevices::pdf(NULL, width = 16, height = 12)
  on.exit(grDevices::dev.off())
  grob <- grid::grid.grabExpr(
    forest(x, layout = "RevMan5", print.Q.subgroup = TRUE,
           label.test.subgroup.common = "INFERENCE_COMMON",
           label.test.subgroup.random = "INFERENCE_RANDOM",
           digits.pval.Q = 4, scientific.pval = FALSE,
            zero.pval = TRUE, JAMA.pval = FALSE, ...),
    width = 16, height = 12)
  labels <- function(g) {
    c(if (inherits(g, "text"))
        paste(deparse(g$label, width.cutoff = 500L), collapse = " "),
      unlist(lapply(g$children, labels), use.names = FALSE))
  }
  result <- labels(grob)
  result[grepl("INFERENCE_COMMON|INFERENCE_RANDOM", result)]
}

test_that("GLMM subgroup coefficients use joint-fit DF in either method order", {
  singles <- lapply(c("classic", "HK"), subgroup_inference_fixture)
  names(singles) <- c("classic", "HK")
  fits <- c(singles, list(subgroup_inference_fixture(c("classic", "HK")),
                         subgroup_inference_fixture(c("HK", "classic"))))
  fields <- c("seTE.random.w", "df.random.w", "statistic.random.w",
              "pval.random.w", "lower.random.w", "upper.random.w")
  for (x in fits) {
    for (method in x$method.random.ci) {
      column <- match(method, x$method.random.ci)
      values <- lapply(x[fields], function(value) {
        if (is.matrix(value)) unname(value[, column]) else unname(value)
      })
      reference <- lapply(singles[[method]][fields], unname)
      expect_equal(values, reference)
      expect_equal(x$TE.random.w, singles[[method]]$TE.random.w)

      # Both coefficients come from one ten-study, two-coefficient fit.
      df <- if (method == "HK") length(x$TE) - 2 else Inf
      expect_equal(values$df.random.w, rep(df, 2))
      statistic <- unname(x$TE.random.w) / values$seTE.random.w
      expect_equal(values$statistic.random.w, statistic)
      expect_equal(values$pval.random.w,
                   2 * pt(abs(statistic), df = values$df.random.w,
                          lower.tail = FALSE))
      margin <- qt((1 + x$level.ma) / 2, df = values$df.random.w) *
        values$seTE.random.w
      expect_equal(values$lower.random.w, unname(x$TE.random.w) - margin)
      expect_equal(values$upper.random.w, unname(x$TE.random.w) + margin)
    }
  }
})

test_that("subgroup P fallback distinguishes Q from F and retains method names", {
  statistic <- c(classic = 4.2, HK = 3.1)
  dfs <- list(classic = 2, HK = c(2, 8))
  expected <- c(classic = pchisq(4.2, 2, lower.tail = FALSE),
                HK = pf(3.1, 2, 8, lower.tail = FALSE))
  expect_equal(meta:::pval_subgroup(statistic, dfs), expected)
  expect_equal(meta:::pval_subgroup(rev(statistic), rev(dfs)), rev(expected))
  expect_equal(meta:::pval_subgroup(c(first = 4.2, second = 3.1), 2),
               c(first = pchisq(4.2, 2, lower.tail = FALSE),
                 second = pchisq(3.1, 2, lower.tail = FALSE)))
  for (denominator in c(6.42, Inf))
    expect_equal(meta:::pval_subgroup(c(HK = 3.1), c(2, denominator)),
                 c(HK = pf(3.1, 2, denominator, lower.tail = FALSE)))
})

test_that("subgroup P fallback rejects unavailable or invalid DF without guessing", {
  for (df in list(NULL, numeric(), NA_real_, 0, -1,
                  c(1, NA_real_), c(0, 8), c(1, 0), c(1, -2), c(1, 2, 3)))
    expect_identical(meta:::pval_subgroup(c(HK = 3), df), c(HK = NA_real_))
  expect_identical(meta:::pval_subgroup(c(HK = NA_real_), c(1, 8)),
                   c(HK = NA_real_))
  expect_identical(meta:::pval_subgroup(numeric(), list()), numeric())
  expect_error(meta:::pval_subgroup(c(classic = 3, HK = 4), c(1, 8)),
               "method-specific degrees of freedom", fixed = TRUE)
  expect_error(meta:::pval_subgroup(c(classic = 3, HK = 4), list(1)),
               "Number of subgroup-test DF entries", fixed = TRUE)
})

test_that("print and real forest recover absent common and random subgroup P", {
  x <- subgroup_inference_fixture(c("HK", "classic"))
  # Exercise a fractional F denominator through both presentation paths too.
  x$df.Q.b.random[[1]] <- c(1, 6.42)
  x$pval.Q.b.random[1] <- pf(x$Q.b.random[1], 1, 6.42, lower.tail = FALSE)
  expect_equal(unname(x$pval.Q.b.common),
               pchisq(x$Q.b.common, x$df.Q.b.common, lower.tail = FALSE))
  baseline.print <- subgroup_inference_print(x)
  baseline.forest <- subgroup_inference_forest(x)
  expect_length(baseline.forest, 2)
  expect_true(any(grepl("INFERENCE_COMMON", baseline.forest)))
  expect_true(any(grepl("INFERENCE_RANDOM", baseline.forest)))
  for (field in c("pval.Q.b.common", "pval.Q.b.random")) {
    missing <- x
    missing[[field]] <- NULL
    expect_identical(subgroup_inference_print(missing), baseline.print)
    expect_identical(subgroup_inference_forest(missing), baseline.forest)
    expect_null(missing[[field]])
  }
})

test_that("print and real forest preserve stored subgroup P including NA", {
  x <- subgroup_inference_fixture()
  for (field in c("pval.Q.b.common", "pval.Q.b.random")) {
    missing <- x
    missing[[field]] <- NULL
    fallback.print <- subgroup_inference_print(missing)
    fallback.forest <- subgroup_inference_forest(missing)
    tag <- if (field == "pval.Q.b.common") "INFERENCE_COMMON" else "INFERENCE_RANDOM"
    for (value in c(0.4321, NA_real_)) {
      stored <- x
      stored[[field]] <- value
      printed <- subgroup_inference_print(stored)
      drawn <- subgroup_inference_forest(stored)
      expect_false(identical(printed, fallback.print))
      expect_false(identical(drawn, fallback.forest))
      label <- drawn[grepl(tag, drawn, fixed = TRUE)]
      expect_length(label, 1)
      expect_match(label, if (is.na(value)) "NA" else "0.4321", fixed = TRUE)
      rows <- printed[grepl("^Between groups", printed)]
      row <- rows[if (field == "pval.Q.b.common") 1 else 2]
      expect_match(row, if (is.na(value)) "--" else "0.4321", fixed = TRUE)
      expect_identical(stored[[field]], value)
    }
  }
})

test_that("common forest annotations extract scalar DF from merged lists", {
  ma <- metagen(c(-0.3, 0.2, 0.1, 0.4, 0.7, 0.9, 0.6, 1.2, 0.8), rep(0.2, 9),
                 sm = "MD", subgroup = rep(c("A", "B", "C"), each = 3),
                 method.tau = "DL")
  merged <- metamerge(ma, ma)
  nested <- metamerge(merged, ma)
  observed <- NULL
  draw <- meta:::forest_meta_internal
  testthat::local_mocked_bindings(
    forest_meta_internal = function(...) {
      observed <<- mget(c("text.subgroup.common", "df.Q.b.common",
                          "hetstat.Q.b.common"), envir = parent.frame())
      draw(...)
    }, .package = "meta")
  grDevices::pdf(NULL, width = 20, height = 16)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (x in list(ma, merged, nested)) {
    before <- x
    expected <- unname(unlist(x$df.Q.b.common))
    for (layout in c("meta", "meta8", "RevMan5", "BMJ", "JAMA")) {
      expect_silent(forest(
        x, layout = layout, random = FALSE, details = FALSE,
        print.Q.subgroup = TRUE, test.subgroup.common = TRUE,
        overall.hetstat = FALSE, resid.hetstat = FALSE))
      expect_true(is.numeric(observed$df.Q.b.common))
      expect_equal(unname(observed$df.Q.b.common), expected)
      expect_length(observed$hetstat.Q.b.common, length(expected))
      label <- observed$text.subgroup.common
      expect_identical(label[[3]][[1]], as.name("^"))
      if (layout %in% c("meta", "meta8", "JAMA")) {
        df <- label[[3]][[2]][[3]]
        expect_true(is.numeric(df))
        expect_length(df, 1L)
        expect_equal(unname(df), expected[1])
      }
      else {
        expect_true(all(mapply(grepl, paste0("df ?= ?", expected, "$"),
                               observed$hetstat.Q.b.common)))
      }
    }
    expect_identical(x, before)
  }
})

test_that("random DF validation is independent of stored P and unused when hidden", {
  ma <- subgroup_inference_fixture(c("HK", "classic"))
  hidden <- subgroup_inference_print(ma, common = FALSE, test.subgroup = FALSE)
  for (df in list(list(), list(1), c(1, 8))) {
    error <- if (is.list(df)) "Number of subgroup-test DF entries" else
      "Multiple subgroup tests require a list"
    for (missing.p in c(FALSE, TRUE)) {
      x <- ma
      x$df.Q.b.random <- df
      if (missing.p)
        x$pval.Q.b.random <- NULL
      before <- x
      expect_error(subgroup_inference_print(x, common = FALSE), error, fixed = TRUE)
      expect_error(subgroup_inference_forest(x, common = FALSE), error, fixed = TRUE)
      expect_silent(output <- subgroup_inference_print(
        x, common = FALSE, test.subgroup = FALSE))
      expect_identical(output, hidden)
      expect_silent(labels <- subgroup_inference_forest(
        x, common = FALSE, test.subgroup = FALSE))
      expect_length(labels, 0L)
      expect_silent(subgroup_inference_print(x, random = FALSE))
      expect_silent(subgroup_inference_forest(x, random = FALSE))
      expect_identical(x, before)
    }
  }
})

test_that("zero random subgroup tests do not index an empty DF list", {
  x <- subgroup_inference_fixture(c("classic", "HK"))
  x$Q.b.random <- numeric()
  x$df.Q.b.random <- list()
  x$pval.Q.b.random <- NULL
  before <- x
  for (visible in c(TRUE, FALSE)) {
    expect_silent(output <- subgroup_inference_print(
      x, common = FALSE, test.subgroup = visible))
    expect_false(any(grepl("^Test for subgroup differences", output)))
    expect_silent(labels <- subgroup_inference_forest(
      x, common = FALSE, test.subgroup = visible))
    expect_length(labels, 0L)
  }
  expect_identical(x, before)
})

test_that("shared scalar random DF and missing-field fallbacks repeat by test", {
  ma <- metagen(c(0.1, 0.4, 0.2, 0.8, 1.2, 1.5, 0.9, 1.8), rep(0.2, 8),
                 sm = "MD", subgroup = rep(c("A", "B"), each = 4),
                 method.random.ci = c("HK", "classic"))
  printed <- subgroup_inference_print(ma, common = FALSE)
  drawn <- subgroup_inference_forest(ma, common = FALSE)
  for (df in list(ma$df.Q.b, NULL)) for (missing.p in c(FALSE, TRUE)) {
    x <- ma
    x$df.Q.b.random <- df
    if (missing.p)
      x$pval.Q.b.random <- NULL
    before <- x
    expect_identical(subgroup_inference_print(x, common = FALSE), printed)
    expect_identical(subgroup_inference_forest(x, common = FALSE), drawn)
    expect_identical(x, before)
  }
})

test_that("random F DF pairs retain repeated, fractional, and infinite values", {
  ma <- subgroup_inference_fixture()
  for (denominator in c(1, 6.42, Inf)) {
    x <- ma
    x$df.Q.b.random <- c(1, denominator)
    x$pval.Q.b.random <- pf(x$Q.b.random, 1, denominator, lower.tail = FALSE)
    printed <- subgroup_inference_print(x, common = FALSE)
    drawn <- subgroup_inference_forest(x, common = FALSE)
    df.label <- paste(1, denominator, sep = ", ")
    expect_true(any(grepl(df.label, printed, fixed = TRUE)))
    expect_length(drawn, 1L)
    expect_match(drawn, '"F"', fixed = TRUE)
    expect_match(drawn, df.label, fixed = TRUE)
    x$pval.Q.b.random <- NULL
    expect_identical(subgroup_inference_print(x, common = FALSE), printed)
    expect_identical(subgroup_inference_forest(x, common = FALSE), drawn)
  }
})

test_that("common DF validation is consistent and hidden metadata is ignored", {
  ma <- metagen(c(0.1, 0.4, 0.2, 0.8, 1.2, 1.5, 0.9, 1.8), rep(0.2, 8),
                 sm = "MD", subgroup = rep(c("A", "B"), each = 4))
  ma <- metamerge(ma, ma)
  hidden <- subgroup_inference_print(ma, random = FALSE, test.subgroup = FALSE)
  for (df in list(list(), list(1))) for (missing.p in c(FALSE, TRUE)) {
    x <- ma
    x$df.Q.b.common <- df
    if (missing.p)
      x$pval.Q.b.common <- NULL
    before <- x
    expect_error(subgroup_inference_print(x, random = FALSE),
                 "Number of subgroup-test DF entries", fixed = TRUE)
    expect_error(subgroup_inference_forest(x, random = FALSE),
                 "Number of subgroup-test DF entries", fixed = TRUE)
    expect_identical(subgroup_inference_print(
      x, random = FALSE, test.subgroup = FALSE), hidden)
    expect_length(subgroup_inference_forest(
      x, random = FALSE, test.subgroup = FALSE), 0L)
    expect_silent(subgroup_inference_print(x, common = FALSE))
    expect_silent(subgroup_inference_forest(x, common = FALSE))
    expect_identical(x, before)
  }
})

test_that("native clustered CR2 fractional subgroup DF reach print and summary", {
  skip_if_not_installed("clubSandwich")
  x <- metagen(
    c(0.1, 0.4, 0.2, 0.8, 1.2, 1.5, 0.9, 1.8),
    c(0.1, 0.2, 0.15, 0.25, 0.1, 0.2, 0.15, 0.25),
    subgroup = rep(c("A", "B"), each = 4), cluster = rep(1:4, 2),
    method.random.ci = "CR2", common = FALSE, tau.common = TRUE,
    method.tau.ci = "")
  df <- x$df.Q.b.random
  expect_length(df, 2L)
  expect_true(is.finite(df[2]))
  expect_gt(abs(df[2] - round(df[2])), 1e-4)
  expect_equal(unname(x$pval.Q.b.random),
               pf(x$Q.b.random, df[1], df[2], lower.tail = FALSE))
  before <- x
  printed <- subgroup_inference_print(x)
  summarized <- subgroup_inference_print(summary(x))
  label <- meta:::format_df(df, lab.NA = "NA")
  for (output in list(printed, summarized)) {
    expect_true(any(grepl(label, output, fixed = TRUE)))
    expect_true(any(grepl("F", output[grepl("d.f.", output, fixed = TRUE)],
                          fixed = TRUE)))
  }
  drawn <- subgroup_inference_forest(x)
  expect_length(drawn, 1L)
  expect_match(drawn, paste(df, collapse = ", "), fixed = TRUE)
  x$pval.Q.b.random <- NULL
  expect_identical(subgroup_inference_print(x), printed)
  expect_identical(subgroup_inference_print(summary(x)), summarized)
  expect_identical(subgroup_inference_forest(x), drawn)
  expect_identical(x[names(x) != "pval.Q.b.random"],
                   before[names(before) != "pval.Q.b.random"])
})

test_that("DF normalization counts subgroup statistics rather than CI methods", {
  x <- subgroup_inference_fixture(c("HK", "classic"))
  x$Q.b.random <- x$Q.b.random[1]
  x$df.Q.b.random <- x$df.Q.b.random[[1]]
  x$pval.Q.b.random <- x$pval.Q.b.random[1]
  before <- x
  printed <- subgroup_inference_print(x, common = FALSE)
  expect_equal(sum(grepl("^Test for subgroup differences", printed)), 1)
  expect_equal(sum(grepl("^Results for subgroups", printed)), 2)
  drawn <- subgroup_inference_forest(x, common = FALSE)
  expect_length(drawn, 1L)
  expect_match(drawn, '"F"', fixed = TRUE)
  expect_identical(x, before)
  x$pval.Q.b.random <- NULL
  expect_identical(subgroup_inference_print(x, common = FALSE), printed)
  expect_identical(subgroup_inference_forest(x, common = FALSE), drawn)
})

test_that("missing DF entries stay unavailable without flattening adjacent pairs", {
  statistic <- c(unavailable = 3, fractional = 4, infinite = 5)
  for (missing in list(NULL, numeric(), NA_real_)) {
    df <- list(missing, c(1, 6.42), c(1, Inf))
    before <- df
    expect_identical(meta:::subgroup_test_df(df, 3L), df)
    expect_equal(meta:::pval_subgroup(statistic, df),
                 c(unavailable = NA_real_,
                   fractional = pf(4, 1, 6.42, lower.tail = FALSE),
                   infinite = pf(5, 1, Inf, lower.tail = FALSE)))
    expect_identical(df, before)
  }
})
