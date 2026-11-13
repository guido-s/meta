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
