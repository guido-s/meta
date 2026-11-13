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
