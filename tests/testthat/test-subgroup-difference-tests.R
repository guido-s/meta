subgroup_difference_fixture <- function(methods = "HK", method = "GLMM") {
  dat <- data.frame(event = c(10, 20, 30, 40, 50, 30, 40, 50, 60, 70),
                    total = 100, group = rep(c("A", "B"), each = 5))
  metaprop(event, total, data = dat, method = method, method.random.ci = methods,
             common = FALSE, subgroup = group, tau.common = TRUE)
}

subgroup_difference_forest_label <- function(x, layout,
                                             tag = "SUBGROUP_DIFFERENCE", ...) {
  observed <- NULL
  draw <- meta:::forest_meta_internal
  testthat::local_mocked_bindings(
    forest_meta_internal = function(...) {
      args <- list(...)
      observed <<- lapply(args$cols$col.studlab$labels, `[[`, "label")
      do.call(draw, args)
    }, .package = "meta")
  grDevices::pdf(NULL, width = 16, height = 12)
  on.exit(grDevices::dev.off())
  forest(x, layout = layout, label.test.subgroup.random = "SUBGROUP_DIFFERENCE",
         resid.hetlab = "RESIDUAL_HETEROGENEITY",
         print.Q.subgroup = TRUE, digits.Q = 2, digits.pval.Q = 4,
         scientific.pval = FALSE, zero.pval = TRUE, JAMA.pval = FALSE, ...)
  Filter(function(label) is.call(label) &&
           identical(label[[2]], tag), observed)
}

subgroup_residual_fixture <- function(methods = "classic", clustered = FALSE) {
  args <- list(
    TE = c(-0.4, -0.1, 0.2, 0.1, 0.6, 0.3, 0.5, 0.2, 0.9, 0.8, 1.1, 0.6),
    seTE = rep(0.2, 12), studlab = paste0("S", 1:12),
    subgroup = rep(c("A", "B"), each = 6), common = FALSE,
    tau.common = TRUE, method.tau = "REML", method.random.ci = methods)
  if (clustered)
    args$cluster <- rep(1:6, each = 2)
  do.call(metagen, args)
}

expect_subgroup_residual_print <- function(x) {
  output <- capture.output(print(
    x, details.methods = FALSE, digits.Q = 2, digits.pval.Q = 4,
    scientific.pval = FALSE, zero.pval = TRUE, JAMA.pval = FALSE))
  residual <- which(output == "Test of residual heterogeneity:")
  expect_length(residual, 1)
  expected <- paste(meta:::formatN(x$Q.resid, digits = 2),
                     meta:::formatN(x$df.Q.resid, digits = 0),
                     meta:::formatPT(x$pval.Q.resid, digits = 4,
                                    scientific = FALSE, zero = TRUE, JAMA = FALSE))
  expect_identical(strsplit(trimws(output[residual + 2]), "[[:space:]]+")[[1]],
                   strsplit(trimws(expected), "[[:space:]]+")[[1]])
  expect_false(any(grepl("^Within groups", output)))
  expect_equal(sum(grepl("^Test for subgroup differences", output)),
                length(x$method.random.ci))
  summary.output <- capture.output(print(summary(x), details.methods = FALSE))
  expect_equal(sum(summary.output == "Test of residual heterogeneity:"), 1)
  expect_false(any(grepl("^Within groups", summary.output)))
}

test_that("single HK subgroup tests print F with numerator and denominator DF", {
  ma <- subgroup_difference_fixture()
  expect_equal(unname(ma$df.Q.b.random), c(1, 8))
  output <- capture.output(print(ma, details.methods = FALSE,
                                  digits.Q = 2, digits.pval.Q = 4,
                                  scientific.pval = FALSE, JAMA.pval = FALSE))
  test <- grep("^Test for subgroup differences", output)
  expect_length(test, 1)
  expect_match(output[test + 1], "F", fixed = TRUE)
  expect_match(output[test + 2], "4.86", fixed = TRUE)
  expect_match(output[test + 2], "1, 8", fixed = TRUE)
  expect_match(output[test + 2], "0.0587", fixed = TRUE)
})

test_that("forest subgroup F labels retain both DF in all layouts", {
  ma <- subgroup_difference_fixture()
  for (layout in c("RevMan5", "meta", "BMJ", "JAMA")) {
    labels <- subgroup_difference_forest_label(ma, layout)
    expect_length(labels, 1)
    label <- labels[[1]]
    if (layout %in% c("RevMan5", "BMJ")) {
      expect_identical(label[[3]], "F")
      expect_match(label[[4]], "df = 1, 8", fixed = TRUE)
    }
    else {
      expect_identical(label[[3]][[2]], quote(F))
      expect_identical(label[[3]][[3]], "1, 8")
    }
    text <- paste(deparse(label), collapse = " ")
    expect_match(text, "4.86", fixed = TRUE)
    expect_match(text, "0.0587", fixed = TRUE)
  }
})

test_that("multiple random subgroup tests keep method-specific DF vectors", {
  ma <- subgroup_difference_fixture(c("classic", "HK"))
  expect_length(ma$df.Q.b.random, 2)
  expect_length(ma$df.Q.b.random[[1]], 1)
  expect_equal(unname(ma$df.Q.b.random[[2]]), c(1, 8))
  output <- capture.output(print(ma, details.methods = FALSE))
  tests <- grep("^Test for subgroup differences", output)
  expect_length(tests, 2)
  expect_match(output[tests[1] + 1], "Q", fixed = TRUE)
  expect_match(output[tests[2] + 1], "F", fixed = TRUE)
  expect_match(output[tests[2] + 2], "1, 8", fixed = TRUE)
  expect_false(grepl("1, 8", output[tests[1] + 2], fixed = TRUE))

  label <- subgroup_difference_forest_label(ma, "RevMan5")[[1]]
  expect_equal(label[[3]], quote("Chi"^2))
  expect_false(grepl("df = 1, 8", label[[4]], fixed = TRUE))
  reversed <- subgroup_difference_fixture(c("HK", "classic"))
  label <- subgroup_difference_forest_label(reversed, "RevMan5")[[1]]
  expect_identical(label[[3]], "F")
  expect_match(label[[4]], "df = 1, 8", fixed = TRUE)
})

test_that("random subgroup tests omit within rows and report residual heterogeneity once", {
  methods <- gs("meth4random.ci")
  for (selected in list(methods, rev(methods))) {
    ma <- subgroup_difference_fixture(selected, method = "Inverse")
    fields <- c("Q.w.random", "df.Q.w.random", "pval.Q.w.random",
                 "Q.resid", "df.Q.resid", "pval.Q.resid")
    stored <- ma[fields]
    expect_length(ma$Q.w.random, 1)
    expect_length(ma$pval.Q.w.random, 1)
    output <- capture.output(print(ma, details.methods = FALSE,
                                    digits.Q = 2, digits.pval.Q = 4,
                                    scientific.pval = FALSE, JAMA.pval = FALSE))
    expect_false(any(grepl("^Within groups", output)))
    tests <- grep("^Test for subgroup differences", output)
    expect_length(tests, length(selected))
    expect_equal(sum(grepl("^Between groups", output)), length(selected))
    expect_true(all(grepl("^Between groups", output[tests + 2])))
    residual <- which(output == "Test of residual heterogeneity:")
    expect_length(residual, 1)
    expect_identical(strsplit(trimws(output[residual + 2]), "[[:space:]]+")[[1]],
                     c("80.45", "8", "<", "0.0001"))
    expect_identical(ma[fields], stored)

    summary.output <- capture.output(print(summary(ma), details.methods = FALSE))
    expect_false(any(grepl("^Within groups", summary.output)))
    expect_equal(sum(summary.output == "Test of residual heterogeneity:"), 1)
  }
})

test_that("common-effect decomposition retains its within-group row", {
  ma <- update(subgroup_difference_fixture(c("classic", "HK"), method = "Inverse"),
                 common = TRUE)
  output <- capture.output(print(ma, details.methods = FALSE))
  tests <- grep("^Test for subgroup differences", output)
  expect_length(tests, 3)
  expect_equal(sum(grepl("^Within groups", output)), 1)
  expect_true(grepl("^Within groups", output[tests[1] + 3]))
  expect_false(any(grepl("^Within groups", output[seq.int(tests[2], length(output))])))
  expect_equal(sum(output == "Test of residual heterogeneity:"), 1)
})

test_that("residual heterogeneity test respects controls and canonical availability", {
  ma <- subgroup_difference_fixture(method = "Inverse")
  for (args in list(list(print.Q = FALSE), list(overall.hetstat = FALSE))) {
    output <- capture.output(do.call(print, c(list(x = ma, details.methods = FALSE), args)))
    expect_false(any(output == "Test of residual heterogeneity:"))
  }
  for (missing in c("Q.resid", "df.Q.resid", "zero DF")) {
    unavailable <- ma
    if (missing == "zero DF")
      unavailable$df.Q.resid <- 0
    else
      unavailable[[missing]] <- NULL
    output <- capture.output(print(unavailable, details.methods = FALSE))
    expect_false(any(output == "Test of residual heterogeneity:"))
    expect_identical(unavailable$Q.w.random, ma$Q.w.random)
    expect_identical(unavailable$pval.Q.w.random, ma$pval.Q.w.random)
  }

  glmm <- subgroup_difference_fixture()
  expect_true(is.na(glmm$Q.w.random))
  expect_true(is.finite(glmm$Q.resid))
  output <- capture.output(print(glmm, details.methods = FALSE))
  expect_equal(sum(output == "Test of residual heterogeneity:"), 1)
})

test_that("forest residual annotation uses canonical rather than legacy fields", {
  for (method in c("Inverse", "GLMM")) {
    ma <- subgroup_difference_fixture(method = method)
    baseline <- subgroup_difference_forest_label(
      ma, "RevMan5", tag = "RESIDUAL_HETEROGENEITY",
      resid.hetstat = TRUE, print.Q = TRUE, print.pval.Q = TRUE)
    expect_length(baseline, 1)
    expect_match(paste(deparse(baseline[[1]]), collapse = " "), "80.45", fixed = TRUE)

    legacy <- c("Q.w.random", "df.Q.w.random", "pval.Q.w.random")
    for (remove in c(FALSE, TRUE)) {
      changed <- ma
      if (remove)
        changed[legacy] <- NULL
      else {
        changed$Q.w.random <- -999
        changed$df.Q.w.random <- 777
        changed$pval.Q.w.random <- 0.9999
      }
      label <- subgroup_difference_forest_label(
        changed, "RevMan5", tag = "RESIDUAL_HETEROGENEITY",
        resid.hetstat = TRUE, print.Q = TRUE, print.pval.Q = TRUE)
      expect_identical(label, baseline)
    }

    changed <- ma
    changed$Q.resid <- 12.34
    changed$df.Q.resid <- 7
    changed$pval.Q.resid <- 0.0123
    label <- subgroup_difference_forest_label(
      changed, "RevMan5", tag = "RESIDUAL_HETEROGENEITY",
      resid.hetstat = TRUE, print.Q = TRUE, print.pval.Q = TRUE)
    expect_length(label, 1)
    text <- paste(deparse(label[[1]]), collapse = " ")
    for (value in c("12.34", "df = 7", "0.0123"))
      expect_match(text, value, fixed = TRUE)
    expect_false(identical(label, baseline))
  }
})

test_that("residual reporting covers univariate, three-level, and robust models", {
  for (clustered in c(FALSE, TRUE)) {
    for (methods in list(c("classic", "HK"), c("CR0", "CR1"))) {
      ma <- subgroup_residual_fixture(methods, clustered)
      expect_identical(ma$three.level, clustered)
      expect_true(is.finite(ma$Q.resid))
      expect_equal(ma$df.Q.resid, 10)
      fields <- c("Q.w.random", "df.Q.w.random", "pval.Q.w.random",
                   "Q.resid", "df.Q.resid", "pval.Q.resid")
      stored <- ma[fields]
      expect_subgroup_residual_print(ma)
      expect_identical(ma[fields], stored)
    }
  }
})

test_that("CR2 residual reporting is separate from robust subgroup inference", {
  skip_if_not_installed("clubSandwich")
  for (clustered in c(FALSE, TRUE)) {
    ma <- subgroup_residual_fixture("CR2", clustered)
    expect_true(is.finite(ma$Q.resid))
    expect_true(is.finite(ma$pval.Q.resid))
    expect_subgroup_residual_print(ma)
  }
})

test_that("separate or preset tau does not invent a canonical residual test", {
  ma <- subgroup_difference_fixture(c("classic", "HK"), method = "Inverse")
  for (args in list(list(tau.common = FALSE), list(tau.preset = 0.5))) {
    changed <- do.call(update, c(list(object = ma), args))
    expect_true(all(is.na(changed$Q.resid)))
    output <- capture.output(print(changed, details.methods = FALSE))
    expect_false(any(output == "Test of residual heterogeneity:"))
    expect_false(any(grepl("^Within groups", output)))
    expect_equal(sum(grepl("^Test for subgroup differences", output)), 2)

    changed$Q.w.random <- 123
    changed$df.Q.w.random <- 8
    changed$pval.Q.w.random <- 0.001
    stored <- changed[c("Q.w.random", "df.Q.w.random", "pval.Q.w.random")]
    output <- capture.output(print(changed, details.methods = FALSE))
    expect_false(any(output == "Test of residual heterogeneity:"))
    expect_false(any(grepl("^Within groups", output)))
    expect_identical(changed[names(stored)], stored)
    summary.output <- capture.output(print(summary(changed), details.methods = FALSE))
    expect_false(any(summary.output == "Test of residual heterogeneity:"))
  }
})
