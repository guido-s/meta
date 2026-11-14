test_that("grouped and individual print/digits options give identical output", {
  ma <- metagen(c(0.23456, 0.34567, 0.45678), c(0.12345, 0.23456, 0.34567),
                 sm = "MD", prediction = TRUE)
  mb <- metabind(ma, update(ma, method.random.ci = "HK", method.predict = "HK",
                            adhoc.hakn.pi = "se"))
  flags <- list(func = TRUE, method = TRUE, method.interval = TRUE,
                sm = TRUE, se = TRUE, pval = TRUE, stat = TRUE,
                level.ma = TRUE, level.predict = TRUE, level.hetstat = TRUE,
                tau2 = TRUE, tau2.ci = TRUE, tau = TRUE, tau.ci = TRUE,
                I2 = TRUE, I2.ci = TRUE, Q = TRUE, pval.Q = TRUE,
                H = TRUE, H.ci = TRUE)
  precision <- list(digits = 4, se = 1, stat = 3, pval = 1, pval.Q = 3,
                    Q = 1, tau2 = 3, tau = 2, I2 = 2)

  grouped <- capture.output(print(mb, print = flags, digits = precision,
                                  details.methods = FALSE))
  individual <- c(
    list(x = mb, digits = precision$digits, details.methods = FALSE),
    setNames(flags, paste0("print.", names(flags))),
    setNames(precision[-1], paste0("digits.", names(precision)[-1]))
  )
  legacy <- capture.output(do.call(print, individual))
  expect_identical(grouped, legacy)
  expect_true(any(grepl(formatC(ma$TE.common, digits = 4, format = "f"),
                        grouped, fixed = TRUE)))
})

test_that("automatic confidence-level selection stays independent by table", {
  ma <- metagen(1:3, 3:1, sm = "MD", common = TRUE, random = FALSE,
                 level.ma = 0.95)
  other <- metagen(1:3, 3:1, sm = "MD", common = FALSE, random = TRUE,
                    level.ma = 0.99)
  third <- metagen(1:3, 3:1, sm = "MD", common = FALSE, random = TRUE,
                    level.ma = 0.90)
  output <- capture.output(print(metabind(ma, other, third),
                                  details.methods = FALSE))
  common <- which(output == "Common effect model")
  random <- which(output == "Random effects model")
  heterogeneity <- which(output == "Information on heterogeneity")
  expect_false(any(grepl("level", output[common:(random - 1)], fixed = TRUE)))
  expect_true(any(grepl("level", output[random:(heterogeneity - 1)],
                       fixed = TRUE)))
})

test_that("shared pooled results use one interval column for selected CI and PI levels", {
  dat <- data.frame(TE = c(0.1, 0.3, 0.5, 0.2, 0.4, 0.6), seTE = rep(0.1, 6),
                    first = rep(c("A", "B"), each = 3),
                    second = rep(c("C", "D"), 3))
  for (level in c(0.95, 0.90)) {
    ma <- metagen(TE, seTE, data = dat, subgroup = first, sm = "MD",
                    level.ma = 0.95, level.predict = level,
                    prediction = TRUE, method.predict = "V", method.tau = "DL")
    mb <- metabind(ma, update(ma, subgroup = second))
    output <- capture.output(print(mb, details.methods = FALSE))
    overall <- output[seq_len(match("Subgroup meta-analysis results", output) - 1)]
    header <- overall[grepl("CI/PI", overall, fixed = TRUE)]
    expect_length(header, 1)
    expect_false(any(grepl("90% PI", header, fixed = TRUE)))
    if (level == 0.95) {
      expect_match(header, "95% CI/PI", fixed = TRUE)
      expect_false(grepl("level", header, fixed = TRUE))
    }
    else {
      expect_false(grepl("95% CI/PI", header, fixed = TRUE))
      expect_true(grepl("level", header, fixed = TRUE))
    }

    output <- capture.output(print(mb, prediction = FALSE, details.methods = FALSE))
    overall <- output[seq_len(match("Subgroup meta-analysis results", output) - 1)]
    expect_true(any(grepl("95% CI", overall, fixed = TRUE)))
    expect_false(any(grepl("CI/PI", overall, fixed = TRUE)))

    output <- capture.output(print(mb, common = FALSE, random = FALSE,
                                    prediction = TRUE, details.methods = FALSE))
    overall <- output[seq_len(match("Subgroup meta-analysis results", output) - 1)]
    expect_true(any(grepl(meta:::percent_label(level, type = "PI"),
                          overall, fixed = TRUE)))
    expect_false(any(grepl("CI/PI", overall, fixed = TRUE)))
  }
})

test_that("single-subgroup printing works for both effect models", {
  ma <- metagen(1:3, rep(0.1, 3), subgroup = rep("A", 3), sm = "MD")
  mb <- metabind(ma, ma)
  for (random in c(FALSE, TRUE)) {
    expect_error(output <- capture.output(print(
      mb, common = !random, random = random, details.methods = FALSE)), NA)
    expect_true(any(grepl("d.f.", output, fixed = TRUE)))
  }
})

test_that("printing multiple CI methods retains their separate subgroup tests", {
  ma <- metagen(c(0.1, 0.3, 0.5, 0.2, 0.4, 0.6), rep(0.1, 6),
                  subgroup = rep(c("A", "B"), each = 3), sm = "MD",
                  method.tau = "DL", method.random.ci = c("classic", "HK"))
  mb <- metabind(ma, ma)
  output <- capture.output(print(
    mb, common = FALSE, print = list(method.interval = FALSE),
    digits = list(Q = 4), details.methods = FALSE))
  for (Q in as.vector(ma$Q.b.random))
    expect_true(any(grepl(trimws(meta:::formatN(Q, digits = 4)), output, fixed = TRUE)))

  output <- capture.output(print(
    mb, common = FALSE, print = list(pval.Q.b = FALSE), details.methods = FALSE))
  expect_false(any(grepl("Q.b|df.Q.b|pval.Q.b|d.f.", output)))
})

test_that("subgroup prediction degrees of freedom honor the missing-value label", {
  ma <- metagen(c(0.1, 0.3, 0.5), rep(0.1, 3), sm = "RR",
                  subgroup = c("A", "A", "B"),
                  prediction = TRUE, method.predict = "V")
  mb <- metabind(ma, ma, common = FALSE, random = FALSE, prediction = TRUE)
  output <- capture.output(print(mb, lab.NA = "MISSING", details.methods = FALSE))
  expect_true(any(grepl("B.*MISSING.*MISSING", output)))
  expect_false(any(grepl("<NA>", output, fixed = TRUE)))
})

test_that("subgroup level formatting follows heading calculation and handles missing levels", {
  ma <- metagen(c(0.1, 0.3, 0.5, 0.2, 0.4, 0.6), rep(0.1, 6),
                  subgroup = rep(c("A", "B"), each = 3), sm = "MD")
  mb <- metabind(ma, ma)
  mb$overall$level.ma[2] <- NA_real_
  for (random in c(FALSE, TRUE)) {
    output <- capture.output(print(
      mb, common = !random, random = random, overall = FALSE,
      print = list(level.ma = TRUE), lab.NA = "MISSING", details.methods = FALSE))
    expect_true(any(grepl("95% CI", output, fixed = TRUE)))
    expect_true(any(grepl("0.9500", output, fixed = TRUE)))
    expect_true(any(grepl("MISSING", output, fixed = TRUE)))
    expect_false(any(grepl("<NA>", output, fixed = TRUE)))
  }
})

test_that("common-only printing shares effects independently of heterogeneity", {
  dat <- data.frame(TE = c(0.1, 0.3, 0.5, 0.2, 0.4, 0.6), seTE = rep(0.1, 6),
                    first = rep(c("A", "B"), each = 3),
                    second = rep(c("C", "D"), 3))
  ma <- metagen(TE, seTE, data = dat, subgroup = first, sm = "MD", method.tau = "DL")
  for (different.heterogeneity in c(FALSE, TRUE)) {
    other <- update(ma, subgroup = second, method.random.ci = "HK",
                      method.tau = if (different.heterogeneity) "REML" else "DL")
    mb <- metabind(ma, other)
    output <- capture.output(print(
      mb, random = FALSE, print = list(tau2 = TRUE, tau = TRUE),
      details.methods = FALSE))
    expect_true(any(output == "Overall results"))
    opening <- output[seq_len(match("Subgroup meta-analysis results", output) - 1)]
    expect_equal(sum(grepl("Common effect model", opening, fixed = TRUE)), 1)
    expect_identical(any(grepl("tau^2", opening, fixed = TRUE)),
                      !different.heterogeneity)

    output <- capture.output(print(mb, random = FALSE, overall = FALSE,
                                    details.methods = FALSE))
    expect_false(any(output == "Overall results"))
    expect_true(any(output == "Subgroup meta-analysis results"))

    output <- capture.output(print(mb, random = FALSE, overall = NULL,
                                    details.methods = FALSE))
    expect_false(any(output == "Overall results"))
  }
})

test_that("print sharing uses only selected model methods and interval levels", {
  dat <- data.frame(TE = c(0.1, 0.3, 0.5, 0.2, 0.4, 0.6), seTE = rep(0.1, 6),
                    first = rep(c("A", "B"), each = 3),
                    second = rep(c("C", "D"), 3))
  ma <- metagen(TE, seTE, data = dat, subgroup = first, sm = "MD",
                  method.tau = "DL", prediction = TRUE, method.predict = "V")
  for (difference in c("random CI", "prediction level")) {
    other <- if (difference == "random CI")
      update(ma, subgroup = second, method.random.ci = "HK")
    else
      update(ma, subgroup = second, level.predict = 0.90)
    mb <- metabind(ma, other)
    output <- capture.output(print(mb, random = FALSE, prediction = FALSE,
                                    details.methods = FALSE))
    expect_true(any(output == "Overall results"), info = difference)

    selected <- capture.output(print(
      mb, random = difference == "random CI",
      prediction = difference == "prediction level", details.methods = FALSE))
    expect_false(any(selected == "Overall results"), info = difference)
  }

  mb <- metabind(ma, update(ma, subgroup = second, level.ma = 0.90))
  output <- capture.output(print(mb, common = FALSE, random = FALSE,
                                  prediction = TRUE, details.methods = FALSE))
  expect_true(any(output == "Overall results"))
  opening <- output[seq_len(match("Subgroup meta-analysis results", output) - 1)]
  expect_true(any(grepl("95% PI", opening, fixed = TRUE)))
  expect_equal(sum(grepl("Prediction interval", opening, fixed = TRUE)), 1)
})

test_that("subgroup analysis names print by default and respect explicit overrides", {
  dat <- data.frame(TE = c(0.1, 0.3, 0.5, 0.2, 0.4, 0.6), seTE = rep(0.1, 6),
                    age = rep(c(55, 65), each = 3),
                    region = rep(c("Europe", "Asia"), 3))
  ma <- metagen(TE, seTE, data = dat, subgroup = age, subgroup.name = "Age group",
                  sm = "MD", prediction = TRUE, method.predict = "V")
  mb <- metabind(ma, update(ma, subgroup = region, subgroup.name = "Region"))
  default <- capture.output(print(mb, details.methods = FALSE))
  expect_equal(sum(grepl("1. Age group", default, fixed = TRUE)), 3)
  expect_equal(sum(grepl("2. Region", default, fixed = TRUE)), 3)

  for (show in list(NULL, TRUE)) {
    explicit <- capture.output(print(mb, print = list(subgroup.name = show),
                                      details.methods = FALSE))
    expect_identical(explicit, default)
  }

  hidden <- capture.output(print(mb, print = list(subgroup.name = FALSE),
                                  details.methods = FALSE))
  expect_false(any(grepl("Age group|Region", hidden)))
  expect_equal(sum(hidden %in% c("1.", "2.")), 6)
  legacy <- capture.output(print(mb, print.subgroup.name = FALSE,
                                  details.methods = FALSE))
  expect_identical(legacy, hidden)
})
