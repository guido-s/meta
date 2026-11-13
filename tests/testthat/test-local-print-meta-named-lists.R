test_that("meta named settings preserve legacy formatting and subgroup modes", {
  ma <- metagen(c(0.12345, 0.65432, -0.23456, 0.45678),
                c(0.1, 0.2, 0.3, 0.15), sm = "MD", method.tau = "DL",
                subgroup = c("A", "A", "B", "B"), prediction = TRUE,
                method.predict = "V")
  flags <- list(subgroup.name = FALSE, CMH = FALSE,
                tau2 = TRUE, tau2.ci = TRUE, tau = TRUE, tau.ci = TRUE,
                Q = TRUE, I2 = TRUE, I2.ci = TRUE, H = TRUE, Rb = TRUE)
  precision <- list(digits = 2, stat = 3, pval = 3, tau2 = 3, tau = 2,
                    Q = 2, df = 1, pval.Q = 3, H = 2, I2 = 2)
  for (overall in c(FALSE, TRUE)) {
    for (random in c(FALSE, TRUE)) {
      extra <- list(overall = overall, common = !random, random = random,
                    prediction.subgroup = c(TRUE, FALSE),
                    details.methods = FALSE, text.tau2 = "variance",
                    text.tau = "SD", text.I2 = "I-squared")
      grouped <- capture.output(do.call(print, c(
        list(ma, print = flags, digits = precision), extra)))
      legacy <- capture.output(do.call(print, c(
        list(ma, digits = precision$digits), extra,
        setNames(flags, paste0("print.", names(flags))),
        setNames(precision[-1], paste0("digits.", names(precision)[-1])))))
      expect_identical(grouped, legacy)
    }
  }
  expect_identical(
    capture.output(print(ma, digits = 2, print.tau = FALSE)),
    capture.output(print(ma, digits = list(digits = 2),
                         print = list(tau = FALSE))))
})

test_that("summary named settings reach individual and pooled results", {
  ma <- metagen(c(0.12345, NA, -0.23456, 0.45678),
                c(0.1, 0.2, 0.3, 0.15), sm = "MD", method.tau = "DL")
  flags <- list(tau2 = TRUE, tau2.ci = TRUE, tau = FALSE, tau.ci = FALSE,
                Q = TRUE, I2 = TRUE, I2.ci = TRUE, H = FALSE, Rb = FALSE)
  precision <- list(digits = 2, se = 3, stat = 1, pval = 3,
                    tau2 = 2, tau = 3, Q = 1, df = 1, pval.Q = 2,
                    H = 1, I2 = 2, prop = 2, weight = 2)
  subgroup <- metagen(c(0.12345, 0.65432, -0.23456, 0.45678),
                      c(0.1, 0.2, 0.3, 0.15), sm = "MD", method.tau = "DL",
                      subgroup = c("A", "A", "B", "B"))
  binary <- metabin(c(2, 4, 6), rep(20, 3), c(1, 3, 4), rep(20, 3),
                    method.tau = "DL")
  for (model in list(ma, metagen(0.12345, 0.1, sm = "MD"), subgroup, binary)) {
    for (ma.results in c(FALSE, TRUE)) {
      extra <- list(details = TRUE, ma = ma.results, details.methods = FALSE)
      grouped <- capture.output(do.call(print, c(
        list(summary(model), print = flags, digits = precision), extra)))
      legacy <- capture.output(do.call(print, c(
        list(summary(model), digits = precision$digits), extra,
        setNames(flags, paste0("print.", names(flags))),
        setNames(precision[-1], paste0("digits.", names(precision)[-1])))))
      expect_identical(grouped, legacy)
    }
  }
  pooled <- capture.output(print(ma, digits = list(stat = 1, pval = 3),
                                 details.methods = FALSE, header = FALSE))
  detailed <- capture.output(print(summary(ma),
                                   digits = list(stat = 1, pval = 3),
                                   details.methods = FALSE))
  expect_true(all(pooled[grepl("effect model", pooled)] %in% detailed))
  expect_identical(capture.output(print(summary(ma), digits = 2)),
                   capture.output(print(summary(ma), digits = list(digits = 2))))
  expect_identical(
    capture.output(print(summary(subgroup), print.byvar = FALSE)),
    capture.output(print(summary(subgroup), print = list(subgroup.name = FALSE))))
  for (show in c(FALSE, TRUE)) {
    expect_identical(
      capture.output(print(binary, print.CMH = show)),
      capture.output(print(binary, print = list(CMH = show))))
  }
})

test_that("cumulative and influence printers forward lists and nonformat settings", {
  ma <- metagen(c(0.12345, 0.65432, -0.23456, 0.45678),
                c(0.1, 0.2, 0.3, 0.15), sm = "MD", method.tau = "DL",
                cid = c(-0.1, 0.1), prediction = TRUE, method.predict = "V")
  flags <- list(stat = TRUE, tau2 = TRUE, tau2.ci = TRUE,
                tau = TRUE, tau.ci = TRUE, I2 = TRUE, I2.ci = TRUE,
                prob = FALSE)
  precision <- list(digits = 2, stat = 1, pval = 3, tau2 = 3,
                    tau = 3, I2 = 2, cid = 3, percent = 2)
  for (model in list(metacum(ma), metainf(ma))) {
    for (overall in c(FALSE, TRUE)) {
      extra <- list(overall = overall, prediction = TRUE, header = FALSE,
                    backtransf = FALSE, lab.NA = "missing",
                    text.tau2 = "variance", text.tau = "SD",
                    text.I2 = "I-squared", details.methods = TRUE)
      grouped <- capture.output(do.call(print, c(
        list(model, print = flags, digits = precision), extra)))
      legacy <- capture.output(do.call(print, c(
        list(model, digits = precision$digits), extra,
        setNames(flags, paste0("print.", names(flags))),
        setNames(precision[-1], paste0("digits.", names(precision)[-1])))))
      expect_identical(grouped, legacy)
      expect_true(any(grepl("variance", grouped, fixed = TRUE)))
    }
    expect_identical(capture.output(print(model, digits = 2)),
                     capture.output(print(model, digits = list(digits = 2))))
  }
  influence <- metainf(ma)
  influence$TE[1] <- NA
  output <- capture.output(print(influence, lab.NA = "missing"))
  expect_true(any(grepl("missing", output, fixed = TRUE)))
  expect_identical(capture.output(print(influence, lab.NA = "missing")),
                   capture.output(meta:::print.metacum(influence,
                                                       lab.NA = "missing")))
})

test_that("named defaults retain method-specific settings", {
  old <- settings.meta(digits = 3, digits.forest = 1,
                       digits.pval = 1, digits.pval.Q = 1)
  withr::defer(settings.meta(old))
  ma <- metagen(c(0.12345, 0.65432, -0.23456, 0.45678),
                c(0.1, 0.2, 0.3, 0.15), sm = "MD", method.tau = "DL")
  for (model in list(ma, summary(ma), metacum(ma), metainf(ma))) {
    cl <- if (inherits(model, "summary.meta")) "summary.meta" else
      if (inherits(model, "metacum")) "metacum" else
        if (inherits(model, "metainf")) "metainf" else "meta"
    x <- if (cl == "summary.meta") ma else model
    flags <- meta:::print_defaults(x, cl, func = "print")
    precision <- meta:::digits_defaults(cl, func = "print")
    expect_identical(capture.output(print(model)),
                     capture.output(print(model, print = flags,
                                          digits = precision)))
    expect_identical(capture.output(print(model)),
                     capture.output(print(model, print = list(),
                                          digits = list(digits = meta:::gs("digits")))))
  }
  expect_true(meta:::print_defaults(ma, "metacum", "print")$prob)
  expect_false(meta:::print_defaults(ma, "metacum", "forest")$prob)
  expect_equal(meta:::digits_defaults("meta", "print")$pval, 2)
  expect_equal(meta:::digits_defaults("metacum", "print")$pval, 1)
  expect_equal(meta:::digits_defaults("metainf", "print")$digits, 3)
  expect_equal(meta:::digits_defaults("metainf", "forest")$digits, 1)
  expect_named(meta:::digits_defaults("summary.meta", "print"),
               c("digits", "stat", "pval", "tau2", "tau", "Q", "df",
                 "pval.Q", "H", "I2", "se", "prop", "weight"))
})

test_that("formatting duplicates and invalid values are rejected", {
  ma <- metagen(c(0.1, 0.3, 0.5), rep(0.1, 3), sm = "MD")
  for (model in list(ma, summary(ma), metacum(ma), metainf(ma))) {
    expect_error(print(model, print = list(tau = FALSE), print.tau = TRUE),
                 "cannot both be specified")
    expect_error(print(model, digits = list(pval = 3), digits.pval = 2),
                 "cannot both be specified")
    expect_error(print(model, print = list(tau = TRUE, tau = FALSE)),
                 "duplicated")
    expect_error(print(model, digits = list(pval = 3, pval = 2)), "duplicated")
    expect_error(print(model, print.tau = TRUE, print.tau = FALSE), "duplicated")
    expect_error(print(model, digits.pval = 3, digits.pval = 2), "duplicated")
    expect_error(print(model, print = TRUE), "named list")
    expect_error(print(model, digits = list(2)), "named list")
    expect_error(print(model, print = list(tau = NA)), "must be a logical")
    expect_error(print(model, digits = list(pval = -1)), "must be")
  }
})
