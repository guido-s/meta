test_that("print and forest pass full option lists to catmeth_metabind", {
  expect_identical(names(formals(meta:::catmeth_metabind))[5:7],
                    c("backtransf", "print", "digits"))
  ma <- metagen(log(2:4), c(0.1, 0.2, 0.3), sm = "RR")
  mb <- metabind(ma, ma)
  calls <- list()
  file <- tempfile(fileext = ".pdf")
  on.exit(unlink(file))

  testthat::with_mocked_bindings(
    {
      invisible(capture.output(print(
        mb, print = list(pval = FALSE, I2 = TRUE),
        digits = list(digits = 3, se = 2, tau = 5, tau2 = 6),
        details.methods = TRUE
      )))
      forest(mb, print = list(pval = FALSE, I2 = TRUE),
             digits = list(digits = 3, se = 2, tau = 5, tau2 = 6),
             details = TRUE, filename = file)
    },
    catmeth_metabind = function(x, common, random, prediction,
                                backtransf, print, digits, ...) {
      calls[[length(calls) + 1L]] <<- list(print = print, digits = digits)
      invisible("")
    },
    .package = "meta"
  )

  expect_length(calls, 2)
  for (call in calls) {
    expect_type(call$print, "list")
    expect_type(call$digits, "list")
    expect_true(all(c("func", "method", "method.interval", "tau2", "pval") %in%
                      names(call$print)))
    expect_false(call$print$pval)
    expect_true(call$print$I2)
    expect_equal(call$digits[c("digits", "se", "tau", "tau2")],
                  list(digits = 3, se = 2, tau = 5, tau2 = 6))
  }
})

test_that("method details use list precisions without modifying caller options", {
  ma <- metagen(log(2:4), c(0.1, 0.2, 0.3), sm = "RR",
                 tau.preset = 0.1234567, null.effect = log(1.23456))
  mb <- metabind(ma, ma)
  flags <- meta:::print_defaults(mb, "metabind", func = "print")
  precision <- meta:::digits_defaults("metabind", func = "print")
  precision$digits <- precision$tau <- 3
  original.flags <- flags
  original.precision <- precision

  details <- meta:::catmeth_metabind(
    mb, common = TRUE, random = TRUE, prediction = FALSE,
    backtransf = TRUE, print = flags, digits = precision,
    pscale = 1, irscale = 1, irunit = "units", big.mark = "", forest = TRUE
  )
  expect_match(details, "tau = 0.123", fixed = TRUE)
  expect_match(details, "Null hypothesis: RR = 1.235", fixed = TRUE)
  expect_identical(flags, original.flags)
  expect_identical(precision, original.precision)
})

test_that("bound proportion analyses print method details without Rb columns", {
  mp2e <- metaprop(Ee, Ne, data = Quan2000, studlab = study, common = FALSE)
  mp2c <- metaprop(Ec, Nc, data = Quan2000, studlab = study, common = FALSE)
  mb <- metabind(mp2e, mp2c, name = c("Experimental", "Control"))
  expect_false("upper.Rb" %in% names(mb$overall))

  expect_silent(output <- capture.output(print(
    mb, digits.tau2 = 2, digits.tau = 2, digits.Q = 1, digits.I2 = 0)))
  expect_true(any(grepl("Experimental", output, fixed = TRUE)))
  expect_true(any(grepl("Control", output, fixed = TRUE)))
  expect_true(any(grepl("Details of meta-analysis methods:", output, fixed = TRUE)))
  expect_true(any(grepl("single proportions", output, fixed = TRUE)))
  expect_silent(summary.output <- capture.output(print(summary(mb))))
  expect_true(any(grepl("Details of meta-analysis methods:", summary.output,
                        fixed = TRUE)))

  grDevices::pdf(NULL, width = 16, height = 12)
  on.exit(grDevices::dev.off())
  expect_silent(forest(mb, details = TRUE))
})
