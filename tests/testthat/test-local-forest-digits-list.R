test_that("numeric-list validation honors named length, minimum and missing overrides", {
  values <- list(digits = 2, pval = 1, addcols = c(1, 3), unused = "not numeric")
  expect_silent(meta:::chknumeric_list(values, min = list(pval = 1),
                                      length = list(addcols = 0), miss = list(unused = TRUE)))
  expect_error(meta:::chknumeric_list(values, length = list(addcols = 0)), "values[.]unused")
  values$pval <- 0
  expect_error(meta:::chknumeric_list(values, min = list(pval = 1),
                                      length = list(addcols = 0), miss = list(unused = TRUE)),
               "values[.]pval")
})

test_that("digits lists preserve formatting and vector precision overrides", {
  dat <- data.frame(TE = c(0.123456, 0.234567, 0.345678), seTE = rep(0.1, 3),
                     first = c(1.1234, 2.2345, 3.3456),
                     second = c(4.4567, 5.5678, 6.6789),
                     third = c(7.7891, 8.8912, 9.9123))
  ma <- metagen(TE, seTE, data = dat, method.tau = "DL")
  settings <- list(digits = 3, se = 5, stat = 4, pval = 3, pval.Q = 4,
                    Q = 3, tau2 = 5, tau = 4, I2 = 2, weight = 1, TE = 4,
                    addcols = c(2, 4), addcols.left = c(1, 3), addcols.right = 2)
  observed <- NULL
  draw <- meta:::forest_meta_internal
  testthat::local_mocked_bindings(
    forest_meta_internal = function(...) {
      observed <<- list(...)
      do.call(draw, observed)
    }, .package = "meta")
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())

  values <- function(column) {
    labels <- lapply(observed$cols[[paste0("col.", column)]]$labels, `[[`, "label")
    unname(unlist(Filter(function(z)
      is.character(z) && length(z) == 1 && grepl("^-?[0-9]", trimws(z)), labels)))
  }
  baseline <- NULL
  for (legacy in c(FALSE, TRUE)) {
    if (legacy) {
      options <- settings
      names(options) <- c("digits", paste0("digits.", names(settings)[-1]))
    }
    else
      options <- list(digits = settings)
    result <- do.call(meta::forest,
                      c(list(x = ma,
                             leftcols = c("studlab", "first", "second", "effect", "ci"),
                             rightcols = "third"), options))
    expect_equal(trimws(values("first")), meta:::formatN(dat$first, digits = 1))
    expect_equal(trimws(values("second")), meta:::formatN(dat$second, digits = 3))
    expect_equal(trimws(values("third")), meta:::formatN(dat$third, digits = 2))
    expect_equal(trimws(result$TE.format[nzchar(result$TE.format)]),
                 meta:::formatN(ma$TE, digits = 4))
    if (!legacy)
      baseline <- result
    else {
      expect_equal(result$effect.format, baseline$effect.format)
      expect_equal(result$ci.format, baseline$ci.format)
      expect_equal(result$seTE.format, baseline$seTE.format)
    }
  }
  expect_warning(forest(ma, digits.zval = 5,
                         leftcols = c("studlab", "statistic", "effect", "ci")), "digits.stat")
  expect_equal(trimws(values("statistic")),
               meta:::formatN(c(ma$statistic.common, ma$statistic.random, ma$statistic),
                               digits = 5))
  expect_error(forest(ma, digits = list(pval = 0)), "digits[.]pval")
  expect_error(forest(ma, digits = list(pval.Q = 0)), "digits[.]pval.Q")
  expect_error(forest(ma, digits = list(mean = 1:2)), "digits[.]mean")
})
