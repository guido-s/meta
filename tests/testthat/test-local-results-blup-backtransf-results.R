test_that("BLUP storage stays on the ordered analysis scale for either display default", {
  for (subgroup in c(FALSE, TRUE)) {
    ma <- blup_bounds_fixture(subgroup = subgroup)
    analysis <- blup(ma, backtransf = FALSE)
    natural <- blup(ma, backtransf = TRUE)
    for (field in c("blup", "se.blup", "lower", "upper"))
      expect_identical(analysis[[field]], natural[[field]])
    expect_true(all(analysis$lower <= analysis$upper))
    margin <- qnorm((1 + attr(analysis, "x")$level) / 2) * analysis$se.blup
    expect_equal(analysis$lower, analysis$blup - margin)
    expect_equal(analysis$upper, analysis$blup + margin)
  }
})

test_that("extracted VE BLUP limits reverse only after back-transformation", {
  for (subgroup in c(FALSE, TRUE)) for (default in c(FALSE, TRUE)) {
    b <- blup(blup_bounds_fixture(subgroup = subgroup), backtransf = default)
    before <- b
    extracted <- estimates(b, writexl = FALSE)
    expect_identical(attr(extracted, "x")$backtransf, default)
    for (backtransf in c(FALSE, TRUE)) {
      res <- estimates(b, backtransf = backtransf, writexl = FALSE)
      expect_length(res$lower, nrow(b))
      expect_length(res$upper, nrow(b))
      expect_true(all(res$lower <= res$upper))
      if (backtransf) {
        expect_equal(res$blup, 100 * (1 - exp(b$blup)))
        expect_equal(res$lower, 100 * (1 - exp(b$upper)))
        expect_equal(res$upper, 100 * (1 - exp(b$lower)))
      }
      else {
        expect_identical(res$blup, b$blup)
        expect_identical(res$lower, b$lower)
        expect_identical(res$upper, b$upper)
      }
      if (backtransf == default) {
        expect_identical(res$lower, extracted$lower)
        expect_identical(res$upper, extracted$upper)
      }
    }
    expect_identical(b, before)
  }
})

test_that("increasing BLUP back-transformations preserve limit order and vector length", {
  for (sm in c("MD", "RR")) {
    b <- blup(blup_bounds_fixture(sm = sm))
    before <- b
    res <- estimates(b, backtransf = TRUE, writexl = FALSE)
    transform <- if (sm == "RR") exp else identity
    expect_length(res$lower, nrow(b))
    expect_equal(res$lower, transform(b$lower))
    expect_equal(res$upper, transform(b$upper))
    expect_true(all(res$lower <= res$upper))
    expect_identical(b, before)
  }
})

test_that("BLUP bound normalization handles infinity and missing endpoints", {
  b <- blup(blup_bounds_fixture())
  b$lower[1] <- -Inf
  b$upper[1] <- Inf
  b$lower[2] <- NA_real_
  b$upper[3] <- NA_real_
  before <- b
  res <- estimates(b, backtransf = TRUE, writexl = FALSE)
  expect_equal(res$lower[1], -Inf)
  expect_equal(res$upper[1], 100)
  expect_true(is.na(res$lower[2]))
  expect_true(is.na(res$upper[3]))
  expect_length(res$lower, nrow(b))
  complete <- !is.na(b$lower) & !is.na(b$upper)
  expect_true(all(res$lower[complete] <= res$upper[complete]))
  expect_identical(b, before)
})

test_that("meta extraction retains full CI vectors with and without reversed transformations", {
  for (sm in c("MD", "RR", "VE")) {
    ma <- blup_bounds_fixture(sm = sm)
    before <- ma
    raw <- estimates(ma, study.results = TRUE, common = TRUE, random = TRUE,
                     prediction = TRUE, backtransf = FALSE, digits = 12,
                     writexl = FALSE)
    res <- estimates(ma, study.results = TRUE, common = TRUE, random = TRUE,
                     prediction = TRUE, backtransf = TRUE, digits = 12,
                     writexl = FALSE)
    transform <- switch(sm, MD = identity, RR = exp, VE = function(x) 100 * (1 - exp(x)))
    lower <- transform(raw$lower)
    upper <- transform(raw$upper)
    expect_equal(res$lower, pmin(lower, upper), tolerance = 1e-9)
    expect_equal(res$upper, pmax(lower, upper), tolerance = 1e-9)
    expect_equal(nrow(res), nrow(raw))
    expect_true(all(res$lower <= res$upper))
    expect_identical(ma, before)
    without.ci <- estimates(ma, backtransf = TRUE, ci = FALSE, writexl = FALSE)
    expect_false(any(c("lower", "upper") %in% names(without.ci)))
  }
})
