metacont_exact_fixture <- function(sm = "SMD", ...) {
  metacont(c(8, 12, 20), c(4, 5, 6), c(1, 1.2, 1.4),
            c(10, 14, 22), c(3, 4, 4.5), c(1.1, 1.3, 1.5),
            sm = sm, method.tau = "DL", ...)
}

test_that("Hedges and Cohen retain exact and approximate formulae", {
  for (method in c("Hedges", "Cohen")) {
    expect_silent(exact <- metacont_exact_fixture(method.smd = method,
                                                  exact.smd = TRUE))
    expect_silent(approx <- metacont_exact_fixture(method.smd = method,
                                                   exact.smd = FALSE))
    expect_identical(exact$exact.smd, TRUE)
    expect_identical(approx$exact.smd, FALSE)
    expect_gt(max(abs(exact$seTE - approx$seTE)), 1e-5)
    if (method == "Cohen")
      expect_equal(exact$TE, approx$TE)
    else
      expect_gt(max(abs(exact$TE - approx$TE)), 1e-5)
    for (value in list(NA, NULL, c(TRUE, FALSE), "invalid"))
      expect_error(metacont_exact_fixture(method.smd = method, exact.smd = value),
                   "Argument 'exact.smd' must be a logical", fixed = TRUE)
  }
})

test_that("Glass accepts NA and warns only for explicit non-NA exact.smd", {
  for (sd in c("control", "experimental")) {
    expect_silent(default <- metacont_exact_fixture(method.smd = "Glass",
                                                    sd.glass = sd))
    expect_silent(explicit <- metacont_exact_fixture(method.smd = "Glass",
                                                     sd.glass = sd, exact.smd = NA))
    expect_identical(default$exact.smd, NA)
    expect_identical(explicit$exact.smd, NA)
    expect_equal(default$TE, explicit$TE)
    expect_equal(default$seTE, explicit$seTE)
    for (value in list(TRUE, FALSE, NULL, NA_real_)) {
      expect_warning(ignored <- metacont_exact_fixture(method.smd = "Glass",
                                                       sd.glass = sd,
                                                       exact.smd = value),
                     "Argument 'exact.smd' is ignored for method.smd = \"Glass\".",
                     fixed = TRUE)
      expect_identical(ignored$exact.smd, NA)
      expect_equal(ignored$TE, default$TE)
      expect_equal(ignored$seTE, default$seTE)
    }
  }
})

test_that("updates to Glass distinguish inherited and explicit exact.smd", {
  fields <- c("TE", "seTE", "TE.common", "seTE.common", "TE.random", "seTE.random")
  for (method in c("Hedges", "Cohen")) {
    for (value in c(TRUE, FALSE)) {
      original <- metacont_exact_fixture(method.smd = method, exact.smd = value)
      expect_silent(glass <- update(original, method.smd = "Glass"))
      expect_identical(glass$method.smd, "Glass")
      expect_identical(glass$exact.smd, NA)
      expect_silent(explicit <- update(original, method.smd = "Glass", exact.smd = NA))
      expect_equal(explicit[fields], glass[fields])
      expect_warning(ignored <- update(original, method.smd = "Glass", exact.smd = value),
                     "Argument 'exact.smd' is ignored for method.smd = \"Glass\".",
                     fixed = TRUE)
      expect_equal(ignored[fields], glass[fields])
      expect_identical(ignored$exact.smd, NA)
      expect_silent(other <- update(glass, method.tau = "PM"))
      expect_identical(other$exact.smd, NA)
      expect_silent(abbreviated <- update(original, method.smd = "G"))
      expect_equal(abbreviated[fields], glass[fields])
    }
  }
})

test_that("updates preserve applicable settings and resolve inherited NA using defaults", {
  old <- settings.meta(exact.smd = TRUE)
  withr::defer(settings.meta(old))
  glass <- metacont_exact_fixture(method.smd = "Glass")
  for (default in c(TRUE, FALSE)) {
    settings.meta(exact.smd = default)
    for (method in c("Hedges", "Cohen")) {
      expect_silent(restored <- update(glass, method.smd = method))
      expect_identical(restored$exact.smd, default)
      for (value in c(TRUE, FALSE)) {
        original <- metacont_exact_fixture(method.smd = method, exact.smd = value)
        expect_silent(unchanged <- update(original, method.tau = "PM"))
        expect_identical(unchanged$exact.smd, value)
        for (sm in c("smd", "S")) {
          expect_silent(normalized <- update(original, sm = sm))
          expect_identical(normalized$sm, "SMD")
          expect_identical(normalized$exact.smd, value)
        }
        expect_silent(explicit <- update(glass, method.smd = method, exact.smd = value))
        expect_identical(explicit$exact.smd, value)
      }
      expect_error(update(glass, method.smd = method, exact.smd = NA),
                   "Argument 'exact.smd' must be a logical", fixed = TRUE)
    }
    md <- metacont_exact_fixture(sm = "MD")
    expect_identical(md$exact.smd, NA)
    expect_silent(smd <- update(md, sm = "smd"))
    expect_identical(smd$exact.smd, default)
    expect_silent(back <- update(smd, sm = "MD"))
    expect_identical(back$exact.smd, NA)
  }
})
