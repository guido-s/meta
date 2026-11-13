test_that("printed VE BLUP limits are ordered for defaults and explicit overrides", {
  for (subgroup in c(FALSE, TRUE)) for (default in c(FALSE, TRUE)) {
    b <- blup(blup_bounds_fixture(subgroup = subgroup), backtransf = default)
    before <- b
    default.output <- capture.output(print(b, digits = 6, details = FALSE))
    for (backtransf in c(FALSE, TRUE)) {
      expect_silent(output <- capture.output(print(
        b, backtransf = backtransf, digits = 6, details = FALSE)))
      lower <- if (backtransf) 100 * (1 - exp(b$upper)) else b$lower
      upper <- if (backtransf) 100 * (1 - exp(b$lower)) else b$upper
      expected <- meta:::formatCI(
        meta:::formatN(round(lower, 6), digits = 6, lab.NA = "NA"),
        meta:::formatN(round(upper, 6), digits = 6, lab.NA = "NA"))
      compact <- gsub("[[:space:]]", "", output)
      for (interval in unique(gsub("[[:space:]]", "", expected)))
        expect_true(any(grepl(interval, compact, fixed = TRUE)))
      if (backtransf == default)
        expect_identical(output, default.output)
    }
    expect_identical(b, before)
  }
})
