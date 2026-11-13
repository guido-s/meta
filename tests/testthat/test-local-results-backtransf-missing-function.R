test_that("missing custom functions use the summary-measure transformation", {
  for (fun in list(NULL, NA, NA_character_, "")) {
    expect_equal(backtransf(0.3, "MD", func = fun), 0.3)
    expect_equal(backtransf(0.3, "RR", func = fun), exp(0.3))
  }
  expect_equal(backtransf(0.3, "MD", func = "exp"), exp(0.3))
  expect_equal(backtransf(0.3, "MD", func = exp), exp(0.3))
})
