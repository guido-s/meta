test_that("missing custom functions use the summary-measure transformation", {
  for (fun in list(NULL, NA, NA_character_, "")) {
    expect_equal(backtransf(0.3, "MD", func = fun), 0.3)
    expect_equal(backtransf(0.3, "RR", func = fun), exp(0.3))
  }
  expect_equal(backtransf(0.3, "MD", func = "exp"), exp(0.3))
  expect_equal(backtransf(0.3, "MD", func = exp), exp(0.3))
})

test_that("mixed custom and default metabind transformations can be printed", {
  ma <- metagen(c(0.2, 0.3, 0.4), rep(0.1, 3), sm = "MD",
                 func.backtransf = "exp")
  other <- metagen(c(0.2, 0.3, 0.4), rep(0.1, 3), sm = "MD")
  expect_output(print(metabind(ma, other)), "Common effect model")
})
