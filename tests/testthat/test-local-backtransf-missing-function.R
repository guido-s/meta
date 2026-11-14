test_that("mixed custom and default metabind transformations can be printed", {
  ma <- metagen(c(0.2, 0.3, 0.4), rep(0.1, 3), sm = "MD",
                 func.backtransf = "exp")
  other <- metagen(c(0.2, 0.3, 0.4), rep(0.1, 3), sm = "MD")
  expect_output(print(metabind(ma, other)), "Common effect model")
})
