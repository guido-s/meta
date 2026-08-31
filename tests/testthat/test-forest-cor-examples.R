test_that("forest uses symmetric default x-axis for correlations", {
  m <- metacor(c(-0.6, -0.2, 0.3), c(30, 40, 50))

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off())
  f <- forest(m)
  f.analysis <- forest(m, backtransf = FALSE)
  f.xlim <- forest(m, xlim = c(-0.9, 0.6))
  dev.off()

  expect_equal(abs(f$xlim[1]), abs(f$xlim[2]))
  expect_equal(abs(f.analysis$xlim[1]), abs(f.analysis$xlim[2]))
  expect_equal(f.xlim$xlim, c(-0.9, 0.6))
})
