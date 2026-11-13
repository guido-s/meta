test_that("funnel accepts user-specified x-axis tick marks", {
  m <- metaprop(1:5, 101:105)

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off())
  expect_silent(funnel(m, at = c(0.02, 0.04, 0.06)))
  expect_silent(funnel(m, backtransf = FALSE, at = c(-4, -3, -2)))
  dev.off()
})

test_that("bubble plot uses stored null effect for transformed correlations", {
  dat <- data.frame(cor = c(0.1, 0.2, 0.3, 0.4, 0.5),
                    n = c(30, 40, 35, 50, 45),
                    x = 1:5)
  m <- metacor(cor, n, data = dat, sm = "ZCOR", null.effect = 0.2)
  mr <- metareg(m, ~ x)

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off())
  expect_invisible(bubble(mr, backtransf = TRUE))
  expect_invisible(bubble(mr, backtransf = FALSE))
  dev.off()
})
