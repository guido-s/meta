test_that("funnel back transforms supported proportion transformations", {
  dat <- data.frame(event = c(2, 5, 8, 12, 18),
                    n = c(30, 40, 35, 50, 45))
  
  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off())
  for (sm in c("PLOGIT", "PLN", "PAS", "PFT")) {
    m <- metaprop(event, n, data = dat, sm = sm, null.effect = 0.2)
    f <- funnel(m, backtransf = TRUE, xlim = c(0.01, 0.5),
                ref.triangle = TRUE)
    fc <- funnel(m, type = "contour", backtransf = TRUE, xlim = c(0.01, 0.5))
    fa <- funnel(m, backtransf = FALSE)
    xlim <- if (sm == "PFT")
      transf(c(0.01, 0.5), sm, n = m$n.harmonic.mean)
    else
      transf(c(0.01, 0.5), sm)
    
    expect_equal(f$xvals, m$TE)
    expect_equal(f$xlim, xlim)
    expect_equal(fa$xvals, m$TE)
    expect_named(fc, c("xvals", "yvals", "pch", "text", "cex", "col", "bg",
                      "cex.studlab", "xlim", "ylim", "contour.levels",
                      "col.contour", "text.contour"))
  }
  dev.off()
})

test_that("funnel back transforms supported rate transformations", {
  dat <- data.frame(event = c(2, 5, 8, 12, 18),
                    time = c(30, 40, 35, 50, 45))
  
  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off())
  for (sm in c("IRLN", "IRS", "IRFT")) {
    m <- metarate(event, time, data = dat, sm = sm, null.effect = 0.2)
    f <- funnel(m, backtransf = TRUE, xlim = c(0.01, 0.5),
                ref.triangle = TRUE)
    fc <- funnel(m, type = "contour", backtransf = TRUE, xlim = c(0.01, 0.5))
    fa <- funnel(m, backtransf = FALSE)
    xlim <- if (sm == "IRFT")
      transf(c(0.01, 0.5), sm, time = m$t.harmonic.mean)
    else
      transf(c(0.01, 0.5), sm)
    
    expect_equal(f$xvals, m$TE)
    expect_equal(f$xlim, xlim)
    expect_equal(fa$xvals, m$TE)
    expect_named(fc, c("xvals", "yvals", "pch", "text", "cex", "col", "bg",
                      "cex.studlab", "xlim", "ylim", "contour.levels",
                      "col.contour", "text.contour"))
  }
  dev.off()
})

test_that("funnel back transforms mean and correlation transformations", {
  dat.mean <- data.frame(n = c(20, 30, 40, 25, 35),
                         mean = c(2, 3, 4, 5, 6),
                         sd = c(0.5, 0.7, 0.8, 1.0, 1.1))
  dat.cor <- data.frame(cor = c(0.1, 0.2, 0.3, 0.4, 0.5),
                        n = c(30, 40, 35, 50, 45))
  
  m.mean <- metamean(n, mean, sd, data = dat.mean, sm = "MLN",
                     null.effect = 3)
  m.cor <- metacor(cor, n, data = dat.cor, sm = "ZCOR", null.effect = 0.2)
  
  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off())
  f.mean <- funnel(m.mean, backtransf = TRUE, xlim = c(1, 7),
                   ref.triangle = TRUE)
  f.cor <- funnel(m.cor, backtransf = TRUE, xlim = c(0, 0.6),
                  ref.triangle = TRUE)
  fa.mean <- funnel(m.mean, backtransf = FALSE)
  fa.cor <- funnel(m.cor, backtransf = FALSE)
  dev.off()
  
  expect_equal(f.mean$xvals, m.mean$TE)
  expect_equal(f.mean$xlim, transf(c(1, 7), "MLN"))
  expect_equal(f.cor$xvals, m.cor$TE)
  expect_equal(f.cor$xlim, transf(c(0, 0.6), "ZCOR"))
  expect_equal(fa.mean$xvals, m.mean$TE)
  expect_equal(fa.cor$xvals, m.cor$TE)
})

test_that("funnel still back transforms relative effects on log axis", {
  dat <- data.frame(event.e = c(2, 5, 8, 12, 18),
                    n.e = c(30, 40, 35, 50, 45),
                    event.c = c(4, 6, 9, 10, 16),
                    n.c = c(32, 42, 38, 48, 44))
  m <- metabin(event.e, n.e, event.c, n.c, data = dat, sm = "RR")

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off())
  f <- funnel(m, backtransf = TRUE)
  fa <- funnel(m, backtransf = FALSE)
  dev.off()

  expect_equal(f$xvals, exp(m$TE))
  expect_equal(fa$xvals, m$TE)
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
