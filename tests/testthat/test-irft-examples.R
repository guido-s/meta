test_that("transf handles Freeman-Tukey transformed incidence rates", {
  rate <- c(0.1, 0.2, 0.3)
  time <- c(10, 20, 30)
  expected <- 0.5 * (sqrt(rate) + sqrt(rate + 1 / time))
  
  expect_equal(ir2asin(rate, time), expected)
  expect_equal(transf(rate, "IRFT", time = time), expected)
  expect_equal(transf(rate, "IRFT"), sqrt(rate))
  expect_equal(asin2ir(ir2asin(rate, time), time), rate)
  expect_equal(asin2ir(sqrt(rate)), rate)
  expect_equal(backtransf(sqrt(rate), "IRFT"), rate)
  expect_equal(asin2ir(ir2asin(rate, 20), 20), rate)
})

test_that("metarate uses Freeman-Tukey transformation and harmonic mean", {
  event <- c(2, 5, 8, 12)
  time <- c(30, 40, 35, 50)
  m <- metarate(event, time, sm = "IRFT", null.effect = 0.2)
  
  expect_equal(m$TE, ir2asin(event / time, time))
  expect_equal(m$t.harmonic.mean, 1 / mean(1 / time))
  expect_equal(m$statistic.random,
               (m$TE.random - ir2asin(0.2, m$t.harmonic.mean)) /
                 m$seTE.random)
})

test_that("cidprop transforms IRFT decision thresholds using harmonic mean", {
  m <- metarate(c(2, 5, 8, 12), c(30, 40, 35, 50),
                sm = "IRFT", null.effect = 0.2, prediction = TRUE)
  pp <- cidprop(m, cid.below.null = 0.1, cid.above.null = 0.4)
  
  cid.below <- ir2asin(0.1, m$t.harmonic.mean)
  cid.above <- ir2asin(0.4, m$t.harmonic.mean)
  
  expect_equal(pp$prop.cid.below.null,
               pt((cid.below - m$TE.random) / m$seTE.predict, m$df.predict))
  expect_equal(pp$prop.cid.above.null,
               pt((cid.above - m$TE.random) / m$seTE.predict, m$df.predict,
                  lower.tail = FALSE))
})

test_that("trimfill IRFT cidprop uses square root threshold transformation", {
  m <- metarate(c(2, 5, 8, 12, 18), c(30, 40, 35, 50, 45),
                sm = "IRFT", null.effect = 0.2, prediction = TRUE)
  tf <- trimfill(m)
  pp <- cidprop(tf, cid.below.null = 0.1, cid.above.null = 0.4)
  
  expect_null(tf$t.harmonic.mean)
  expect_equal(pp$prop.cid.below.null,
               pt((sqrt(0.1) - tf$TE.random) / tf$seTE.predict,
                  tf$df.predict))
  expect_equal(pp$prop.cid.above.null,
               pt((sqrt(0.4) - tf$TE.random) / tf$seTE.predict,
                  tf$df.predict, lower.tail = FALSE))
})

test_that("bubble and funnel plots back transform IRFT results", {
  dat <- data.frame(event = c(2, 5, 8, 12, 18),
                    time = c(30, 40, 35, 50, 45),
                    x = 1:5,
                    group = c("a", "a", "b", "b", "b"))
  m <- metarate(event, time, data = dat, sm = "IRFT", null.effect = 0.2)
  mr1 <- metareg(m, ~ x)
  mr2 <- metareg(m, ~ group)
  
  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off())
  expect_invisible(bubble(mr1, backtransf = TRUE, irscale = 100))
  expect_invisible(bubble(mr1, backtransf = FALSE))
  expect_invisible(bubble(mr2, backtransf = TRUE))
  f <- funnel(m, backtransf = TRUE, xlim = c(0, 0.5), ref.triangle = TRUE)
  fc <- funnel(m, type = "contour", backtransf = TRUE, xlim = c(0, 0.5))
  fa <- funnel(m, backtransf = FALSE)
  dev.off()
  
  expect_equal(f$xvals, m$TE)
  expect_equal(f$xlim, transf(c(0, 0.5), "IRFT", time = m$t.harmonic.mean))
  expect_equal(fa$xvals, m$TE)
  expect_named(fc, c("xvals", "yvals", "pch", "text", "cex", "col", "bg",
                    "cex.studlab", "xlim", "ylim", "contour.levels",
                    "col.contour", "text.contour"))
})

test_that("IRFT restrictions avoid ambiguous Freeman-Tukey times", {
  m <- metarate(c(2, 5, 8, 12), c(30, 40, 35, 50), sm = "IRFT")
  
  expect_error(metaadd(m, type = "random", TE = 0.2, lower = 0.1, upper = 0.3),
               "not available")
  expect_error(metagen(c(0.1, 0.2), c(0.01, 0.02), sm = "IRFT",
                       transf = FALSE),
               "not available")
})
