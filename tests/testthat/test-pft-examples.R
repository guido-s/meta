test_that("transf handles Freeman-Tukey transformed proportions", {
  p <- c(0.1, 0.2, 0.3)
  n <- c(10, 20, 30)
  expected <-
    0.5 * (asin(sqrt(p * n / (n + 1))) +
             asin(sqrt((p * n + 1) / (n + 1))))
  
  expect_equal(transf(p, "PFT", n = n), expected)
  expect_equal(transf(p, "PFT"), p2asin(p))
  expect_equal(asin2p(p2asin(p, 20), 20), p)
})

test_that("metaprop uses Freeman-Tukey transformation and harmonic mean", {
  event <- c(1, 5, 10, 20)
  n <- c(10, 30, 50, 100)
  m <- metaprop(event, n, sm = "PFT")
  
  expect_equal(m$TE, p2asin(event / n, n))
  expect_equal(m$n.harmonic.mean, 1 / mean(1 / n))
})

test_that("cidprop transforms PFT decision thresholds using harmonic mean", {
  m <- metaprop(c(1, 5, 10, 20), c(10, 30, 50, 100),
                sm = "PFT", prediction = TRUE)
  pp <- cidprop(m, cid.below.null = 0.1, cid.above.null = 0.4)
  
  cid.below <- p2asin(0.1, m$n.harmonic.mean)
  cid.above <- p2asin(0.4, m$n.harmonic.mean)
  
  expect_equal(pp$prop.cid.below.null,
               pt((cid.below - m$TE.random) / m$seTE.predict, m$df.predict))
  expect_equal(pp$prop.cid.above.null,
               pt((cid.above - m$TE.random) / m$seTE.predict, m$df.predict,
                  lower.tail = FALSE))
})

test_that("trimfill PFT cidprop uses arcsine threshold transformation", {
  m <- metaprop(c(1, 5, 10, 20, 25), c(10, 30, 50, 100, 120),
                sm = "PFT", prediction = TRUE)
  tf <- trimfill(m)
  pp <- cidprop(tf, cid.below.null = 0.1, cid.above.null = 0.4)
  
  expect_null(tf$n.harmonic.mean)
  expect_equal(pp$prop.cid.below.null,
               pt((p2asin(0.1) - tf$TE.random) / tf$seTE.predict,
                  tf$df.predict))
  expect_equal(pp$prop.cid.above.null,
               pt((p2asin(0.4) - tf$TE.random) / tf$seTE.predict,
                  tf$df.predict, lower.tail = FALSE))
})

test_that("bubble plot back transforms PFT results", {
  dat <- data.frame(event = c(2, 5, 8, 12, 18),
                    n = c(30, 40, 35, 50, 45),
                    x = 1:5,
                    group = c("a", "a", "b", "b", "b"))
  m <- metaprop(event, n, data = dat, sm = "PFT", null.effect = 0.2)
  mr1 <- metareg(m, ~ x)
  mr2 <- metareg(m, ~ group)
  
  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off())
  expect_invisible(bubble(mr1, backtransf = TRUE, pscale = 100))
  expect_invisible(bubble(mr1, backtransf = FALSE))
  expect_invisible(bubble(mr2, backtransf = TRUE))
  dev.off()
})

test_that("funnel plot shows PFT tick marks on original scale", {
  dat <- data.frame(event = c(2, 5, 8, 12, 18),
                    n = c(30, 40, 35, 50, 45))
  m <- metaprop(event, n, data = dat, sm = "PFT", null.effect = 0.2)
  
  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off())
  f <- funnel(m, backtransf = TRUE, ref.triangle = TRUE)
  fc <- funnel(m, type = "contour", backtransf = TRUE)
  fa <- funnel(m, backtransf = FALSE)
  dev.off()
  
  expect_equal(f$xvals, m$TE)
  expect_equal(fa$xvals, m$TE)
  expect_named(fc, c("xvals", "yvals", "pch", "text", "cex", "col", "bg",
                    "cex.studlab", "xlim", "ylim", "contour.levels",
                    "col.contour", "text.contour"))
})
