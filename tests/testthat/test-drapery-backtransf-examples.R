test_that("drapery shows back-transformed tick marks for single proportions", {
  m <- metaprop(1:5, 101:105)

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off())
  expect_silent(res <- drapery(m, xlim = c(0.001, 0.1)))
  expect_silent(drapery(m))
  expect_silent(drapery(m, axes = FALSE))
  expect_silent(drapery(m, at = c(0.02, 0.04, 0.06)))
  dev.off()

  expect_equal(attr(res, "xlim"), transf(c(0.001, 0.1), m$sm))
})

test_that("drapery labels x-axis for single proportions", {
  m <- metaprop(1:5, 101:105)
  assign(".drapery_xlab", NULL, envir = .GlobalEnv)
  on.exit(rm(".drapery_xlab", envir = .GlobalEnv), add = TRUE)

  suppressMessages(capture.output(
    trace(graphics::plot, quote({
      .drapery_xlab <<- list(...)$xlab
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(graphics::plot))),
          add = TRUE)
  #
  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  expect_silent(drapery(m))
  dev.off()

  expect_equal(get(".drapery_xlab", envir = .GlobalEnv), "Proportion")
})

test_that("drapery transforms original-scale x-limits for axis-only measures", {
  dat.rate <- data.frame(event = c(2, 5, 8, 12, 18),
                         time = c(30, 40, 35, 50, 45))
  dat.cor <- data.frame(cor = c(0.1, 0.2, 0.3, 0.4, 0.5),
                        n = c(30, 40, 35, 50, 45))
  m.irft <- metarate(event, time, data = dat.rate, sm = "IRFT")
  m.zcor <- metacor(cor, n, data = dat.cor, sm = "ZCOR")

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off())
  res.irft <- drapery(m.irft, xlim = c(0.01, 0.5))
  res.zcor <- drapery(m.zcor, xlim = c(0, 0.6))
  dev.off()

  expect_equal(attr(res.irft, "xlim"),
               transf(c(0.01, 0.5), "IRFT", time = m.irft$t.harmonic.mean))
  expect_equal(attr(res.zcor, "xlim"), transf(c(0, 0.6), "ZCOR"))
})
