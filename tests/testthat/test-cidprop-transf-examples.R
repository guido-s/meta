test_that("cidprop transforms thresholds for supported summary measures", {
  dat.prop <- data.frame(event = c(2, 5, 8, 12, 18),
                         n = c(30, 40, 35, 50, 45))
  dat.rate <- data.frame(event = c(2, 5, 8, 12, 18),
                         time = c(30, 40, 35, 50, 45))
  dat.mean <- data.frame(n = c(20, 30, 40, 25, 35),
                         mean = c(2, 3, 4, 5, 6),
                         sd = c(0.5, 0.7, 0.8, 1.0, 1.1))
  dat.cor <- data.frame(cor = c(0.1, 0.2, 0.3, 0.4, 0.5),
                        n = c(30, 40, 35, 50, 45))
  
  tests <- list(
    list(m = metaprop(event, n, data = dat.prop, sm = "PLOGIT",
                      prediction = TRUE), lower = 0.1, upper = 0.4),
    list(m = metaprop(event, n, data = dat.prop, sm = "PLN",
                      prediction = TRUE), lower = 0.1, upper = 0.4),
    list(m = metaprop(event, n, data = dat.prop, sm = "PAS",
                      prediction = TRUE), lower = 0.1, upper = 0.4),
    list(m = metarate(event, time, data = dat.rate, sm = "IRLN",
                      prediction = TRUE), lower = 0.1, upper = 0.4),
    list(m = metarate(event, time, data = dat.rate, sm = "IRS",
                      prediction = TRUE), lower = 0.1, upper = 0.4),
    list(m = metarate(event, time, data = dat.rate, sm = "IRFT",
                      prediction = TRUE), lower = 0.1, upper = 0.4),
    list(m = metamean(n, mean, sd, data = dat.mean, sm = "MLN",
                      prediction = TRUE), lower = 2, upper = 5),
    list(m = metacor(cor, n, data = dat.cor, sm = "ZCOR",
                     prediction = TRUE), lower = 0.1, upper = 0.4)
  )
  
  for (x in tests) {
    pp <- cidprop(x$m, cid.below.null = x$lower, cid.above.null = x$upper)
    cid.below <- transf(x$lower, x$m$sm, time = x$m$t.harmonic.mean)
    cid.above <- transf(x$upper, x$m$sm, time = x$m$t.harmonic.mean)
    
    expect_equal(pp$prop.cid.below.null,
                 pt((cid.below - x$m$TE.random) / x$m$seTE.predict,
                    x$m$df.predict))
    expect_equal(pp$prop.cid.above.null,
                 pt((cid.above - x$m$TE.random) / x$m$seTE.predict,
                    x$m$df.predict, lower.tail = FALSE))
  }
})
