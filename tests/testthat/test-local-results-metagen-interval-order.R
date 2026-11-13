test_that("metagen orders monotone transformed CIs before deriving SEs", {
  cases <- list(
    list(sm = "VE", lower = c(20, 40), upper = c(60, 80),
         forward = function(x) log(1 - x / 100),
         inverse = function(x) 100 * (1 - exp(x)), args = list(),
         decreasing = TRUE),
    list(sm = "VE", lower = c(1, 4), upper = c(9, 16),
         forward = log, inverse = exp,
         args = list(func.transf = "log", args.transf = list(base = exp(1)),
                     func.backtransf = "exp"), decreasing = FALSE),
    list(sm = "MD", lower = c(0.1, 0.2), upper = c(0.8, 0.95),
         forward = function(x) qnorm(x, lower.tail = FALSE),
         inverse = function(x) pnorm(x, lower.tail = FALSE),
         args = list(func.transf = "qnorm",
                     args.transf = list(lower.tail = FALSE),
                     func.backtransf = "pnorm",
                     args.backtransf = list(lower.tail = FALSE)),
         decreasing = TRUE)
  )

  for (case in cases) {
    for (n in 1:2) {
      lower <- case$lower[seq_len(n)]
      upper <- case$upper[seq_len(n)]
      lo <- if (case$decreasing) case$forward(upper) else case$forward(lower)
      hi <- if (case$decreasing) case$forward(lower) else case$forward(upper)
      te <- (lo + hi) / 2
      dat <- data.frame(TE = case$inverse(te), lower = lower, upper = upper)
      original <- dat
      for (supply.TE in c(TRUE, FALSE)) {
        args <- c(list(lower = lower, upper = upper, sm = case$sm,
                       transf = FALSE, level = 0.95, level.ci = 0.95,
                       method.tau = "DL", keepdata = TRUE), case$args)
        if (supply.TE) {
          args$TE <- quote(TE)
          args$data <- dat
        }
        m <- do.call(metagen, args)

        expect_equal(m$TE, te)
        expect_equal(m$lower, lo)
        expect_equal(m$upper, hi)
        expect_equal(m$seTE, (hi - lo) / (2 * qnorm(0.975)))
        expect_true(all(m$seTE > 0))
        expect_equal(m$data$.lower, lo)
        expect_equal(m$data$.upper, hi)
        expect_equal(m$data$.lower.orig, if (case$decreasing) upper else lower)
        expect_equal(m$data$.upper.orig, if (case$decreasing) lower else upper)
        if (supply.TE)
          expect_equal(m$data$.TE.orig, dat$TE)
      }
      expect_identical(dat, original)
    }
  }
})

test_that("metagen swaps only comparable reversed endpoints including infinities", {
  cases <- list(
    list(sm = "VE", forward = function(x) log(1 - x / 100),
         args = list(), te = 60),
    list(sm = "VE", forward = identity,
         args = list(func.transf = "identity", func.backtransf = "identity"),
         te = 50),
    list(sm = "MD", forward = function(x) -x,
         args = list(func.transf = "-", func.backtransf = "-"), te = 50)
  )

  for (case in cases) {
    dat <- data.frame(TE = case$te, seTE = 0.2,
                      lower = c(20, NA, 20, NA, -Inf, 20, -Inf, 50),
                      upper = c(80, 80, NA, NA, 80, 100, 100, 50))
    finite <- is.finite(dat$lower) & is.finite(dat$upper)
    if (length(case$args))
      dat$TE[finite] <- (dat$lower[finite] + dat$upper[finite]) / 2
    original <- dat
    # Include scalar missing endpoints as well as mixed vector masks.
    for (rows in c(list(seq_len(nrow(dat))), as.list(2:8))) {
      d <- dat[rows, , drop = FALSE]
      lo <- case$forward(d$lower)
      hi <- case$forward(d$upper)
      swap <- !is.na(lo) & !is.na(hi) & lo > hi
      raw.lo <- d$lower
      raw.hi <- d$upper
      lo[swap] <- case$forward(d$upper[swap])
      hi[swap] <- case$forward(d$lower[swap])
      raw.lo[swap] <- d$upper[swap]
      raw.hi[swap] <- d$lower[swap]
      m <- do.call(metagen,
                   c(list(TE = quote(TE), seTE = quote(seTE),
                          lower = quote(lower), upper = quote(upper), data = d,
                          sm = case$sm, transf = FALSE, method.tau = "DL",
                          level = 0.95, level.ci = 0.95, keepdata = TRUE),
                     case$args))

      expect_equal(m$TE, case$forward(d$TE))
      expect_equal(m$seTE, d$seTE)
      expect_equal(m$data$.lower, lo)
      expect_equal(m$data$.upper, hi)
      expect_equal(m$data$.lower.orig, raw.lo)
      expect_equal(m$data$.upper.orig, raw.hi)
      expect_equal(case$forward(m$data$.lower.orig), lo)
      expect_equal(case$forward(m$data$.upper.orig), hi)
      expected.lo <- case$forward(d$TE) - qnorm(0.975) * d$seTE
      expected.hi <- case$forward(d$TE) + qnorm(0.975) * d$seTE
      expected.lo[!is.na(lo)] <- lo[!is.na(lo)]
      expected.hi[!is.na(hi)] <- hi[!is.na(hi)]
      expect_equal(m$lower, expected.lo)
      expect_equal(m$upper, expected.hi)
    }
    expect_identical(dat, original)
  }
})
