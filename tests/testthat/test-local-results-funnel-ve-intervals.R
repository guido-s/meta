test_that("VE funnel envelopes and contours stay ordered in reflected coordinates", {
  m <- metagen(c(-0.9, -0.5, -0.2, 0.1), c(0.12, 0.2, 0.25, 0.3),
               sm = "VE", method.tau = "DL")
  original <- m
  real.ci <- meta::ci
  observed <- list()
  calls <- 0L
  testthat::local_mocked_bindings(
    ci = function(...) {
      calls <<- calls + 1L
      frame <- parent.frame()
      # The first contour call sees reflected envelopes; subsequent calls
      # see the preceding contour after reflection and polygon construction.
      if (calls >= 3L) {
        observed[[length(observed) + 1L]] <<- list(
          TE = frame$TE, common = frame$TE.common, random = frame$TE.random,
          ref = frame$ref, pooled = frame$ciTE, reference = frame$ciTE.ref,
          contour = if (calls > 3L) frame$ciContour else NULL,
          contour.level = if (calls > 3L) frame$contour.levels[calls - 3L]
                          else NULL)
      }
      real.ci(...)
    }, .package = "meta")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  for (backtransf in c(TRUE, FALSE)) {
    for (ref in c(0, 0.35)) {
      for (random.only in c(FALSE, TRUE)) {
        calls <- 0L
        observed <- list()
        # A final extra contour call lets the mock observe both 90% and 95%.
        f <- funnel(m, backtransf = backtransf, ref = ref,
                    common = !random.only, random = random.only,
                    ref.triangle = TRUE, level = 0.95,
                    contour.levels = c(0.90, 0.95, 0.99),
                    col.contour = c("gray80", "gray70", "gray60"))
        sign <- if (backtransf) -1 else 1
        centre <- sign * if (random.only) m$TE.random else m$TE.common
        se <- seq(0, max(m$seTE), length.out = 500)
        expect_identical(calls, 5L)
        expect_length(observed, 3L)
        expect_equal(f$xvals, sign * m$TE)
        expect_equal(f$yvals, m$seTE)
        for (state in observed) {
          expect_equal(state$TE, sign * m$TE)
          expect_equal(state$common, sign * m$TE.common)
          expect_equal(state$random, sign * m$TE.random)
          expect_equal(state$ref, sign * ref)
          expect_equal(state$pooled$lower, centre - qnorm(0.975) * se)
          expect_equal(state$pooled$upper, centre + qnorm(0.975) * se)
          expect_equal(state$reference$lower, sign * ref - qnorm(0.975) * se)
          expect_equal(state$reference$upper, sign * ref + qnorm(0.975) * se)
          expect_true(all(state$pooled$lower <= state$pooled$upper))
          expect_true(all(state$reference$lower <= state$reference$upper))
          if (!is.null(state$contour)) {
            z <- qnorm((1 + state$contour.level) / 2)
            expect_equal(state$contour$TE, sign * ref)
            expect_equal(state$contour$lower, sign * ref - z * rev(se))
            expect_equal(state$contour$upper, sign * ref + z * rev(se))
            expect_equal((state$contour$lower + state$contour$upper) / 2,
                         rep(sign * ref, length(se)))
            expect_true(all(state$contour$lower <= state$contour$upper))
          }
        }
        expect_identical(m, original)
      }
    }
  }
})
