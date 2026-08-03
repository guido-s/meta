test_that("metagen Fleiss1993bin example matches rma.uni", {
  data(Fleiss1993bin)

  ma <- metabin(d.asp, n.asp, d.plac, n.plac, study,
                data = Fleiss1993bin, sm = "RR", method = "I")
  ma.gen <- metagen(TE, seTE, studlab, n.e = n.e + n.c,
                    data = ma, sm = "RR")

  fit.common <- metafor::rma.uni(yi = ma$TE, sei = ma$seTE,
                                 method = "FE", test = "z",
                                 level = ma.gen$level.ma * 100)
  fit.random <- metafor::rma.uni(yi = ma$TE, sei = ma$seTE,
                                 method = ma.gen$method.tau, test = "z",
                                 level = ma.gen$level.ma * 100)

  expect_equal(ma.gen$TE.common, as.numeric(coef(fit.common)))
  expect_equal(ma.gen$seTE.common, fit.common$se)
  expect_equal(ma.gen$lower.common, fit.common$ci.lb)
  expect_equal(ma.gen$upper.common, fit.common$ci.ub)
  expect_equal(ma.gen$pval.common, fit.common$pval)
  expect_equal(ma.gen$statistic.common, fit.common$zval)

  expect_equal(ma.gen$TE.random, as.numeric(coef(fit.random)))
  expect_equal(ma.gen$seTE.random, fit.random$se)
  expect_equal(ma.gen$lower.random, fit.random$ci.lb)
  expect_equal(ma.gen$upper.random, fit.random$ci.ub)
  expect_equal(ma.gen$pval.random, fit.random$pval)
  expect_equal(ma.gen$statistic.random, fit.random$zval)
  expect_true(is.infinite(ma.gen$df.random))
  expect_true(is.na(fit.random$ddf))
  expect_equal(ma.gen$tau2, fit.random$tau2)
})
