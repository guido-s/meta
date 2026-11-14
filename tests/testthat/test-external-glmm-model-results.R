test_that("six-study binomial GLMM matches direct ML and FE fits", {
  dat <- data.frame(event = c(8, 18, 31, 44, 57, 69),
                    total = c(80, 95, 100, 110, 105, 120))
  ma <- metaprop(event, total, data = dat, sm = "PLOGIT", method = "GLMM",
                 method.tau = "ML", method.random.ci = "classic",
                 method.I2 = "tau2", common = TRUE, null.effect = 0.5)
  random <- metafor::rma.glmm(measure = "PLO", xi = event, ni = total,
                              data = dat, method = "ML", test = "z",
                              level = ma$level.ma * 100)
  common <- metafor::rma.glmm(measure = "PLO", xi = event, ni = total,
                              data = dat, method = "FE", test = "z",
                              level = ma$level.ma * 100)
  for (model in c("common", "random")) {
    reference <- if (model == "common") common else random
    expect_equal(ma[[paste0("TE.", model)]], as.numeric(coef(reference)))
    expect_equal(ma[[paste0("seTE.", model)]], reference$se)
    expect_equal(ma[[paste0("lower.", model)]], reference$ci.lb)
    expect_equal(ma[[paste0("upper.", model)]], reference$ci.ub)
    expect_equal(ma[[paste0("statistic.", model)]], reference$zval)
    expect_equal(ma[[paste0("pval.", model)]], reference$pval)
  }
  expect_identical(ma$method.tau, "ML")
  expect_identical(ma$method.I2, "tau2")
  expect_true(is.infinite(ma$df.random))
  expect_equal(ma$tau2, random$tau2)
  expect_equal(ma$se.tau2, random$se.tau2)
  expect_equal(unname(ma$Q), c(random$QE.Wld, random$QE.LRT))
  expect_equal(unname(ma$df.Q), rep(random$QE.df, 2))
  expect_equal(unname(ma$pval.Q),
               pchisq(c(random$QE.Wld, random$QE.LRT), random$QE.df,
                      lower.tail = FALSE))
  expect_equal(ma$I2, random$I2 / 100)
})
