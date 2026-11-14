test_that("small three-level metagen matches direct robust inference", {
  dat <- external_multilevel_fixture()
  ma <- metagen(yi, sqrt(vi), data = dat, cluster = district,
                method.tau = "REML", method.tau.ci = "",
                method.random.ci = c("classic", "CR1"),
                detail.tau = c("district", "district/effect"))
  V <- metafor::vcalc(vi = dat$vi, cluster = dat$district,
                      obs = dat$idx, rho = 0)
  direct <- metafor::rma.mv(yi, V = V, random = ~ 1 | district / idx,
                            data = dat, method = "REML", test = "z",
                            level = ma$level.ma * 100)
  robust <- metafor::robust(direct, cluster = dat$district, adjust = TRUE,
                            clubSandwich = FALSE)

  expect_true(ma$three.level)
  expect_identical(ma$method.random.ci, c("classic", "CR1"))
  expect_equal(ma$detail.tau, c("district", "district/effect"))
  expect_true(all(direct$sigma2 > 0))
  expect_equal(unname(ma$tau2), direct$sigma2)
  expect_equal(unname(ma$TE.random), rep(as.numeric(coef(direct)), 2))
  expect_equal(unname(ma$seTE.random), c(direct$se, robust$se))
  expect_equal(unname(ma$lower.random), c(direct$ci.lb, robust$ci.lb))
  expect_equal(unname(ma$upper.random), c(direct$ci.ub, robust$ci.ub))
  expect_equal(unname(ma$statistic.random), c(direct$zval, robust$zval))
  expect_equal(unname(ma$pval.random), c(direct$pval, robust$pval))
  expect_equal(unname(ma$df.random), c(direct$ddf, robust$ddf))
  expect_true(all(is.na(c(ma$TE.common, ma$seTE.common, ma$lower.common,
                          ma$upper.common, ma$pval.common, ma$statistic.common))))
})
