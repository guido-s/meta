test_that("small three-level metareg matches direct CR2 inference", {
  skip_if_not_installed("clubSandwich")
  dat <- external_multilevel_fixture()
  ma <- metagen(yi, sqrt(vi), data = dat, cluster = district,
                method.tau = "REML", method.tau.ci = "")
  mr <- metareg(ma, ~ group, method.random.ci = "CR2")
  V <- metafor::vcalc(vi = dat$vi, cluster = dat$district,
                      obs = dat$idx, rho = 0)
  direct <- metafor::rma.mv(yi, V = V, mods = ~ group,
                            random = ~ 1 | district / idx, data = dat,
                            method = "REML", test = "z",
                            level = ma$level.ma * 100)
  reference <- suppressPackageStartupMessages(
    metafor::robust(direct, cluster = dat$district,
                    adjust = TRUE, clubSandwich = TRUE))
  expect_s3_class(mr, "metareg")
  expect_identical(mr$.meta$method.random.ci, "CR2")
  expect_true(all(direct$sigma2 > 0))
  expect_equal(as.numeric(coef(mr)), as.numeric(coef(reference)))
  fields <- c("sigma2", "se", "vb", "ci.lb", "ci.ub", "zval", "pval",
              "ddf", "QM", "QMdf", "QMp")
  expect_equal(mr[fields], reference[fields])
  expect_equal(mr$V, V)
})
