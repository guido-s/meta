test_that("small common-tau subgroup fit filters rho-aware covariance", {
  dat <- external_multilevel_fixture()
  dat <- dat[dat$district %in% c(1, 2, 4, 5), ]
  dat$idx <- seq_len(nrow(dat))
  dat$exclude <- dat$idx == 2L
  retained <- dat[!dat$exclude, ]
  ma <- metagen(yi, sqrt(vi), data = dat, cluster = district, rho = 0.7,
                exclude = exclude, subgroup = group, tau.common = TRUE,
                common = FALSE, method.tau = "REML", method.tau.ci = "",
                method.random.ci = "HK",
                keeprma = TRUE)
  V <- metafor::vcalc(vi = retained$vi, cluster = retained$district,
                      obs = retained$idx, rho = 0.7)
  reference <- metafor::rma.mv(yi, V = V, mods = ~ group - 1,
                              random = ~ 1 | district / idx, data = retained,
                              method = "REML", test = "t",
                              level = ma$level.ma * 100)
  expect_true(ma$three.level)
  expect_true(ma$tau.common)
  expect_equal(ma$rma.three.level$V, V)
  expect_equal(as.numeric(ma$rma.three.level$yi), retained$yi)
  expect_identical(ma$data$.idx[!dat$exclude], retained$idx)
  expect_identical(names(ma$TE.random.w), levels(retained$group))
  expect_equal(unname(ma$TE.random.w), as.numeric(coef(reference)))
  expect_equal(unname(ma$seTE.random.w), reference$se)
  expect_equal(unname(ma$lower.random.w), reference$ci.lb)
  expect_equal(unname(ma$upper.random.w), reference$ci.ub)
  expect_equal(unname(ma$statistic.random.w), reference$zval)
  expect_equal(unname(ma$pval.random.w), reference$pval)
  expect_equal(unname(ma$df.random.w), unname(reference$ddf))
  expect_equal(unname(ma$tau2.w), rep(sum(reference$sigma2), 2))
})
