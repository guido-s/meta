test_that("metagen three-level robust example runs", {
  skip_if_not_installed("metadat")

  data(dat.konstantopoulos2011, package = "metadat")
  dat <- dat.konstantopoulos2011
  dat$idx <- seq_len(nrow(dat))

  ma3 <- metagen(yi, sqrt(vi), studlab = study,
                 data = dat,
                 sm = "SMD",
                 cluster = district,
                 detail.tau = c("district", "district/school"))

  expect_true(ma3$three.level)
  expect_equal(ma3$detail.tau, c("district", "district/school"))

  ma3.cr <- update(ma3, method.random.ci = c("classic", "CR0", "CR1"))
  expect_equal(ma3.cr$method.random.ci, c("classic", "CR0", "CR1"))
  expect_true(all(is.finite(ma3.cr$seTE.random)))

  expect_true(is.na(ma3$TE.common))
  expect_true(is.na(ma3$seTE.common))
  expect_true(is.na(ma3$lower.common))
  expect_true(is.na(ma3$upper.common))
  expect_true(is.na(ma3$pval.common))
  expect_true(is.na(ma3$statistic.common))

  V <- metafor::vcalc(vi = dat$vi, cluster = dat$district,
                      obs = dat$idx, rho = 0)
  fit.z <- metafor::rma.mv(yi, V = V, random = ~ 1 | district / idx,
                           data = dat, method = ma3$method.tau,
                           test = "z", level = ma3$level.ma * 100)
  fit.cr0 <- metafor::robust(fit.z, cluster = dat$district,
                             adjust = FALSE, clubSandwich = FALSE)
  fit.cr1 <- metafor::robust(fit.z, cluster = dat$district,
                             adjust = TRUE, clubSandwich = FALSE)

  expect_equal(unname(ma3.cr$TE.random), rep(as.numeric(coef(fit.z)), 3))
  expect_equal(unname(ma3.cr$seTE.random),
               c(fit.z$se, fit.cr0$se, fit.cr1$se))
  expect_equal(unname(ma3.cr$lower.random),
               c(fit.z$ci.lb, fit.cr0$ci.lb, fit.cr1$ci.lb))
  expect_equal(unname(ma3.cr$upper.random),
               c(fit.z$ci.ub, fit.cr0$ci.ub, fit.cr1$ci.ub))
  expect_equal(unname(ma3.cr$pval.random),
               c(fit.z$pval, fit.cr0$pval, fit.cr1$pval))
  expect_equal(unname(ma3.cr$statistic.random),
               c(fit.z$zval, fit.cr0$zval, fit.cr1$zval))
  expect_equal(unname(ma3.cr$df.random),
               c(fit.z$ddf, fit.cr0$ddf, fit.cr1$ddf))
  expect_equal(unname(ma3.cr$tau2), unname(fit.z$sigma2))

  skip_if_not(suppressPackageStartupMessages(
    requireNamespace("clubSandwich", quietly = TRUE)),
    "clubSandwich not installed")

  ma3.cr2 <- update(ma3, method.random.ci = c("classic", "CR0", "CR1", "CR2"))
  fit.cr2 <-
    suppressPackageStartupMessages(
      metafor::robust(fit.z, cluster = dat$district,
                      adjust = TRUE, clubSandwich = TRUE))

  expect_equal(ma3.cr2$method.random.ci, c("classic", "CR0", "CR1", "CR2"))
  expect_true(all(is.finite(ma3.cr2$seTE.random)))
  expect_true(is.finite(ma3.cr2$df.random["CR2"]))
  expect_equal(unname(ma3.cr2$TE.random), rep(as.numeric(coef(fit.z)), 4))
  expect_equal(unname(ma3.cr2$seTE.random),
               c(fit.z$se, fit.cr0$se, fit.cr1$se, fit.cr2$se))
  expect_equal(unname(ma3.cr2$lower.random),
               c(fit.z$ci.lb, fit.cr0$ci.lb, fit.cr1$ci.lb, fit.cr2$ci.lb))
  expect_equal(unname(ma3.cr2$upper.random),
               c(fit.z$ci.ub, fit.cr0$ci.ub, fit.cr1$ci.ub, fit.cr2$ci.ub))
  expect_equal(unname(ma3.cr2$pval.random),
               c(fit.z$pval, fit.cr0$pval, fit.cr1$pval, fit.cr2$pval))
  expect_equal(unname(ma3.cr2$statistic.random),
               c(fit.z$zval, fit.cr0$zval, fit.cr1$zval, fit.cr2$zval))
  expect_equal(unname(ma3.cr2$df.random),
               c(fit.z$ddf, fit.cr0$ddf, fit.cr1$ddf, fit.cr2$ddf))
  expect_equal(unname(ma3.cr2$tau2), unname(fit.z$sigma2))
})
