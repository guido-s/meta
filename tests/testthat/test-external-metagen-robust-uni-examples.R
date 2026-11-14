test_that("metagen robust methods work without cluster", {
  skip_if_not_installed("metadat")

  data(dat.konstantopoulos2011, package = "metadat")
  dat <- dat.konstantopoulos2011

  ma <- metagen(yi, sqrt(vi), studlab = study,
                data = dat,
                sm = "SMD")
  ma_cr01 <- update(ma, method.random.ci = c("classic", "CR0", "CR1"))

  rma <- metafor::rma.uni(yi = dat$yi, sei = sqrt(dat$vi),
                          method = ma$method.tau, test = "z",
                          level = ma$level.ma * 100)
  rma_cr0 <- metafor::robust(rma, cluster = dat$study,
                             adjust = FALSE, clubSandwich = FALSE)
  rma_cr1 <- metafor::robust(rma, cluster = dat$study,
                             adjust = TRUE, clubSandwich = FALSE)

  expect_false(ma_cr01$three.level)
  expect_equal(ma_cr01$method.random.ci, c("classic", "CR0", "CR1"))

  expect_equal(unname(ma_cr01$TE.random), rep(as.numeric(coef(rma)), 3))
  expect_equal(unname(ma_cr01$seTE.random),
               c(rma$se, rma_cr0$se, rma_cr1$se))
  expect_equal(unname(ma_cr01$lower.random),
               c(rma$ci.lb, rma_cr0$ci.lb, rma_cr1$ci.lb))
  expect_equal(unname(ma_cr01$upper.random),
               c(rma$ci.ub, rma_cr0$ci.ub, rma_cr1$ci.ub))
  expect_equal(unname(ma_cr01$pval.random),
               c(rma$pval, rma_cr0$pval, rma_cr1$pval))
  expect_equal(unname(ma_cr01$statistic.random),
               c(rma$zval, rma_cr0$zval, rma_cr1$zval))
  expect_equal(unname(ma_cr01$df.random),
               c(Inf, rma_cr0$ddf, rma_cr1$ddf))
  expect_equal(unname(ma_cr01$tau2), unname(rma$tau2))

  skip_if_not(suppressPackageStartupMessages(
    requireNamespace("clubSandwich", quietly = TRUE)),
    "clubSandwich not installed")

  ma_cr2 <- update(ma, method.random.ci = "CR2")
  rma_cr2 <-
    suppressPackageStartupMessages(
      metafor::robust(rma, cluster = dat$study,
                      adjust = TRUE, clubSandwich = TRUE))

  expect_equal(ma_cr2$method.random.ci, "CR2")
  expect_equal(unname(ma_cr2$TE.random), as.numeric(coef(rma)))
  expect_equal(unname(ma_cr2$seTE.random), rma_cr2$se)
  expect_equal(unname(ma_cr2$lower.random), rma_cr2$ci.lb)
  expect_equal(unname(ma_cr2$upper.random), rma_cr2$ci.ub)
  expect_equal(unname(ma_cr2$pval.random), rma_cr2$pval)
  expect_equal(unname(ma_cr2$statistic.random), rma_cr2$zval)
  expect_equal(unname(ma_cr2$df.random), rma_cr2$ddf)
  expect_equal(unname(ma_cr2$tau2), unname(rma$tau2))
})
