test_that("metagen three-level robust example runs", {
  skip_if_not_installed("metadat")
  
  data(dat.konstantopoulos2011, package = "metadat")
  dat <- dat.konstantopoulos2011
  dat$idx <- seq_len(nrow(dat))
  
  ma <- metagen(yi, sqrt(vi), studlab = study,
                data = dat,
                sm = "SMD",
                cluster = district,
                detail.tau = c("district", "district/school"))
  
  expect_true(ma$three.level)
  expect_equal(ma$detail.tau, c("district", "district/school"))
  
  ma_cr01 <- update(ma, method.random.ci = c("classic", "CR0", "CR1"))
  expect_equal(ma_cr01$method.random.ci, c("classic", "CR0", "CR1"))
  expect_true(all(is.finite(ma_cr01$seTE.random)))
  
  expect_true(is.na(ma$TE.common))
  expect_true(is.na(ma$seTE.common))
  expect_true(is.na(ma$lower.common))
  expect_true(is.na(ma$upper.common))
  expect_true(is.na(ma$pval.common))
  expect_true(is.na(ma$statistic.common))
  
  V <- metafor::vcalc(vi = dat$vi, cluster = dat$district,
                      obs = dat$idx, rho = 0)
  rma <- metafor::rma.mv(yi, V = V, random = ~ 1 | district / idx,
                         data = dat, method = ma$method.tau,
                         test = "z", level = ma$level.ma * 100)
  rma_cr0 <- metafor::robust(rma, cluster = dat$district,
                             adjust = FALSE, clubSandwich = FALSE)
  rma_cr1 <- metafor::robust(rma, cluster = dat$district,
                             adjust = TRUE, clubSandwich = FALSE)
  
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
               c(rma$ddf, rma_cr0$ddf, rma_cr1$ddf))
  expect_equal(unname(ma_cr01$tau2), unname(rma$sigma2))
  
  skip_if_not(suppressPackageStartupMessages(
    requireNamespace("clubSandwich", quietly = TRUE)),
    "clubSandwich not installed")
  
  ma_cr2 <- update(ma, method.random.ci = "CR2")
  
  rma_cr2 <-
    suppressPackageStartupMessages(
      metafor::robust(rma, cluster = dat$district,
                      adjust = TRUE, clubSandwich = TRUE))
  
  expect_equal(ma_cr2$method.random.ci, "CR2")
  expect_true(all(is.finite(ma_cr2$seTE.random)))
  expect_true(is.finite(ma_cr2$df.random))
  
  expect_equal(unname(ma_cr2$tau2), unname(rma$sigma2))
  
  expect_equal(unname(ma_cr2$TE.random), as.numeric(coef(rma)))
  expect_equal(unname(ma_cr2$seTE.random), rma_cr2$se)
  expect_equal(unname(ma_cr2$lower.random), rma_cr2$ci.lb)
  expect_equal(unname(ma_cr2$upper.random), rma_cr2$ci.ub)
  expect_equal(unname(ma_cr2$pval.random), rma_cr2$pval)
  expect_equal(unname(ma_cr2$statistic.random), rma_cr2$zval)
  expect_equal(unname(ma_cr2$df.random), rma_cr2$ddf)
})
