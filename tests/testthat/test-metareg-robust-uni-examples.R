test_that("metareg robust methods work without cluster", {
  skip_if_not_installed("metadat")

  data(dat.konstantopoulos2011, package = "metadat")
  dat <- dat.konstantopoulos2011
  dat$group <-
    factor(ifelse(dat$district <= median(dat$district), "A", "B"))

  ma <- metagen(yi, sqrt(vi), studlab = study,
                data = dat,
                sm = "SMD")

  mr.cr0 <- metareg(ma, ~ group, method.random.ci = "CR0")
  mr.cr1 <- metareg(ma, ~ group, method.random.ci = "CR1")

  rma <- metafor::rma.uni(yi = yi, sei = sqrt(vi), mods = ~ group,
                          data = dat, method = ma$method.tau, test = "z",
                          level = ma$level.ma * 100)
  rma_cr0 <- metafor::robust(rma, cluster = dat$study,
                             adjust = FALSE, clubSandwich = FALSE)
  rma_cr1 <- metafor::robust(rma, cluster = dat$study,
                             adjust = TRUE, clubSandwich = FALSE)

  expect_s3_class(mr.cr0, "metareg")
  expect_s3_class(mr.cr1, "metareg")

  expect_equal(mr.cr0$.meta$method.random.ci, "CR0")
  expect_equal(mr.cr1$.meta$method.random.ci, "CR1")

  expect_equal(as.numeric(coef(mr.cr0)), as.numeric(coef(rma_cr0)))
  expect_equal(mr.cr0$se, rma_cr0$se)
  expect_equal(mr.cr0$ci.lb, rma_cr0$ci.lb)
  expect_equal(mr.cr0$ci.ub, rma_cr0$ci.ub)
  expect_equal(mr.cr0$pval, rma_cr0$pval)

  expect_equal(as.numeric(coef(mr.cr1)), as.numeric(coef(rma_cr1)))
  expect_equal(mr.cr1$se, rma_cr1$se)
  expect_equal(mr.cr1$ci.lb, rma_cr1$ci.lb)
  expect_equal(mr.cr1$ci.ub, rma_cr1$ci.ub)
  expect_equal(mr.cr1$pval, rma_cr1$pval)

  skip_if_not(suppressPackageStartupMessages(
    requireNamespace("clubSandwich", quietly = TRUE)),
    "clubSandwich not installed")

  mr.cr2 <- metareg(ma, ~ group, method.random.ci = "CR2")
  rma_cr2 <-
    suppressPackageStartupMessages(
      metafor::robust(rma, cluster = dat$study,
                      adjust = TRUE, clubSandwich = TRUE))

  expect_s3_class(mr.cr2, "metareg")

  expect_equal(mr.cr2$.meta$method.random.ci, "CR2")

  expect_equal(as.numeric(coef(mr.cr2)), as.numeric(coef(rma_cr2)))
  expect_equal(mr.cr2$se, rma_cr2$se)
  expect_equal(mr.cr2$ci.lb, rma_cr2$ci.lb)
  expect_equal(mr.cr2$ci.ub, rma_cr2$ci.ub)
  expect_equal(mr.cr2$pval, rma_cr2$pval)
})
