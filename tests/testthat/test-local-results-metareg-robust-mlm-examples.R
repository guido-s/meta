test_that("metareg three-level robust example runs", {
  skip_if_not_installed("metadat")

  data(dat.konstantopoulos2011, package = "metadat")
  dat <- dat.konstantopoulos2011
  dat$idx <- seq_len(nrow(dat))
  dat$group <-
    factor(ifelse(dat$district <= median(dat$district), "A", "B"))

  ma <- metagen(yi, sqrt(vi), studlab = study,
                data = dat,
                cluster = district)

  mr <- metareg(ma, ~ group, method.random.ci = "CR1")

  rma <- metafor::rma.mv(yi = yi, V = vi, random = ~ 1 | .id / .idx,
                         mods = ~ group,
                         data = transform(dat, .id = district, .idx = idx),
                         method = ma$method.tau, test = "z",
                         level = ma$level.ma * 100)

  rma_cr1 <- metafor::robust(rma, cluster = dat$district,
                             adjust = TRUE, clubSandwich = FALSE)

  expect_s3_class(mr, "metareg")

  expect_equal(mr$.meta$method.random.ci, "CR1")

  expect_equal(as.numeric(coef(mr)), as.numeric(coef(rma_cr1)))
  expect_equal(mr$se, rma_cr1$se)
  expect_equal(mr$ci.lb, rma_cr1$ci.lb)
  expect_equal(mr$ci.ub, rma_cr1$ci.ub)
  expect_equal(mr$pval, rma_cr1$pval)

  skip_if_not(suppressPackageStartupMessages(
    requireNamespace("clubSandwich", quietly = TRUE)),
    "clubSandwich not installed")

  mr.cr2 <- metareg(ma, ~ group, method.random.ci = "CR2")

  rma_cr2 <-
    suppressPackageStartupMessages(
      metafor::robust(rma, cluster = dat$district,
                      adjust = TRUE, clubSandwich = TRUE))

  expect_s3_class(mr.cr2, "metareg")

  expect_equal(mr.cr2$.meta$method.random.ci, "CR2")

  expect_equal(as.numeric(coef(mr.cr2)), as.numeric(coef(rma_cr2)))
  expect_equal(mr.cr2$se, rma_cr2$se)
  expect_equal(mr.cr2$ci.lb, rma_cr2$ci.lb)
  expect_equal(mr.cr2$ci.ub, rma_cr2$ci.ub)
  expect_equal(mr.cr2$pval, rma_cr2$pval)
})
