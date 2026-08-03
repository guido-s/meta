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
  fit <- metafor::rma.mv(yi = yi, V = vi, random = ~ 1 | .id / .idx,
                         mods = ~ group,
                         data = transform(dat, .id = district, .idx = idx),
                         method = ma$method.tau, test = "z",
                         level = ma$level.ma * 100)
  fit.cr1 <- metafor::robust(fit, cluster = dat$district,
                             adjust = TRUE, clubSandwich = FALSE)

  expect_s3_class(mr, "metareg")
  expect_equal(mr$.meta$method.random.ci, "CR1")
  expect_equal(as.numeric(coef(mr)), as.numeric(coef(fit.cr1)))
  expect_equal(mr$se, fit.cr1$se)
  expect_equal(mr$ci.lb, fit.cr1$ci.lb)
  expect_equal(mr$ci.ub, fit.cr1$ci.ub)
  expect_equal(mr$pval, fit.cr1$pval)

  skip_if_not(suppressPackageStartupMessages(
    requireNamespace("clubSandwich", quietly = TRUE)),
    "clubSandwich not installed")

  mr.cr2 <- metareg(ma, ~ group, method.random.ci = "CR2")
  fit.cr2 <-
    suppressPackageStartupMessages(
      metafor::robust(fit, cluster = dat$district,
                      adjust = TRUE, clubSandwich = TRUE))

  expect_s3_class(mr.cr2, "metareg")
  expect_equal(mr.cr2$.meta$method.random.ci, "CR2")
  expect_equal(as.numeric(coef(mr.cr2)), as.numeric(coef(fit.cr2)))
  expect_equal(mr.cr2$se, fit.cr2$se)
  expect_equal(mr.cr2$ci.lb, fit.cr2$ci.lb)
  expect_equal(mr.cr2$ci.ub, fit.cr2$ci.ub)
  expect_equal(mr.cr2$pval, fit.cr2$pval)
})
