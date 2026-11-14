test_that("compact moderator regression matches rma.uni", {
  dat <- data.frame(yi = c(-0.4, 0.2, -0.1, 0.6, 0.3, 0.9),
                    se = c(0.2, 0.3, 0.25, 0.2, 0.35, 0.3),
                    region = factor(rep(c("A", "B"), each = 3)))
  ma <- metagen(yi, se, data = dat, method.tau = "REML")
  mr <- metareg(ma, ~ region)
  reference <- metafor::rma.uni(yi, sei = se, mods = ~ region, data = dat,
                               method = "REML", test = "z",
                               level = ma$level.ma * 100)
  expect_s3_class(mr, "metareg")
  expect_equal(as.numeric(coef(mr)), as.numeric(coef(reference)))
  fields <- c("tau2", "se", "ci.lb", "ci.ub", "zval", "pval", "QM", "QMp")
  expect_equal(mr[fields], reference[fields])
})
