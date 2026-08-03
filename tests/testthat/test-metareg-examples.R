test_that("metareg examples run", {
  data(Fleiss1993cont)

  Fleiss1993cont$age <- c(55, 65, 55, 65, 55)
  Fleiss1993cont$region <- c("Europe", "Europe", "Asia", "Asia", "Europe")

  ma1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
                  data = Fleiss1993cont, sm = "SMD")

  ma2 <- update(ma1, subgroup = region, tau.common = TRUE, common = FALSE)
  mr2 <- metareg(ma2)
  expect_s3_class(mr2, "metareg")
  fit2 <- metafor::rma.uni(yi = ma2$TE, sei = ma2$seTE,
                           mods = ~ .subgroup, data = ma2$data,
                           method = ma2$method.tau, test = "z",
                           level = ma2$level.ma * 100)
  expect_equal(as.numeric(coef(mr2)), as.numeric(coef(fit2)))
  expect_equal(mr2$se, fit2$se)
  expect_equal(mr2$ci.lb, fit2$ci.lb)
  expect_equal(mr2$ci.ub, fit2$ci.ub)
  expect_equal(mr2$pval, fit2$pval)

  mr2.no.intercept <- metareg(ma2, intercept = FALSE)
  expect_s3_class(mr2.no.intercept, "metareg")
  fit2.no.intercept <- metafor::rma.uni(yi = ma2$TE, sei = ma2$seTE,
                                        mods = ~ .subgroup - 1, data = ma2$data,
                                        method = ma2$method.tau, test = "z",
                                        level = ma2$level.ma * 100)
  expect_equal(as.numeric(coef(mr2.no.intercept)),
               as.numeric(coef(fit2.no.intercept)))
  expect_equal(mr2.no.intercept$se, fit2.no.intercept$se)
  expect_equal(mr2.no.intercept$ci.lb, fit2.no.intercept$ci.lb)
  expect_equal(mr2.no.intercept$ci.ub, fit2.no.intercept$ci.ub)
  expect_equal(mr2.no.intercept$pval, fit2.no.intercept$pval)

  mr1.region <- metareg(ma1, region)
  expect_s3_class(mr1.region, "metareg")
  fit1.region <- metafor::rma.uni(yi = ma1$TE, sei = ma1$seTE,
                                  mods = ~ region, data = ma1$data,
                                  method = ma1$method.tau, test = "z",
                                  level = ma1$level.ma * 100)
  expect_equal(as.numeric(coef(mr1.region)), as.numeric(coef(fit1.region)))
  expect_equal(mr1.region$se, fit1.region$se)
  expect_equal(mr1.region$ci.lb, fit1.region$ci.lb)
  expect_equal(mr1.region$ci.ub, fit1.region$ci.ub)
  expect_equal(mr1.region$pval, fit1.region$pval)

  ma3 <- update(ma1, subgroup = region)
  mr3 <- metareg(ma3, region + age)
  expect_s3_class(mr3, "metareg")
  fit3 <- metafor::rma.uni(yi = ma3$TE, sei = ma3$seTE,
                           mods = ~ region + age, data = ma3$data,
                           method = ma3$method.tau, test = "z",
                           level = ma3$level.ma * 100)
  expect_equal(as.numeric(coef(mr3)), as.numeric(coef(fit3)))
  expect_equal(mr3$se, fit3$se)
  expect_equal(mr3$ci.lb, fit3$ci.lb)
  expect_equal(mr3$ci.ub, fit3$ci.ub)
  expect_equal(mr3$pval, fit3$pval)

  mr1.formula <- metareg(ma1, ~ region)
  expect_s3_class(mr1.formula, "metareg")
  expect_equal(as.numeric(coef(mr1.formula)), as.numeric(coef(fit1.region)))
  expect_equal(mr1.formula$se, fit1.region$se)
  expect_equal(mr1.formula$ci.lb, fit1.region$ci.lb)
  expect_equal(mr1.formula$ci.ub, fit1.region$ci.ub)
  expect_equal(mr1.formula$pval, fit1.region$pval)

  mr3.formula <- metareg(ma3, ~ region + age)
  expect_s3_class(mr3.formula, "metareg")
  expect_equal(as.numeric(coef(mr3.formula)), as.numeric(coef(fit3)))
  expect_equal(mr3.formula$se, fit3$se)
  expect_equal(mr3.formula$ci.lb, fit3$ci.lb)
  expect_equal(mr3.formula$ci.ub, fit3$ci.ub)
  expect_equal(mr3.formula$pval, fit3$pval)
})
