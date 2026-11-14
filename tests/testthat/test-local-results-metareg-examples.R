test_that("metareg examples run", {
  data(Fleiss1993cont)
  
  Fleiss1993cont$age <- c(55, 65, 55, 65, 55)
  Fleiss1993cont$region <- c("Europe", "Europe", "Asia", "Asia", "Europe")
  
  ma1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
                  data = Fleiss1993cont, sm = "SMD")
  
  ma2 <- update(ma1, subgroup = region, tau.common = TRUE, common = FALSE)
  mr2 <- metareg(ma2)
  
  expect_s3_class(mr2, "metareg")
  
  rma2 <- metafor::rma.uni(yi = ma2$TE, sei = ma2$seTE,
                           mods = ~ .subgroup, data = ma2$data,
                           method = ma2$method.tau, test = "z",
                           level = ma2$level.ma * 100)
  
  expect_equal(as.numeric(coef(mr2)), as.numeric(coef(rma2)))
  expect_equal(mr2$se, rma2$se)
  expect_equal(mr2$ci.lb, rma2$ci.lb)
  expect_equal(mr2$ci.ub, rma2$ci.ub)
  expect_equal(mr2$pval, rma2$pval)
  
  mr2_noint <- metareg(ma2, intercept = FALSE)
  
  expect_s3_class(mr2_noint, "metareg")
  
  rma2.no.intercept <- metafor::rma.uni(yi = ma2$TE, sei = ma2$seTE,
                                        mods = ~ .subgroup - 1, data = ma2$data,
                                        method = ma2$method.tau, test = "z",
                                        level = ma2$level.ma * 100)
  
  expect_equal(as.numeric(coef(mr2_noint)), as.numeric(coef(rma2.no.intercept)))
  expect_equal(mr2_noint$se, rma2.no.intercept$se)
  expect_equal(mr2_noint$ci.lb, rma2.no.intercept$ci.lb)
  expect_equal(mr2_noint$ci.ub, rma2.no.intercept$ci.ub)
  expect_equal(mr2_noint$pval, rma2.no.intercept$pval)
  
  mr1_region <- metareg(ma1, region)
  
  expect_s3_class(mr1_region, "metareg")
  
  rma1_region <- metafor::rma.uni(yi = ma1$TE, sei = ma1$seTE,
                                  mods = ~ region, data = ma1$data,
                                  method = ma1$method.tau, test = "z",
                                  level = ma1$level.ma * 100)
  
  expect_equal(as.numeric(coef(mr1_region)), as.numeric(coef(rma1_region)))
  expect_equal(mr1_region$se, rma1_region$se)
  expect_equal(mr1_region$ci.lb, rma1_region$ci.lb)
  expect_equal(mr1_region$ci.ub, rma1_region$ci.ub)
  expect_equal(mr1_region$pval, rma1_region$pval)
  
  ma3 <- update(ma1, subgroup = region)
  mr3 <- metareg(ma3, region + age)
  
  expect_s3_class(mr3, "metareg")
  
  rma3 <- metafor::rma.uni(yi = ma3$TE, sei = ma3$seTE,
                           mods = ~ region + age, data = ma3$data,
                           method = ma3$method.tau, test = "z",
                           level = ma3$level.ma * 100)
  
  expect_equal(as.numeric(coef(mr3)), as.numeric(coef(rma3)))
  expect_equal(mr3$se, rma3$se)
  expect_equal(mr3$ci.lb, rma3$ci.lb)
  expect_equal(mr3$ci.ub, rma3$ci.ub)
  expect_equal(mr3$pval, rma3$pval)
  
  mr1_formula <- metareg(ma1, ~ region)
  
  expect_s3_class(mr1_formula, "metareg")
  
  expect_equal(as.numeric(coef(mr1_formula)), as.numeric(coef(rma1_region)))
  expect_equal(mr1_formula$se, rma1_region$se)
  expect_equal(mr1_formula$ci.lb, rma1_region$ci.lb)
  expect_equal(mr1_formula$ci.ub, rma1_region$ci.ub)
  expect_equal(mr1_formula$pval, rma1_region$pval)
  
  mr3.formula <- metareg(ma3, ~ region + age)
  
  expect_s3_class(mr3.formula, "metareg")
  
  expect_equal(as.numeric(coef(mr3.formula)), as.numeric(coef(rma3)))
  expect_equal(mr3.formula$se, rma3$se)
  expect_equal(mr3.formula$ci.lb, rma3$ci.lb)
  expect_equal(mr3.formula$ci.ub, rma3$ci.ub)
  expect_equal(mr3.formula$pval, rma3$pval)
})
