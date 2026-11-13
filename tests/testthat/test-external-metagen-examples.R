test_that("metagen Fleiss1993bin example matches rma.uni", {
  data(Fleiss1993bin)
  
  ma_bin <- metabin(d.asp, n.asp, d.plac, n.plac, study,
                    data = Fleiss1993bin, sm = "RR", method = "I")
  ma <- metagen(TE, seTE, studlab, n.e = n.e + n.c,
                data = ma_bin, sm = "RR")
  
  rma_c <- metafor::rma.uni(yi = ma$TE, sei = ma$seTE,
                            method = "FE", test = "z",
                            level = ma$level.ma * 100)
  rma_r <- metafor::rma.uni(yi = ma$TE, sei = ma$seTE,
                            method = ma$method.tau, test = "z",
                            level = ma$level.ma * 100)
  
  expect_equal(ma$TE.common, as.numeric(coef(rma_c)))
  expect_equal(ma$seTE.common, rma_c$se)
  expect_equal(ma$lower.common, rma_c$ci.lb)
  expect_equal(ma$upper.common, rma_c$ci.ub)
  expect_equal(ma$pval.common, rma_c$pval)
  expect_equal(ma$statistic.common, rma_c$zval)
  
  expect_equal(ma$tau2, rma_r$tau2)
  
  expect_equal(ma$TE.random, as.numeric(coef(rma_r)))
  expect_equal(ma$seTE.random, rma_r$se)
  expect_equal(ma$lower.random, rma_r$ci.lb)
  expect_equal(ma$upper.random, rma_r$ci.ub)
  expect_equal(ma$pval.random, rma_r$pval)
  expect_equal(ma$statistic.random, rma_r$zval)
  
  expect_true(is.infinite(ma$df.random))
  expect_true(is.na(rma_r$ddf))
})
