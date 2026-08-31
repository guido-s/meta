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

test_that("metagen recalculates supplied confidence limits for changed level", {
  catch_warnings <- function(expr) {
    warnings <- character(0)
    res <- withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(res = res, warnings = warnings)
  }

  dat <- data.frame(
    study = paste("study", 1:5),
    MD = 1:5,
    seMD = c(NA, 2, NA, 4, NA),
    mylower = c(0, NA, 1, NA, 0),
    myupper = c(2, NA, 5, NA, 8),
    mylevel = c(0.95, NA, 0.95, NA, 0.9)
  )

  res1 <- catch_warnings(
    metagen(MD, seMD, study, lower = mylower, upper = myupper,
            data = dat, level = 0.95, level.ci = 0.95)
  )
  ma1 <- res1$res
  res2 <- catch_warnings(
    metagen(MD, seMD, study, lower = mylower, upper = myupper,
            data = dat, level = 0.90, level.ci = 0.95)
  )
  ma2 <- res2$res
  res3 <- catch_warnings(
    metagen(MD, seMD, study, lower = mylower, upper = myupper,
            data = dat, level = 0.95, level.ci = mylevel)
  )
  ma3 <- res3$res

  warn.asymmetric <- grep("not halfway", res1$warnings, value = TRUE)
  expect_length(warn.asymmetric, 1)
  expect_true(grepl("'study 5'", warn.asymmetric))
  expect_false(grepl("'study 1'|'study 3'", warn.asymmetric))
  expect_true(any(grepl("not halfway.*'study 5'", res2$warnings)))
  expect_true(any(grepl("not halfway.*'study 5'", res3$warnings)))
  expect_true(any(grepl("following studies.*'study 1'.*'study 3'.*'study 5'",
                        res2$warnings)))
  expect_true(any(grepl("following study.*'study 5'", res3$warnings)))

  expect_equal(ma1$lower[c(1, 3, 5)], dat$mylower[c(1, 3, 5)])
  expect_equal(ma1$upper[c(1, 3, 5)], dat$myupper[c(1, 3, 5)])

  expect_false(
    isTRUE(all.equal(ma1$lower[c(1, 3, 5)], ma2$lower[c(1, 3, 5)]))
  )
  expect_false(
    isTRUE(all.equal(ma1$upper[c(1, 3, 5)], ma2$upper[c(1, 3, 5)]))
  )

  expect_equal(ma1$lower[c(1, 3)], ma3$lower[c(1, 3)])
  expect_equal(ma1$upper[c(1, 3)], ma3$upper[c(1, 3)])
  expect_false(isTRUE(all.equal(ma1$lower[5], ma3$lower[5])))
  expect_false(isTRUE(all.equal(ma1$upper[5], ma3$upper[5])))
})

test_that("metagen rejects n.c for single-arm summary measures", {
  TE <- c(0.1, 0.2, 0.3)
  seTE <- c(0.01, 0.02, 0.03)
  n <- c(10, 20, 30)

  for (sm in c("PRAW", "PLOGIT", "PLN", "PAS", "PFT",
               "IR", "IRLN", "IRS", "IRFT",
               "MRAW", "MLN",
               "COR", "ZCOR")) {
    expect_error(
      metagen(TE, seTE, sm = sm, n.c = n),
      "Argument 'n.c' must not be provided"
    )
  }
})

test_that("metagen uses n.e as sample sizes for untransformed PFT input", {
  p <- c(0.1, 0.2, 0.3)
  lower <- c(0.05, 0.10, 0.20)
  upper <- c(0.20, 0.30, 0.40)
  seTE <- c(0.01, 0.02, 0.03)
  n <- c(10, 20, 30)

  m <- metagen(p, seTE, lower = lower, upper = upper,
               sm = "PFT", transf = FALSE, n.e = n, warn = FALSE)

  expect_equal(m$TE, p2asin(p, n))
  expect_equal(m$lower, p2asin(lower, n))
  expect_equal(m$upper, p2asin(upper, n))
  expect_error(
    metagen(p, seTE, sm = "PFT", transf = FALSE, n.e = n, n.c = n),
    "Argument 'n.c' must not be provided"
  )
})
