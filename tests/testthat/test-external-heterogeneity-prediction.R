test_that("univariate heterogeneity and normal prediction match metafor", {
  TE <- c(-0.8, -0.35, 0.1, 0.55, 1.05, 1.4)
  se <- c(0.12, 0.21, 0.16, 0.28, 0.18, 0.24)
  fit <- metafor::rma.uni(yi = TE, sei = se, method = "REML",
                          test = "z", level = 90)

  for (method in c("QP", "BJ", "J")) {
    ma <- metagen(TE, se, common = FALSE, method.tau = "REML",
                  method.tau.ci = method, method.I2 = "tau",
                  level.hetstat = 0.90, method.random.ci = "classic",
                  prediction = TRUE, method.predict = "S",
                  level.predict = 0.90)
    ci.fit <- if (method == "QP") fit else
      metafor::rma.uni(yi = TE, sei = se, method = "GENQ",
                       weights = if (method == "BJ") 1 / se^2 else 1 / se,
                       test = "z", level = 90)
    ci <- metafor::confint.rma.uni(ci.fit, level = 90)$random

    expect_identical(ma$method.tau, "REML")
    expect_identical(ma$method.tau.ci, method)
    # "tau" is accepted by partial matching to the canonical "tau2".
    expect_identical(ma$method.I2, "tau2")
    expect_equal(ma$tau2, fit$tau2)
    expect_equal(ma$se.tau2, fit$se.tau2)
    expect_equal(ma$tau, sqrt(fit$tau2))
    expect_equal(c(ma$lower.tau2, ma$upper.tau2),
                 unname(ci["tau^2", c("ci.lb", "ci.ub")]))
    expect_equal(c(ma$lower.tau, ma$upper.tau),
                 unname(ci["tau", c("ci.lb", "ci.ub")]))
    expect_equal(ma$I2, fit$I2 / 100)
    expect_equal(ma$H, sqrt(fit$H2))
    expect_equal(c(ma$lower.I2, ma$upper.I2),
                 unname(ci["I^2(%)", c("ci.lb", "ci.ub")]) / 100)
    expect_equal(c(ma$lower.H, ma$upper.H),
                 sqrt(unname(ci["H^2", c("ci.lb", "ci.ub")])))

    if (method == "QP") {
      # predict.rma is imported but unused by meta: this is an independent oracle.
      prediction <- metafor::predict.rma(fit, level = 90)
      expect_identical(ma$method.predict, "S")
      expect_equal(ma$df.predict, Inf)
      expect_equal(ma$lower.predict, as.numeric(prediction$pi.lb))
      expect_equal(ma$upper.predict, as.numeric(prediction$pi.ub))
    }
  }
})

test_that("clustered profile-likelihood heterogeneity matches metafor", {
  dat <- data.frame(
    yi = c(-1.25, -0.65, -0.80, -0.55, 0.15, -0.05,
           0.10, 0.80, 0.60, 0.85, 1.55, 1.30),
    se = rep(c(0.12, 0.16, 0.14), 4),
    cluster = rep(seq_len(4), each = 3),
    idx = seq_len(12))
  ma <- metagen(yi, se, data = dat, cluster = cluster, rho = 0,
                common = FALSE, method.tau = "REML", method.tau.ci = "PL",
                method.random.ci = "classic", level.hetstat = 0.90)
  V <- metafor::vcalc(vi = dat$se^2, cluster = dat$cluster,
                      obs = dat$idx, rho = 0)
  fit <- metafor::rma.mv(yi, V = V, random = ~ 1 | cluster / idx,
                         data = dat, method = "REML", test = "z", level = 90)
  ci <- metafor::confint.rma.mv(fit, level = 90)

  expect_true(ma$three.level)
  expect_identical(ma$method.tau.ci, "PL")
  expect_true(all(fit$sigma2 > 0))
  expect_equal(unname(ma$tau2), fit$sigma2)
  expect_equal(unname(ma$tau), sqrt(fit$sigma2))
  expect_equal(unname(ma$lower.tau2),
               c(ci[[1]]$random["sigma^2.1", "ci.lb"],
                 ci[[2]]$random["sigma^2.2", "ci.lb"]))
  expect_equal(unname(ma$upper.tau2),
               c(ci[[1]]$random["sigma^2.1", "ci.ub"],
                 ci[[2]]$random["sigma^2.2", "ci.ub"]))
  expect_equal(unname(ma$lower.tau),
               c(ci[[1]]$random["sigma.1", "ci.lb"],
                 ci[[2]]$random["sigma.2", "ci.lb"]))
  expect_equal(unname(ma$upper.tau),
               c(ci[[1]]$random["sigma.1", "ci.ub"],
                 ci[[2]]$random["sigma.2", "ci.ub"]))
})
