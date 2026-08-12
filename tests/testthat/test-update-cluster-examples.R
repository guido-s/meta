test_that("update.meta resets tau methods when adding cluster", {
  dat <- data.frame(events = 1:5, N = 5:1 * 100, grp = c(1, 1, 2, 2, 3))

  ma1 <- metaprop(events, N, data = dat)
  ma2 <- metaprop(events, N, data = dat, cluster = grp)
  ma3 <- update(ma1, cluster = grp)

  expect_false(ma1$three.level)
  expect_true(ma2$three.level)
  expect_true(ma3$three.level)
  expect_equal(ma3$method.tau, ma2$method.tau)
  expect_equal(ma3$method.tau, "ML")
  expect_equal(ma3$method.tau.ci, ma2$method.tau.ci)
  expect_equal(ma3$method.tau.ci, "PL")

  ma4 <- update(ma1, cluster = grp, method.tau.ci = "")
  expect_equal(ma4$method.tau.ci, "")

  expect_warning(ma5 <- update(ma1, cluster = grp, method.tau = "DL"),
                 "argument 'method.tau' set to")
  expect_equal(ma5$method.tau, "REML")

  ma6 <- metaprop(events, N, data = dat, method = "Inverse",
                  method.tau = "REML")
  ma8 <- metaprop(events, N, data = dat, method = "Inverse",
                  method.tau = "REML", cluster = grp)
  ma7 <- update(ma6, cluster = grp)
  expect_equal(ma7$method.predict, ma8$method.predict)
  expect_equal(ma7$method.predict, "V")
  expect_equal(ma7$method.tau, ma8$method.tau)

  ma9 <- update(ma1, method = "Inverse", cluster = grp)
  ma10 <- metaprop(events, N, data = dat, method = "Inverse", cluster = grp)
  expect_equal(ma9$method.tau, ma10$method.tau)
})
