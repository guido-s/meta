test_that("update.meta resets method.tau.ci when adding cluster", {
  dat <- data.frame(events = 1:5, N = 5:1 * 100, grp = c(1, 1, 2, 2, 3))

  ma1 <- metaprop(events, N, data = dat)
  ma2 <- metaprop(events, N, data = dat, cluster = grp)
  ma3 <- update(ma1, cluster = grp)

  expect_false(ma1$three.level)
  expect_true(ma2$three.level)
  expect_true(ma3$three.level)
  expect_equal(ma3$method.tau.ci, ma2$method.tau.ci)
  expect_equal(ma3$method.tau.ci, "PL")

  ma4 <- update(ma1, cluster = grp, method.tau.ci = "")
  expect_equal(ma4$method.tau.ci, "")
})
