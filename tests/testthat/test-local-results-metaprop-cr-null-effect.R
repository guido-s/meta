test_that("cluster-robust proportion tests require a null effect", {
  data(Fleiss1993bin, package = "meta")
  ma <- metaprop(d.asp, n.asp, data = Fleiss1993bin, method = "Inverse",
                 method.tau = "REML", method.tau.ci = "QP")

  for (method in c("CR0", "CR1", "CR2")) {
    result <- update(ma, method.random.ci = method)
    expect_true(is.na(result$null.effect))
    expect_true(is.na(result$statistic.random))
    expect_true(is.na(result$pval.random))
    expect_true(is.finite(result$lower.random))
    expect_true(is.finite(result$upper.random))

    with.null <- update(result, null.effect = 0.1)
    expect_true(is.finite(with.null$statistic.random))
    expect_true(is.finite(with.null$pval.random))
    expect_equal(result$lower.random, with.null$lower.random)
    expect_equal(result$upper.random, with.null$upper.random)
  }
})
