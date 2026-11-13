test_that("line-type lists preserve renderer interfaces and legacy arguments", {
  ma <- metagen(1:3, 3:1)
  observed <- NULL
  draw <- meta:::forest_meta_internal
  testthat::local_mocked_bindings(
    forest_meta_internal = function(...) {
      observed <<- list(...)
      do.call(draw, observed)
    }, .package = "meta")
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())

  for (legacy in c(FALSE, TRUE)) {
    if (legacy)
      forest(ma, lty.common = 2, lty.random = 3, lty.cid = 4)
    else
      forest(ma, lty = list(common = 2, random = 3, cid = 4))
    expect_equal(observed$lty.common, 2)
    expect_equal(observed$lty.random, 3)
    expect_equal(observed$lty.cid, 4)
  }
  forest(ma, lty = list(common = NULL, random = NULL))
  expect_null(observed$lty.common)
  expect_null(observed$lty.random)
  expect_warning(forest(ma, lty.equi = 5), "Use argument 'lty.cid'")
  expect_equal(observed$lty.cid, 5)
  expect_warning(forest(ma, lty = list(cid = 4), lty.equi = 5), "ignored")
  expect_equal(observed$lty.cid, 4)
  expect_error(forest(ma, lty = list(common = -1)), "lty[.]common")
})
