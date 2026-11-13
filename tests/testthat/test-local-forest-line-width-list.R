test_that("line-width lists preserve drawing interfaces and legacy arguments", {
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
      forest(ma, lwd = 2, lwd.square = 3, lwd.diamond = 4)
    else
      forest(ma, lwd = list(lwd = 2, square = 3, diamond = 4))
    expect_equal(observed$lwd, 2)
    expect_equal(observed$col.forest$lwd, 2)
    expect_equal(observed$col.forest$lwd.square, 3)
    expect_equal(observed$col.forest$lwd.diamond, 4)
  }
  expect_error(forest(ma, lwd = list(square = 1:2)), "lwd[.]square")
})
