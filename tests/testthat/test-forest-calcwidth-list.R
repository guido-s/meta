test_that("calcwidth lists retain pooled inheritance and legacy controls", {
  ma <- metagen(seq(0.1, 0.6, 0.1), rep(0.1, 6),
                 subgroup = rep(c("A", "B"), each = 3), method.tau = "DL")
  observed <- logical()
  calculate <- meta:::calc_pooled_prefix_extra
  testthat::local_mocked_bindings(
    calc_pooled_prefix_extra = function(lines, pooled.cols, do.calc, ...) {
      observed <<- c(observed, do.calc)
      calculate(lines, pooled.cols, do.calc, ...)
    }, .package = "meta")
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())

  settings <- list(list(pooled = FALSE),
                    list(pooled = FALSE, random = TRUE),
                    list(pooled = TRUE, common = FALSE))
  expected <- list(c(FALSE, FALSE), c(FALSE, TRUE), c(FALSE, TRUE))
  for (legacy in c(FALSE, TRUE)) {
    for (i in seq_along(settings)) {
      if (legacy)
        options <- setNames(settings[[i]], paste0("calcwidth.", names(settings[[i]])))
      else
        options <- list(calcwidth = settings[[i]])
      observed <- logical()
      do.call(meta::forest, c(list(x = ma), options))
      expect_identical(observed, expected[[i]])
    }
  }
  observed <- logical()
  forest(ma, common = FALSE, calcwidth = list(pooled = TRUE))
  expect_identical(observed, c(FALSE, TRUE))
  expect_error(forest(ma, calcwidth = list(common = "invalid")), "calcwidth[.]common")
})
