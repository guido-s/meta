test_that("resolved font-size lists retain formatting and renderer interfaces", {
  ma <- metagen(seq(0.1, 0.6, 0.1), rep(0.1, 6),
                 subgroup = rep(c("A", "B"), each = 3),
                 prediction = TRUE, method.predict = "V", method.tau = "DL")
  sizes <- list(main = 16, heading = 13, common = 9, random = 10,
                 predict = 11, study = 7, study.labels = 8,
                 axis = 6, xlab = 12, lr = 14, rob.symbols = 15)
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
      options <- setNames(sizes, paste0("fs.", names(sizes)))
    else
      options <- list(fs = sizes)
    do.call(meta::forest, c(list(x = ma, main = "Font sizes", details = TRUE), options))

    expected <- list(fs.main = 16, fs.head = 13, fs.axis = 6,
                      fs.xlab = 12, fs.lr = 14, fs.rob.symbols = 15)
    for (i in names(expected))
      expect_equal(observed[[i]], expected[[i]], info = i)

    effect.sizes <- vapply(observed$cols$col.effect$labels,
                           function(grob) grob$gp$fontsize, numeric(1))
    expect_true(all(c(7, 9, 10, 11, 13) %in% effect.sizes))
    label.sizes <- vapply(observed$cols$col.studlab$labels,
                          function(grob) grob$gp$fontsize, numeric(1))
    expect_true(8 %in% label.sizes)
  }
})
