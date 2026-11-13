test_that("font-face lists preserve mixed numeric and character faces", {
  ma <- metagen(seq(0.1, 0.6, 0.1), rep(0.1, 6),
                 subgroup = rep(c("A", "B"), each = 3),
                 prediction = TRUE, method.predict = "V", method.tau = "DL")
  faces <- list(main = 3, heading = "bold.italic", common = 2,
                 random = "italic", predict = 4, hetstat = 3,
                 axis = 1, xlab = "bold", lr = 4, rob.symbols = 2)
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
      options <- setNames(faces, paste0("ff.", names(faces)))
    else
      options <- list(ff = faces)
    do.call(meta::forest, c(list(x = ma, main = "Font faces", details = TRUE), options))
    expected <- list(ff.main = 3, ff.head = "bold.italic", ff.axis = 1,
                      ff.xlab = "bold", ff.lr = 4, ff.rob.symbols = 2)
    for (i in names(expected))
      expect_identical(observed[[i]], expected[[i]], info = i)
    expect_equal(unname(observed$cols$col.effect$labels[[2]]$gp$font), 2)
    expect_equal(unname(observed$cols$col.effect$labels[[3]]$gp$font), 3)
    expect_equal(unname(observed$cols$col.effect$labels[[4]]$gp$font), 4)
  }

  for (layout in c("meta", "BMJ", "JAMA")) {
    forest(ma, layout = layout, ff = list(common = NULL, random = NULL,
                                          predict = NULL, hetstat = NULL, lr = NULL))
    expect_equal(unname(observed$cols$col.effect$labels[[2]]$gp$font),
                 if (layout == "meta") 2 else 1)
    expect_identical(observed$ff.lr, if (layout == "meta") "plain" else "bold")
  }
  expect_error(forest(ma, ff = list(common = "invalid")), "ff[.]common")
})
