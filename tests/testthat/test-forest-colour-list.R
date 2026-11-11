test_that("colour lists preserve inheritance, vectors and legacy arguments", {
  ma <- metagen(seq(0.1, 0.6, 0.1), rep(0.1, 6),
                 subgroup = rep(c("A", "B"), each = 3), method.tau = "DL",
                 method.random.ci = c("classic", "HK"),
                 prediction = TRUE, method.predict = c("V", "HK"))
  colours <- c("red", "blue", "green", "orange", "purple", "cyan")
  settings <- list(common = "navy", random = "purple", diamond = "gold",
                    diamond.common = "red", diamond.lines = "brown",
                    square = colours, circle = colours,
                    study.text = rev(colours), inside = "white",
                    lines = "gray30", predict = c("cyan", "yellow"),
                    subgroup = "darkblue", main = "steelblue", label = "navy")
  observed <- NULL
  draw <- meta:::forest_meta_internal
  testthat::local_mocked_bindings(
    forest_meta_internal = function(...) {
      observed <<- list(...)
      do.call(draw, observed)
    }, .package = "meta")
  grDevices::pdf(NULL, width = 18, height = 14)
  withr::defer(grDevices::dev.off())

  for (legacy in c(FALSE, TRUE)) {
    if (legacy)
      options <- setNames(settings, paste0("col.", names(settings)))
    else
      options <- list(col = settings)
    do.call(meta::forest, c(list(x = ma, sortvar = -ma$TE,
                                 type = list(study = "circle")), options))
    expect_identical(observed$col.common, "navy")
    expect_identical(observed$col.random, "purple")
    expect_identical(observed$col.header.line, "gray30")
    expect_identical(observed$col.jama.line, "darkblue")
    expect_identical(observed$col.main, "steelblue")
    expect_equal(unname(observed$col.forest$col.square[1:5]),
                 c("red", "gold", "gold", "cyan", "yellow"))
    expect_equal(unname(observed$col.forest$col.square.lines[1:3]), rep("brown", 3))
    order <- c(3, 2, 1, 6, 5, 4)
    expect_equal(unname(utils::tail(observed$col.forest$col.square, 6)), colours[order])
    expect_equal(unname(utils::tail(observed$col.forest$col.circle.lines, 6)), colours[order])
    expect_length(observed$col.forest$col.circle.lines, length(observed$col.forest$col.square))
    labels <- utils::tail(observed$cols$col.studlab$labels, 6)
    expect_equal(unname(vapply(labels, function(grob) grob$gp$col, "")), rev(colours)[order])
  }
  expect_silent(forest(ma, col = list(common = NULL, diamond = NULL,
                                      square = NULL, circle = NULL, predict = NULL)))
  expect_error(forest(ma, col = list(square = c("red", "blue"))), "col[.]square")
})

test_that("named scalar font values retain their intended fallback keys", {
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
  expect_silent(forest(ma, fontsize = c(size = 12), fs = list(common = c(size = 14)),
                        ff = list(common = c(face = "italic"))))
  expect_equal(unname(observed$fs.main), 15)
  expect_equal(unname(observed$cols$col.effect$labels[[2]]$gp$fontsize), 14)
  expect_equal(unname(observed$cols$col.effect$labels[[3]]$gp$fontsize), 14)
  expect_equal(unname(observed$cols$col.effect$labels[[2]]$gp$font), 3)
  expect_equal(unname(observed$cols$col.effect$labels[[3]]$gp$font), 3)
})
