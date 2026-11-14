test_that("setchar_list normalizes choices and supports named length constraints", {
  choices <- c("left", "center", "right")
  input <- list(general = c("l", "r"), title = "C", absent = NULL)
  result <- meta:::setchar_list(input, choices, length = list(title = 1))
  expect_named(result, names(input))
  expect_identical(result$general, c("left", "right"))
  expect_identical(result$title, "center")
  expect_null(result$absent)
  expect_error(meta:::setchar_list(input, choices, length = 1), "input[.]general")
  input$title <- c("left", "right")
  expect_error(meta:::setchar_list(input, choices, length = list(title = 1)), "input[.]title")
})

test_that("justification lists preserve vector settings and legacy interfaces", {
  dat <- data.frame(TE = 1:3, seTE = 3:1, first = 4:6, second = 7:9, third = 10:12)
  ma <- metagen(TE, seTE, data = dat)
  settings <- list(just = "c", main = "r", studlab = "l",
                    addcols = c("r", "l"), addcols.left = c("l", "r"),
                    addcols.right = "c", label.e = "c", label.c = "l", rob = "r")
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
    if (legacy) {
      options <- settings
      names(options) <- c("just", paste0("just.", names(settings)[-1]))
    }
    else
      options <- list(just = settings)
    result <- do.call(meta::forest,
                      c(list(x = ma, main = "Justification",
                             leftcols = c("studlab", "first", "second", "effect", "ci"),
                             rightcols = "third"), options))
    expect_identical(observed$just.main, "right")
    expect_identical(observed$just.label.e, "center")
    expect_identical(observed$just.label.c, "left")
    expect_identical(observed$just.addcols.left, c("left", "right"))
    expect_identical(observed$just.addcols.right, "center")
    expect_equal(result$leftjust, c("left", "left", "right", "center", "center"))
    expect_equal(result$rightjust, "center")
  }
  forest(ma, layout = "BMJ", just = list(main = NULL))
  expect_identical(observed$just.main, "center")
  forest(ma, layout = "meta", just = list(main = NULL))
  expect_identical(observed$just.main, "left")
  expect_error(forest(ma, just = list(main = c("left", "right"))), "just[.]main")
  expect_error(forest(ma, just = list(studlab = "invalid")), "just[.]studlab")
})
