test_that("setunit_list converts strings and preserves other list values", {
  values <- list(mm = "3mm", cm = "2cm", inch = "1inch",
                  existing = grid::unit(4, "mm"), number = 2, missing = NULL)
  result <- meta:::setunit_list(values)
  expect_named(result, names(values))
  expect_equal(as.numeric(result$mm), 3)
  expect_equal(grid::unitType(result$mm), "mm")
  expect_equal(as.numeric(result$cm), 2)
  expect_equal(grid::unitType(result$cm), "cm")
  expect_equal(grid::unitType(result$inch), "inches")
  expect_identical(result$existing, values$existing)
  expect_identical(result$number, values$number)
  expect_null(result$missing)
})

test_that("colgap lists retain inheritance and legacy result interfaces", {
  ma <- metagen(1:3, 3:1)
  gaps <- list(colgap = "3mm", left = "5mm", forest = "8mm",
                forest.left = "9mm", forest.right = "10mm")
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  for (legacy in c(FALSE, TRUE)) {
    if (legacy) {
      options <- gaps
      names(options) <- c("colgap", "colgap.left", "colgap.forest",
                           "colgap.forest.left", "colgap.forest.right")
    }
    else
      options <- list(colgap = gaps)
    result <- do.call(meta::forest, c(list(x = ma), options))
    expected <- c(colgap = 3, colgap.left = 5, colgap.right = 3,
                   colgap.studlab = 5, colgap.forest = 8,
                   colgap.forest.left = 9, colgap.forest.right = 10)
    for (i in names(expected)) {
      expect_true(grid::is.unit(result[[i]]), info = i)
      expect_equal(grid::convertWidth(result[[i]], "mm", valueOnly = TRUE),
                   unname(expected[[i]]), info = i)
    }
  }
})
