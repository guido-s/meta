test_that("chksetVar resolves defaults and recycles only scalars", {
  colour <- "red"
  expect_identical(meta:::chksetVar(colour, FALSE, "black", meta:::chkcolor, n = 3),
                   rep("red", 3))
  expect_identical(meta:::chksetVar(NULL, TRUE, "blue", meta:::chkcolor, n = 3),
                   rep("blue", 3))
  expect_identical(meta:::chksetVar(NULL, TRUE, c("red", "blue"),
                                    meta:::chkcolor, n = 2), c("red", "blue"))
  expect_identical(meta:::chksetVar(stop("input forced"), TRUE, "red",
                                    meta:::chkcolor), "red")
  expect_identical(meta:::chksetVar("red", FALSE, stop("default forced"),
                                    meta:::chkcolor), "red")
  expect_identical(meta:::chksetVar(NULL, TRUE, "red", meta:::chkcolor,
                                    oldname = stop("oldname forced"),
                                    warn = stop("warn forced")), "red")
})

test_that("chksetVar preserves names, checker options and legacy values", {
  colour <- list("invalid")
  expect_error(meta:::chksetVar(colour, FALSE, "black", meta:::chkcolor),
               "Argument 'colour'")
  expect_error(meta:::chksetVar(NULL, TRUE, FALSE, meta:::chkcolor),
               "character or numeric")
  expect_identical(meta:::chksetVar(2, FALSE, 0, meta:::chknumeric_strict,
                                    n = 2, min = 0), c(2, 2))
  expect_identical(meta:::chksetVar(NULL, TRUE, "black", meta:::chkcolor,
                                    args = list(old.colour = "blue"),
                                    oldname = "old.colour", warn = FALSE), "blue")
  colour <- "red"
  expect_warning(result <- meta:::chksetVar(colour, FALSE, "black", meta:::chkcolor,
                                             args = list(old.colour = "blue"),
                                             oldname = "old.colour", warn = TRUE),
                 "'colour' is also provided")
  expect_identical(result, "red")
})

test_that("chksetVar_list preserves settings and validates named diagnostics", {
  colours <- list(first = NULL, second = "green", untouched = "white")
  result <- meta:::chksetVar_list(
    colours, list(first = TRUE, second = FALSE),
    list(first = c("red", "blue"), second = "black"),
    meta:::chkcolor, n = list(first = 2, second = 3))
  expect_identical(result$first, c("red", "blue"))
  expect_identical(result$second, rep("green", 3))
  expect_identical(result$untouched, "white")
  expect_error(meta:::chksetVar_list(colours, list(first = TRUE), list(first = TRUE),
                                      meta:::chkcolor), "colours[.]first")
  expect_error(meta:::chksetVar_list(colours, list(first = TRUE), list(first = TRUE),
                                      meta:::chkcolor, name = list(first = "custom")),
               "Argument 'custom'")
  expect_identical(meta:::chksetVar_list(
    colours, list(), list(second = "black"), meta:::chkcolor)$second, "green")
  result <- meta:::chksetVar_list(
    colours, list(first = TRUE), list(first = "black"), meta:::chkcolor,
    args = list(old.first = "blue"), oldname = list(first = "old.first"), warn = FALSE)
  expect_identical(result$first, "blue")
})
