test_that("chknumeric vector mode accepts numeric vectors and preserves numeric checks", {
  for (value in list(1:3, c(.1, .2), c(1, NA_real_)))
    expect_identical(withVisible(meta:::chknumeric(value, vector = TRUE)),
                     list(value = NULL, visible = FALSE))
  expect_error(meta:::chknumeric("1", vector = TRUE), "Non-numeric value")
  expect_error(meta:::chknumeric(c(1, NA), vector = TRUE, NA.ok = FALSE),
               "Missing values not allowed")
  expect_error(meta:::chknumeric(-1, min = 0, vector = TRUE),
               "must be larger equal 0")
  expect_error(meta:::chknumeric(1.5, integer = TRUE, vector = TRUE),
               "must be an integer")
  expect_error(meta:::chknumeric(1:3, length = 2, vector = TRUE),
               "must be a numeric of length 2")
})

test_that("chknumeric checks original vector shape before filtering missing values", {
  for (value in list(matrix(1:4, 2), matrix(c(1, NA, 2, NA), 2),
                     matrix(NA_real_, 2, 2), list(1, 2), list(NA, NA), sqrt))
    expect_error(meta:::chknumeric(value, vector = TRUE),
                 "Argument 'value' must be a vector", fixed = TRUE)
  expect_error(meta:::chknumeric(matrix(1:4, 2), vector = TRUE, name = "TE"),
               "Argument 'TE' must be a vector", fixed = TRUE)
})

test_that("chknumeric vector mode respects NULL.ok", {
  expect_null(meta:::chknumeric(NULL, vector = TRUE))
  expect_null(meta:::chknumeric(numeric(), vector = TRUE))
  expect_null(meta:::chknumeric(NULL, vector = TRUE, NULL.ok = TRUE))
  expect_error(meta:::chknumeric(NULL, vector = TRUE, NULL.ok = FALSE),
               "An object of class 'NULL' was supplied.", fixed = TRUE)
  expect_error(meta:::chknumeric(numeric(), vector = TRUE, NULL.ok = FALSE),
               "must not have length 0")
})

test_that("chknumeric preserves matrix acceptance unless vector mode is requested", {
  for (value in list(matrix(1:4, 2), matrix(c(1, NA, 2, NA), 2))) {
    expect_silent(meta:::chknumeric(value))
    expect_silent(meta:::chknumeric(value, vector = FALSE))
  }
})
