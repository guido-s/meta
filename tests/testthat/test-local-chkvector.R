test_that("chkvector accepts atomic vectors and returns invisibly", {
  for (value in list(numeric(), 1:3, c(0.1, NA_real_), c("a", "b"),
                     c(TRUE, FALSE), factor(c("b", "a")),
                     as.Date("2026-01-01") + 0:2))
    expect_identical(withVisible(meta:::chkvector(value)),
                     list(value = NULL, visible = FALSE))
})

test_that("chkvector only accepts NULL when requested", {
  expect_error(meta:::chkvector(NULL), "An object of class 'NULL' was supplied.",
               fixed = TRUE)
  expect_identical(withVisible(meta:::chkvector(NULL, NULL.ok = TRUE)),
                   list(value = NULL, visible = FALSE))
  expect_error(meta:::chkvector(list(), NULL.ok = TRUE),
               "An object of class 'list' was supplied.", fixed = TRUE)
})

test_that("chkvector reports argument names and rejected classes without a call", {
  invalid <- list(data.frame(rank = 1:3), matrix(1:3), array(1:3, c(3, 1, 1)),
                  list(1, 2, 3), sqrt, new.env())
  for (value in invalid) {
    error <- tryCatch(meta:::chkvector(value), error = identity)
    expect_s3_class(error, "error")
    expect_identical(conditionMessage(error),
                     paste0("Argument 'value' must be a vector ",
                            "(for example, numbers, text, or a factor).\n",
                            "An object of class '",
                            paste(class(value), collapse = "/"),
                            "' was supplied."))
    expect_null(conditionCall(error))
  }
  expect_error(meta:::chkvector(list(), name = "sortvar"),
               "Argument 'sortvar' must be a vector", fixed = TRUE)
})
