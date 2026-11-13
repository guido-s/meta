test_that("condense defaults to legacy flattening and name handling", {
  x <- list(First = list(result = c(left = 1, right = 2)),
            Null = list(result = NULL), Missing = list(other = 9),
            Nested = list(result = list(inner = c(a = 3, b = 4))))
  expected <- unlist(lapply(x, "[[", "result"))

  expect_identical(meta:::condense(x, "result"), expected)
  expect_identical(meta:::condense(x, "result", value = NULL), expected)
  expect_identical(meta:::condense(unname(x), "result"),
                   unlist(lapply(unname(x), "[[", "result")))
  expect_null(meta:::condense(list(list(result = NULL), list()), "result"))
  expect_null(meta:::condense(list(), "result"))
})

test_that("typed condense returns unnamed scalar vectors and typed empties", {
  cases <- list(character = c("Age group", "Region"),
                logical = c(TRUE, FALSE), numeric = c(1.5, 2.5))
  for (values in cases) {
    x <- list(First = list(result = c(inner = values[1])),
              Second = list(result = c(inner = values[2])))
    prototype <- values[1]

    expect_identical(meta:::condense(x, "result", prototype), values)
    expect_identical(meta:::condense(x, value = prototype, var = "result"),
                     values)
    expect_identical(meta:::condense(list(), "result", prototype), values[FALSE])
  }
})

test_that("typed condense deliberately retains vapply validation", {
  invalid <- list(list(result = NULL), list(), list(result = character()),
                  list(result = c("a", "b")))
  for (element in invalid)
    expect_error(meta:::condense(list(element), "result", character(1)),
                   "values must be length 1")

  expect_error(meta:::condense(list(list(result = 1)), "result", character(1)),
                 "values must be type 'character'")
  expect_error(meta:::condense(list(list(result = "TRUE")), "result", logical(1)),
                 "values must be type 'logical'")
  expect_error(meta:::condense(list(list(result = "1")), "result", numeric(1)),
                 "values must be type 'double'")
})
