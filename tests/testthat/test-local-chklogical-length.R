test_that("chklogical length zero disables only the length check", {
  expect_silent(meta:::chklogical(TRUE, length = 0))
  expect_silent(meta:::chklogical(c(TRUE, FALSE, TRUE), length = 0))
  expect_silent(meta:::chklogical(logical(), length = 0))
  expect_silent(meta:::chklogical(matrix(c(TRUE, FALSE), 1), length = 0))
  expect_silent(meta:::chklogical(c(0, 1), length = 0))
  expect_error(meta:::chklogical(c(TRUE, NA), length = 0), "logical vector")
  expect_error(meta:::chklogical(c("TRUE", "FALSE"), length = 0), "logical vector")
  expect_error(meta:::chklogical(NULL, length = 0), "logical vector")
  expect_silent(meta:::chklogical(NULL, length = 0, NULL.ok = TRUE))
})

test_that("chklogical retains scalar and multiple allowed length checks", {
  expect_error(meta:::chklogical(c(TRUE, FALSE)), "logical")
  expect_error(meta:::chklogical(logical()), "logical")
  expect_silent(meta:::chklogical(c(TRUE, FALSE), length = c(1, 2)))
  expect_error(meta:::chklogical(rep(TRUE, 3), length = c(1, 2)), "whose length is")
  expect_silent(meta:::chklogical(logical(), length = c(0, 1)))
  expect_error(meta:::chklogical(c(TRUE, FALSE), length = c(0, 1)), "whose length is")
})

test_that("forest subgroup selection accepts one value or one per subgroup", {
  ma <- metagen(seq(0.1, 0.6, 0.1), rep(0.1, 6),
                 subgroup = rep(c("A", "B"), each = 3), method.tau = "DL")
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  expect_silent(forest(ma, subgroup = TRUE))
  expect_silent(forest(ma, subgroup = c(TRUE, FALSE)))
  expect_error(forest(ma, subgroup = rep(TRUE, 6)), "whose length is 1 or 2")
  expect_error(forest(ma, subgroup = c(TRUE, NA)), "logical vector")
})
