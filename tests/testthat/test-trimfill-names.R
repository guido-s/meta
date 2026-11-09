test_that("trimfill results have unique names and can be printed", {
  ma <- metagen(1:5, 1:5 / 10, level = 0.9)
  result <- trimfill(ma)

  expect_false(anyDuplicated(names(result)) > 0)
  expect_equal(result$level, ma$level)
  expect_length(result[names(result) == "level"], 1)
  expect_output(print(result), "Details")
})
