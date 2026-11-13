test_that("list2mat expands shared estimates using each source's own CI width", {
  values <- list(Left = c(A = 11, B = 12), Right = c(A = 21, B = 22))
  limits <- list(Left = c(A = 1, B = 2),
                  Right = matrix(1:4, nrow = 2,
                                 dimnames = list(c("A", "B"), c("classic", "HK"))))
  expected <- matrix(c(11, 12, 21, 22, 21, 22), nrow = 2,
                      dimnames = list(c("A", "B"), c("Left", "Right", "Right-HK")))
  expect_identical(meta:::list2mat(values, limits), expected)
  expect_identical(meta:::list2mat(rev(values), rev(limits)), expected[, c(2, 3, 1)])
  expect_identical(colnames(meta:::list2mat(values, limits)),
                   colnames(meta:::list2mat(limits)))
  expect_identical(colnames(meta:::list2mat(rev(values), rev(limits))),
                   colnames(meta:::list2mat(rev(limits))))
})

test_that("list2mat preserves method-specific matrices and their column order", {
  values <- list(
    Left = matrix(c(11, 12, 13, 14), nrow = 2,
                   dimnames = list(c("A", "B"), c("classic", "HK"))),
    Right = matrix(c(21, 22, 23, 24), nrow = 2,
                    dimnames = list(c("A", "B"), c("HK", "classic"))))
  expected <- matrix(c(11, 12, 13, 14, 21, 22, 23, 24), nrow = 2,
                      dimnames = list(c("A", "B"),
                                      c("Left", "Left-HK", "Right-HK", "Right")))
  expect_identical(meta:::list2mat(values), expected)
  expect_identical(meta:::list2mat(values, values), expected)
  expect_identical(meta:::list2mat(rev(values), rev(values)), expected[, c(3, 4, 1, 2)])

  mixed <- list(Left = values$Left, Right = c(A = 31, B = 32))
  expect_identical(meta:::list2mat(mixed),
                   matrix(c(11, 12, 13, 14, 31, 32), nrow = 2,
                          dimnames = list(c("A", "B"), c("Left", "Left-HK", "Right"))))
})

test_that("list2mat handles single subgroups, absent blocks, and native inputs", {
  reference <- matrix(c(1, 2), nrow = 1, dimnames = list("A", c("classic", "HK")))
  expect_identical(meta:::list2mat(list(Source = c(A = 10)), list(Source = reference)),
                   matrix(c(10, 10), nrow = 1,
                          dimnames = list("A", c("Source", "Source-HK"))))
  expect_identical(meta:::list2mat(list(Absent = NULL, Present = c(A = 1, B = 2))),
                   matrix(c(1, 2), ncol = 1,
                          dimnames = list(c("A", "B"), "Present")))
  expect_null(meta:::list2mat(list(NULL, NULL)))
  expect_identical(meta:::list2mat(c(A = 1, B = 2)), matrix(c(1, 2), ncol = 1))
  expect_identical(meta:::list2mat(reference), reference)
})

test_that("list2mat rejects incompatible reference blocks rather than recycling", {
  values <- list(Left = c(A = 1, B = 2), Right = c(A = 3, B = 4))
  expect_error(meta:::list2mat(values, matrix(1:4, nrow = 2)),
                "Reference blocks must correspond")
  expect_error(meta:::list2mat(values, list(1:2)), "Reference blocks must correspond")
  expect_error(meta:::list2mat(values, list(1:3, 1:2)), "different numbers of rows")
  expect_error(meta:::list2mat(list(matrix(1:4, nrow = 2)),
                              list(matrix(1:6, nrow = 2))), "incompatible dimensions")
})
