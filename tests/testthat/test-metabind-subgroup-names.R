test_that("subgroup names provide default component names", {
  ma <- metagen(c(-0.4, 0.8, 0.5, -0.3, 0.9, 0.2),
                 c(0.1, 0.3, 0.2, 0.2, 0.4, 0.15),
                 subgroup = rep(c("A", "B"), each = 3),
                 subgroup.name = "Age group", method.tau = "DL")
  other <- update(ma, subgroup = rep(c("C", "D"), 3),
                   subgroup.name = "Region")
  expected <- c("Age group", "Region")
  expect_identical(metabind(ma, other)$overall$name, expected)
  expect_identical(metabind(list(ma, other))$overall$name, expected)
  expect_identical(metabind(ma, other, name = c("First", "Second"))$overall$name,
                   c("First", "Second"))

  ma <- metagen(c(0.2, 0.5, 0.8), rep(0.2, 3), method.tau = "DL")
  expect_identical(metabind(ma, ma)$overall$name, c("meta1", "meta2"))
})
