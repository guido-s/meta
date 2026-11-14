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
  expect_identical(metabind(list(ma, ma, ma))$overall$name,
                   c("meta1", "meta2", "meta3"))
})

test_that("generated subgroup names retain labels and avoid existing suffixes", {
  ma <- metagen(c(-0.4, 0.8, 0.5, -0.3, 0.9, 0.2),
                 c(0.1, 0.3, 0.2, 0.2, 0.4, 0.15),
                 subgroup = rep(c("A", "B"), each = 3),
                 subgroup.name = "Age group", method.tau = "DL")
  suffixed <- update(ma, subgroup.name = "Age group.1")
  expected <- c("Age group", "Age group.2", "Age group.1", "Age group.3")

  expect_identical(metabind(ma, ma, suffixed, ma)$overall$name, expected)
  expect_identical(metabind(list(ma, ma, suffixed, ma))$overall$name,
                   expected)
  expect_identical(ma$subgroup.name, "Age group")
  expect_identical(suffixed$subgroup.name, "Age group.1")
})

test_that("explicit names are neither overwritten nor made unique", {
  ma <- metagen(c(-0.4, 0.8, 0.5, -0.3, 0.9, 0.2),
                 c(0.1, 0.3, 0.2, 0.2, 0.4, 0.15),
                 subgroup = rep(c("A", "B"), each = 3),
                 subgroup.name = "Age group", method.tau = "DL")
  other <- update(ma, method.tau = "REML", method.random.ci = "HK",
                    tau.preset = 0.5, tau.common = TRUE)
  labels <- c("metagen", "meta")

  expect_identical(metabind(ma, other, name = labels)$overall$name, labels)
  expect_identical(metabind(list(ma, other), name = labels)$overall$name,
                   labels)
  expect_error(metabind(ma, other, name = rep("Age group", 2)),
                 "must be unique")
  expect_error(metabind(list(ma, other), name = rep("meta", 2)),
                 "must be unique")
})
