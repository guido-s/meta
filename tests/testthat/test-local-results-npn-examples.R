test_that("metacont gives no weight to non-positive or missing sample sizes", {
  n.e <- c(10, 0, NA, 20)
  n.c <- c(10, 10, 10, 20)
  mean.e <- c(2, 3, 5, 6)
  mean.c <- c(1, 1, 1, 2)
  sd.e <- rep(1, 4)
  sd.c <- rep(1, 4)

  expect_warning(
    metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c),
    "non-positive values for n.e and / or n.c"
  )

  m <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, warn = FALSE)

  expect_equal(m$TE[c(1, 4)], mean.e[c(1, 4)] - mean.c[c(1, 4)])
  expect_true(all(is.na(m$TE[2:3])))
  expect_true(all(is.na(m$seTE[2:3])))
})

test_that("metamean gives no weight to non-positive or missing sample sizes", {
  n <- c(10, 0, NA, 20)
  mean <- c(2, 3, 5, 6)
  sd <- rep(1, 4)

  expect_warning(
    metamean(n, mean, sd),
    "non-positive sample size"
  )

  m <- metamean(n, mean, sd, warn = FALSE)

  expect_equal(m$TE[c(1, 4)], mean[c(1, 4)])
  expect_true(all(is.na(m$TE[2:3])))
  expect_true(all(is.na(m$seTE[2:3])))
})

test_that("metamean log means give no weight to non-positive or missing means", {
  n <- rep(10, 4)
  mean <- c(2, 0, NA, 6)
  sd <- rep(1, 4)

  expect_warning(
    metamean(n, mean, sd, sm = "MLN", null.effect = 1),
    "negative or zero mean"
  )

  m <- metamean(n, mean, sd, sm = "MLN", null.effect = 1, warn = FALSE)

  expect_equal(m$TE[c(1, 4)], log(mean[c(1, 4)]))
  expect_true(all(is.na(m$TE[2:3])))
  expect_true(all(is.na(m$seTE[2:3])))
})
