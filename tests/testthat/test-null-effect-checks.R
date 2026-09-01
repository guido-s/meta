test_that("metacor checks null effect on correlation scale", {
  cor <- c(-0.2, 0.1, 0.4)
  n <- c(20, 30, 40)

  expect_error(metacor(cor, n, sm = "ZCOR", null.effect = -1), NA)
  expect_error(metacor(cor, n, sm = "ZCOR", null.effect = 1), NA)
  expect_equal(metacor(cor, n, sm = "ZCOR", null.effect = -1)$null.effect,
               -1)
  expect_equal(metacor(cor, n, sm = "ZCOR", null.effect = 1)$null.effect,
               1)

  expect_error(metacor(cor, n, null.effect = -1.1),
               "Argument 'null.effect'")
  expect_error(metacor(cor, n, null.effect = 1.1),
               "Argument 'null.effect'")
})

test_that("metarate checks null effect on rate scale", {
  event <- c(2, 5, 8)
  time <- c(30, 40, 50)

  expect_error(metarate(event, time, sm = "IR", null.effect = 0), NA)
  expect_error(metarate(event, time, sm = "IRLN", null.effect = 0), NA)
  expect_equal(metarate(event, time, sm = "IR", null.effect = 0)$null.effect,
               0)
  expect_equal(metarate(event, time, sm = "IRLN", null.effect = 0)$null.effect,
               0)

  expect_error(metarate(event, time, null.effect = -0.1),
               "Argument 'null.effect'")
})

test_that("metamean checks null effect for log means", {
  n <- c(20, 30, 40)
  mean <- c(1, 2, 3)
  sd <- c(0.5, 0.6, 0.7)

  expect_error(metamean(n, mean, sd, sm = "MRAW", null.effect = -1), NA)
  expect_equal(metamean(n, mean, sd, sm = "MRAW", null.effect = -1)$null.effect,
               -1)

  expect_error(metamean(n, mean, sd, sm = "MLN", null.effect = 1), NA)
  expect_equal(metamean(n, mean, sd, sm = "MLN", null.effect = 1)$null.effect,
               1)

  expect_error(metamean(n, mean, sd, sm = "MLN", null.effect = 0),
               "Argument 'null.effect'")
  expect_error(metamean(n, mean, sd, sm = "MLN", null.effect = -1),
               "Argument 'null.effect'")
})

test_that("metagen checks null effect for single-arm summary measures", {
  TE <- c(0.1, 0.2, 0.3)
  seTE <- c(0.01, 0.02, 0.03)

  expect_error(metagen(TE, seTE, sm = "PRAW", null.effect = 0,
                       transf = FALSE), NA)
  expect_error(metagen(TE, seTE, sm = "PRAW", null.effect = 1,
                       transf = FALSE), NA)
  expect_error(metagen(TE, seTE, sm = "PLOGIT", null.effect = -0.1,
                       transf = FALSE),
               "Argument 'null.effect'")
  expect_error(metagen(TE, seTE, sm = "PLOGIT", null.effect = 1.1,
                       transf = FALSE),
               "Argument 'null.effect'")

  expect_error(metagen(TE, seTE, sm = "IR", null.effect = 0,
                       transf = FALSE), NA)
  expect_error(metagen(TE, seTE, sm = "IRLN", null.effect = 0,
                       transf = FALSE), NA)
  expect_error(metagen(TE, seTE, sm = "IR", null.effect = -0.1,
                       transf = FALSE),
               "Argument 'null.effect'")

  expect_error(metagen(TE, seTE, sm = "MRAW", null.effect = -1,
                       transf = FALSE), NA)
  expect_error(metagen(TE, seTE, sm = "MLN", null.effect = 1,
                       transf = FALSE), NA)
  expect_error(metagen(TE, seTE, sm = "MLN", null.effect = 0,
                       transf = FALSE),
               "Argument 'null.effect'")
  expect_error(metagen(TE, seTE, sm = "MLN", null.effect = -1,
                       transf = FALSE),
               "Argument 'null.effect'")

  expect_error(metagen(TE, seTE, sm = "COR", null.effect = -1,
                       transf = FALSE), NA)
  expect_error(metagen(TE, seTE, sm = "ZCOR", null.effect = 1,
                       transf = FALSE), NA)
  expect_error(metagen(TE, seTE, sm = "COR", null.effect = -1.1,
                       transf = FALSE),
               "Argument 'null.effect'")
  expect_error(metagen(TE, seTE, sm = "ZCOR", null.effect = 1.1,
                       transf = FALSE),
               "Argument 'null.effect'")
})
