test_that("constructors without data do not use caller environments as masks", {
  old <- options(lifecycle_verbosity = "error")
  on.exit(options(old))

  estimates <- 1:3
  errors <- 3:1 / 10
  calls <- list(
    function() metagen(estimates, errors),
    function() metabin(1:3, 101:103, 4:6, 201:203),
    function() metacont(101:103, 1:3, 3:1, 101:103, 2:4, 3:1),
    function() metacor(c(0.1, 0.2, 0.3), 101:103),
    function() metainc(1:3, 101:103, 4:6, 201:203),
    function() metamean(101:103, 1:3, 3:1),
    function() metaprop(1:3, 101:103, method = "Inverse"),
    function() metarate(1:3, 101:103)
  )

  for (call in calls)
    expect_s3_class(call(), "meta")

  ma <- metagen(estimates, errors)
  expect_equal(ma$TE, estimates)
  expect_equal(ma$seTE, errors)
})

test_that("data columns take precedence and missing labels retain their defaults", {
  old <- options(lifecycle_verbosity = "error")
  on.exit(options(old))

  TE <- 99
  seTE <- 99
  dat <- data.frame(TE = 1:3, seTE = 3:1 / 10)
  wrapper <- function(...) metagen(...)
  ma <- wrapper(TE = TE, seTE = seTE, data = dat, sm = "MD")
  expect_equal(ma$TE, dat$TE)
  expect_equal(ma$seTE, dat$seTE)
  expect_length(ma$studlab, 3)
})

test_that("metabind updates do not generate environment-mask deprecations", {
  old <- options(lifecycle_verbosity = "error")
  on.exit(options(old))

  ma <- metabin(1:3, 101:103, 4:6, 201:203, sm = "OR", method = "Inverse")
  other <- update(ma, method = "GLMM")
  mb <- metabind(ma, update(other, method.predict = c("S", "V")))
  expect_s3_class(mb, "metabind")
})
